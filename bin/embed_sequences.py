#!/usr/bin/env python3
"""
Compute per-protein embeddings and save to embeddings.npz.

--model none     : amino-acid composition (mean-pooled one-hot), 21-dimensional
--model esm2     : ESM2 650M, mean-pooled over residues (dim 1280)
--model prot_t5  : ProtT5-XL, mean-pooled over residues (dim 1024)

Each entry in the NPZ file is keyed by protein ID and holds a 1-D float32 array.

--require-gpu (passed under -profile gpu) makes a missing CUDA device an error
instead of a fall back to CPU.
"""

import argparse
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import read_fasta

_AA_VOCAB = "ACDEFGHIKLMNPQRSTVWY"
_AA_IDX = {aa: i for i, aa in enumerate(_AA_VOCAB)}

# Both neural models truncate here, which also caps a batch's padded width.
_MAX_LEN = 1024
# Padded tokens and sequences per forward pass.
_BATCH_TOKENS = 16384
_MAX_BATCH = 64


def embed_one_hot(seqs):
    dim = len(_AA_VOCAB) + 1  # 20 standard AAs + 1 for unknown/non-standard
    embeddings = {}
    for pid, seq in seqs.items():
        vec = np.zeros(dim, dtype=np.float32)
        for aa in seq:
            vec[_AA_IDX.get(aa, len(_AA_VOCAB))] += 1
        if seq:
            vec /= len(seq)
        embeddings[pid] = vec
    return embeddings


def resolve_device(require_gpu):
    """Pick cuda when it is genuinely usable, else error or warn per require_gpu."""
    import torch

    built_for = torch.version.cuda or "cpu-only"
    if torch.cuda.is_available():
        name = torch.cuda.get_device_name(0)
        print(f"  device: cuda ({name}); torch {torch.__version__} built for CUDA {built_for}", file=sys.stderr)
        return torch.device("cuda")

    detail = f"torch {torch.__version__} (built for CUDA {built_for}) reports no usable CUDA device"
    if require_gpu:
        sys.exit(
            f"ERROR: {detail}, but this run was launched with -profile gpu.\n"
            "  The torch wheel likely targets a newer CUDA than the node driver (see .command.err).\n"
            "  To embed on CPU deliberately, drop -profile gpu."
        )
    print(f"WARNING: {detail}; embedding on CPU at roughly 2.5 s per sequence.", file=sys.stderr)
    return torch.device("cpu")


def batches(tokens, max_tokens, max_batch):
    """Group ids into batches under a padded-token budget, longest sequence first (ties by id)."""
    batch = []
    for pid in sorted(tokens, key=lambda p: (-tokens[p], p)):
        width = tokens[batch[0]] if batch else tokens[pid]
        if batch and ((len(batch) + 1) * width > max_tokens or len(batch) >= max_batch):
            yield batch
            batch = []
        batch.append(pid)
    if batch:
        yield batch


def residue_mask(attention_mask):
    """attention_mask with the CLS (first) and EOS (last real) positions cleared."""
    import torch

    mask = attention_mask.clone()
    mask[:, 0] = 0
    eos = attention_mask.sum(1) - 1
    mask[torch.arange(mask.size(0), device=mask.device), eos] = 0
    return mask


def mean_pool(hidden, mask):
    """Mean over the positions mask keeps; zero vector if it keeps none."""
    m = mask.unsqueeze(-1).to(hidden.dtype)
    return (hidden * m).sum(1) / m.sum(1).clamp(min=1)


def _report(done, total, reported, step=1000):
    if done - reported >= step or done == total:
        print(f"  {done}/{total} embedded", file=sys.stderr)
        return done
    return reported


def embed_esm2(seqs, device, max_tokens, max_batch):
    import torch
    from transformers import AutoModel, AutoTokenizer

    tokenizer = AutoTokenizer.from_pretrained("facebook/esm2_t33_650M_UR50D")
    model = AutoModel.from_pretrained("facebook/esm2_t33_650M_UR50D").to(device).eval()

    # +2 for the CLS and EOS the tokenizer wraps each sequence in.
    tokens = {pid: min(len(seq), _MAX_LEN - 2) + 2 for pid, seq in seqs.items()}
    embeddings = {}
    done = reported = 0
    for batch in batches(tokens, max_tokens, max_batch):
        enc = tokenizer(
            [seqs[pid] for pid in batch],
            return_tensors="pt",
            padding=True,
            truncation=True,
            max_length=_MAX_LEN,
        ).to(device)
        with torch.no_grad():
            hidden = model(**enc).last_hidden_state
        pooled = mean_pool(hidden, residue_mask(enc["attention_mask"]))
        for pid, vec in zip(batch, pooled):
            embeddings[pid] = vec.cpu().float().numpy()
        done += len(batch)
        reported = _report(done, len(seqs), reported)
    return embeddings


def embed_prot_t5(seqs, device, max_tokens, max_batch):
    import torch
    from transformers import T5EncoderModel, T5Tokenizer

    tokenizer = T5Tokenizer.from_pretrained("Rostlab/prot_t5_xl_half_uniref50-enc", do_lower_case=False)
    model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_half_uniref50-enc").to(device).eval()

    # ProtT5 expects space-separated AAs; replace non-standard residues with X
    fmt = {pid: " ".join(re.sub(r"[UZOB]", "X", seq)) for pid, seq in seqs.items()}
    # +1 for the trailing </s>; T5 prepends no CLS.
    tokens = {pid: min(len(seq), _MAX_LEN - 1) + 1 for pid, seq in seqs.items()}
    embeddings = {}
    done = reported = 0
    for batch in batches(tokens, max_tokens, max_batch):
        enc = tokenizer(
            [fmt[pid] for pid in batch],
            return_tensors="pt",
            padding=True,
            truncation=True,
            max_length=_MAX_LEN,
        ).to(device)
        with torch.no_grad():
            hidden = model(**enc).last_hidden_state
        # Pools every real position, the trailing </s> included.
        pooled = mean_pool(hidden, enc["attention_mask"])
        for pid, vec in zip(batch, pooled):
            embeddings[pid] = vec.cpu().float().numpy()
        done += len(batch)
        reported = _report(done, len(seqs), reported)
    return embeddings


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", nargs="+", required=True)
    ap.add_argument("--model", choices=["none", "esm2", "prot_t5"], default="esm2")
    ap.add_argument(
        "--require-gpu",
        action="store_true",
        help="abort instead of falling back to CPU when no CUDA device is usable",
    )
    args = ap.parse_args()

    seqs = {}
    for path in args.fasta:
        seqs.update(read_fasta(path))
    print(f"Embedding {len(seqs)} unique proteins with {args.model}", file=sys.stderr)

    if args.model == "none":
        embeddings = embed_one_hot(seqs)
    else:
        device = resolve_device(args.require_gpu)
        embed = embed_esm2 if args.model == "esm2" else embed_prot_t5
        embeddings = embed(seqs, device, _BATCH_TOKENS, _MAX_BATCH)

    # Sorted, so the archive is byte-reproducible.
    np.savez("embeddings.npz", **{pid: embeddings[pid] for pid in sorted(embeddings)})
    print(f"Saved {len(embeddings)} embeddings to embeddings.npz", file=sys.stderr)


if __name__ == "__main__":
    main()
