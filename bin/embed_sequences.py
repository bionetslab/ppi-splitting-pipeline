#!/usr/bin/env python3
"""
Compute per-protein embeddings and save to embeddings.npz.

--model none     : amino-acid composition (mean-pooled one-hot), 21-dimensional
--model esm2     : ESM2 650M, mean-pooled over residues (dim 1280)
--model prot_t5  : ProtT5-XL, mean-pooled over residues (dim 1024)

Each entry in the NPZ file is keyed by protein ID and holds a 1-D float32 array.
"""

import argparse
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import read_fasta

_AA_VOCAB = "ACDEFGHIKLMNPQRSTVWY"
_AA_IDX   = {aa: i for i, aa in enumerate(_AA_VOCAB)}


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


def embed_esm2(seqs):
    import torch
    from transformers import AutoModel, AutoTokenizer

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    tokenizer = AutoTokenizer.from_pretrained("facebook/esm2_t33_650M_UR50D")
    model = AutoModel.from_pretrained("facebook/esm2_t33_650M_UR50D").to(device).eval()

    embeddings = {}
    for i, (pid, seq) in enumerate(seqs.items()):
        inputs = tokenizer(seq, return_tensors="pt", truncation=True, max_length=1024).to(device)
        with torch.no_grad():
            out = model(**inputs)
        # skip CLS (0) and EOS (-1) tokens, mean-pool the rest
        embeddings[pid] = out.last_hidden_state[0, 1:-1].mean(0).cpu().float().numpy()
        if (i + 1) % 50 == 0:
            print(f"  {i + 1}/{len(seqs)} embedded", file=sys.stderr)
    return embeddings


def embed_prot_t5(seqs):
    import torch
    from transformers import T5EncoderModel, T5Tokenizer

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    tokenizer = T5Tokenizer.from_pretrained(
        "Rostlab/prot_t5_xl_half_uniref50-enc", do_lower_case=False
    )
    model = T5EncoderModel.from_pretrained(
        "Rostlab/prot_t5_xl_half_uniref50-enc"
    ).to(device).eval()

    embeddings = {}
    for i, (pid, seq) in enumerate(seqs.items()):
        # ProtT5 expects space-separated AAs; replace non-standard residues with X
        seq_fmt = " ".join(re.sub(r"[UZOB]", "X", seq))
        inputs = tokenizer(seq_fmt, return_tensors="pt", truncation=True, max_length=1024).to(device)
        with torch.no_grad():
            out = model(**inputs)
        embeddings[pid] = out.last_hidden_state[0].mean(0).cpu().float().numpy()
        if (i + 1) % 50 == 0:
            print(f"  {i + 1}/{len(seqs)} embedded", file=sys.stderr)
    return embeddings


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", nargs="+", required=True)
    ap.add_argument("--model", choices=["none", "esm2", "prot_t5"], default="esm2")
    args = ap.parse_args()

    seqs = {}
    for path in args.fasta:
        seqs.update(read_fasta(path))
    print(f"Embedding {len(seqs)} unique proteins with {args.model}", file=sys.stderr)

    if args.model == "none":
        embeddings = embed_one_hot(seqs)
    elif args.model == "esm2":
        embeddings = embed_esm2(seqs)
    else:
        embeddings = embed_prot_t5(seqs)

    np.savez("embeddings.npz", **embeddings)
    print(f"Saved {len(embeddings)} embeddings to embeddings.npz", file=sys.stderr)


if __name__ == "__main__":
    main()
