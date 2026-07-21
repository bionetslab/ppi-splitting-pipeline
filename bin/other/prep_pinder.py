"""Before any of this, download pinder with gsutil -m cp gs://pinder/2024-02/ your_pinder_dir/pinder"""

import json
import os
import re
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

from bin.fetch_data import (
    BATCH_SIZE,
    InvalidAccessionBatch,
    _ssl_context,
    build_canonical_map,
    fetch_batch,
    parse_batch,
    write_go_tsv,
    write_species_tsv,
)

PATH_TO_PINDER = "~/Downloads/pinder/2024-02"
RCSB_FASTA_URL = "https://www.rcsb.org/fasta/entry/{pdb_id}"
OUT_DIR = "../../data"
SEQ_OUT = f"{OUT_DIR}/pinder_sequences.fasta"
SPECIES_OUT = f"{OUT_DIR}/pinder_species.tsv"
GO_OUT = f"{OUT_DIR}/pinder_go_annotations.tsv"
CSV_OUT = f"{OUT_DIR}/pinder.csv"
RCSB_CACHE = f"{OUT_DIR}/pinder_rcsb_cache.jsonl"
RCSB_WORKERS = 16

pinder_index = pd.read_parquet(f"{PATH_TO_PINDER}/index.parquet")
print(f"Loaded {len(pinder_index)} entries from {PATH_TO_PINDER}/index.parquet")

# drop everything with uniprot_L or uniprot_R is UNDEFINED
pinder_index = pinder_index[(pinder_index["uniprot_L"] != "UNDEFINED") & (pinder_index["uniprot_R"] != "UNDEFINED")]
print(f"Remaining entries after dropping UNDEFINED: {len(pinder_index)}")
# retain only the important columns: pdb_id, uniprot_L, uniprot_R, chain_L, chain_R
pinder_index = pinder_index[["pdb_id", "uniprot_L", "uniprot_R", "chain_L", "chain_R"]]
# chain looks like A37, B11, or AN7, we only want to keep the letters
pinder_index["chain_L"] = pinder_index["chain_L"].str.extract(r"([A-Z]+)")
pinder_index["chain_R"] = pinder_index["chain_R"].str.extract(r"([A-Z]+)")
# drop duplicates
pinder_index = pinder_index.drop_duplicates()
print(f"Remaining entries after dropping duplicates: {len(pinder_index)}")

_CHAIN_TOKEN_RE = re.compile(r"^([A-Za-z0-9]+)(?:\[auth [A-Za-z0-9]+])?$")


def fetch_entry_fasta(pdb_id, retries=3):
    """Fetch the per-entry FASTA for a PDB ID, which annotates every chain."""
    url = RCSB_FASTA_URL.format(pdb_id=pdb_id)
    req = urllib.request.Request(url)
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=30, context=_ssl_context()) as resp:
                return resp.read().decode("utf-8")
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return None
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
        except urllib.error.URLError:
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
    return None


def load_rcsb_cache(path):
    """Return {pdb_id: fasta_text} fetched by a previous run, if any.

    Reruns can then skip the network entirely for any pdb_id already here,
    since RCSB's per-entry FASTA is immutable once published.
    """
    cache = {}
    if not os.path.exists(path):
        return cache
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            entry = json.loads(line)
            cache[entry["pdb_id"]] = entry["fasta"]
    return cache


def parse_entry_fasta(fasta_text):
    """Return {chain_id: sequence} for every chain in an RCSB entry FASTA.

    A header's chain field lists one or more chains, each either a bare id
    ("Chain C") or an id followed by its author-assigned alias in brackets
    ("Chains RQ[auth C], SQ[auth D]"). pinder's chain_L/chain_R match the id
    before the bracket, so that's the id used as the dict key.
    """
    chain_to_seq = {}
    header, seq_lines = None, []

    def flush():
        if header is None:
            return
        seq = "".join(seq_lines)
        chains_field = header.split("|")[1].split(maxsplit=1)[1]  # drop "Chain"/"Chains"
        for token in chains_field.split(", "):
            match = _CHAIN_TOKEN_RE.match(token.strip())
            if match:
                chain_to_seq[match.group(1)] = seq

    for line in fasta_text.splitlines():
        if line.startswith(">"):
            flush()
            header, seq_lines = line, []
        elif line.strip():
            seq_lines.append(line.strip())
    flush()
    return chain_to_seq


grouped = pinder_index.groupby("pdb_id")
pdb_ids = list(grouped.groups)
print(f"Fetching sequences for {len(pdb_ids)} unique PDB entries...")

rcsb_cache = load_rcsb_cache(RCSB_CACHE)
print(f"Loaded {len(rcsb_cache)} cached PDB FASTA entries from {RCSB_CACHE}")
to_fetch = [pdb_id for pdb_id in pdb_ids if pdb_id not in rcsb_cache]

if to_fetch:
    print(f"Fetching {len(to_fetch)} PDB entries not in the cache "
          f"(using {RCSB_WORKERS} concurrent workers)...")
    with open(RCSB_CACHE, "a") as cache_fh, ThreadPoolExecutor(max_workers=RCSB_WORKERS) as pool:
        futures = {pool.submit(fetch_entry_fasta, pdb_id): pdb_id for pdb_id in to_fetch}
        for i, future in enumerate(as_completed(futures), start=1):
            pdb_id = futures[future]
            fasta_text = future.result()
            if fasta_text is not None:
                rcsb_cache[pdb_id] = fasta_text
                # flushed immediately so a crash mid-run doesn't lose already-fetched entries
                cache_fh.write(json.dumps({"pdb_id": pdb_id, "fasta": fasta_text}) + "\n")
                cache_fh.flush()
            if i % 100 == 0:
                print(f"  ...{i}/{len(to_fetch)} new PDB entries fetched")

sequences = {}      # id (uniprot:chain) -> sequence
seq_to_id = {}      # (uniprot, sequence) -> id, so identical chains share one id
id_to_uniprot = {}  # id -> the uniprot accession it was resolved from
protein1, protein2 = {}, {}  # row index -> id
missing_pdb, missing_chain = set(), set()


def resolve_chain_id(uniprot, chain, chain_to_seq):
    """Return the id for a uniprot/chain pair, reusing an existing id if this
    uniprot already has an identical sequence under a different chain letter
    (e.g. hemoglobin's identical alpha chains A and C in 1a00), so that one id
    only ever corresponds to one unique sequence.

    Chain letters are only unique within a single PDB entry, not across the
    whole dataset, so the same uniprot can have two genuinely different
    sequences (e.g. different crystallized fragments) that happen to land on
    the same chain letter in two different entries. Detect that and append a
    numeric suffix instead of silently colliding the two into one id.
    """
    seq = chain_to_seq.get(chain)
    if seq is None:
        missing_chain.add(f"{uniprot}:{chain}")
        return f"{uniprot}:{chain}"
    key = (uniprot, seq)
    if key in seq_to_id:
        return seq_to_id[key]
    candidate, suffix = f"{uniprot}:{chain}", 2
    while candidate in sequences and sequences[candidate] != seq:
        candidate = f"{uniprot}:{chain}{suffix}"
        suffix += 1
    seq_to_id[key] = candidate
    sequences[candidate] = seq
    id_to_uniprot[candidate] = uniprot
    return candidate


for i, (pdb_id, rows) in enumerate(grouped, start=1):
    fasta_text = rcsb_cache.get(pdb_id)
    if fasta_text is None:
        missing_pdb.add(pdb_id)
        continue
    chain_to_seq = parse_entry_fasta(fasta_text)
    for row_idx, row in rows.iterrows():
        protein1[row_idx] = resolve_chain_id(row["uniprot_L"], row["chain_L"], chain_to_seq)
        protein2[row_idx] = resolve_chain_id(row["uniprot_R"], row["chain_R"], chain_to_seq)
    if i % 100 == 0:
        print(f"  ...{i}/{len(pdb_ids)} PDB entries processed")

pinder_index["protein1"] = pinder_index.index.map(protein1)
pinder_index["protein2"] = pinder_index.index.map(protein2)

print(f"Collected {len(sequences)} unique sequences")
if missing_pdb:
    print(f"Warning: could not fetch a FASTA for {len(missing_pdb)} PDB entries")
if missing_chain:
    print(f"Warning: {len(missing_chain)} chain ids had no matching sequence in their entry's FASTA")

with open(SEQ_OUT, "w") as fh:
    for uid, seq in sequences.items():
        fh.write(f">{uid}\n{seq}\n")
print(f"Wrote {len(sequences)} sequences to {SEQ_OUT}")

# drop rows where either side's sequence fetch failed, so every id in pinder.csv
# is guaranteed to have a matching entry in the sequence/species/GO files
has_sequence = pinder_index["protein1"].isin(sequences) & pinder_index["protein2"].isin(sequences)
if (~has_sequence).any():
    print(f"Dropping {(~has_sequence).sum()} rows with a missing sequence on protein1 or protein2")
pinder_index = pinder_index[has_sequence]

# redundancy-reduce: protein1-protein2 and protein2-protein1 are the same pair, keep one
pair_key = pinder_index[["protein1", "protein2"]].min(axis=1) + "\t" + pinder_index[["protein1", "protein2"]].max(axis=1)
pinder_index = pinder_index[~pair_key.duplicated(keep="first")]
print(f"Remaining entries after redundancy-reducing protein1/protein2 pairs: {len(pinder_index)}")

pinder_index[["protein1", "protein2", "uniprot_L", "uniprot_R", "pdb_id", "chain_L", "chain_R"]].to_csv(CSV_OUT, index=False)
print(f"Wrote {len(pinder_index)} rows to {CSV_OUT}")

protein_ids = sorted(set(pinder_index["protein1"]) | set(pinder_index["protein2"]))
uniprot_for_id = {pid: id_to_uniprot[pid] for pid in protein_ids}
uniprot_accessions = sorted(set(uniprot_for_id.values()))
print(f"Fetching taxonomy and GO annotations for {len(uniprot_accessions)} unique UniProt accessions "
      f"({len(protein_ids)} ids)...")

canonical_map = build_canonical_map(uniprot_accessions)
canonicals = sorted(canonical_map.keys())

canon_species, canon_go = {}, {}
for i in range(0, len(canonicals), BATCH_SIZE):
    batch = canonicals[i : i + BATCH_SIZE]
    try:
        _, go, sp = parse_batch(fetch_batch(batch, include_sequence=False), include_sequence=False)
    except InvalidAccessionBatch:
        print(f"  Batch rejected (contains an accession UniProt doesn't recognize); "
              f"retrying the {len(batch)} accessions individually...")
        go, sp = {}, {}
        for acc in batch:
            try:
                _, g, s = parse_batch(fetch_batch([acc], include_sequence=False), include_sequence=False)
                go.update(g)
                sp.update(s)
            except InvalidAccessionBatch:
                print(f"Warning: UniProt does not recognize accession {acc}")
            time.sleep(0.1)
    canon_species.update(sp)
    canon_go.update(go)
    print(f"  {min(i + BATCH_SIZE, len(canonicals))}/{len(canonicals)} canonical accessions fetched")
    if i + BATCH_SIZE < len(canonicals):
        time.sleep(0.5)

all_species, all_go = {}, {}
missing_uniprot = set()
for canon, originals in canonical_map.items():
    for acc in originals:
        if canon in canon_species:
            all_species[acc] = canon_species[canon]
            all_go[acc] = canon_go.get(canon, {"BP": set(), "MF": set(), "CC": set()})
        else:
            missing_uniprot.add(acc)
if missing_uniprot:
    print(f"Warning: no taxonomy/GO data found for {len(missing_uniprot)} accessions: {sorted(missing_uniprot)}")

# re-key from uniprot accession to the new ids, so protein_id matches pinder.csv/
# the sequence file: every id sharing a uniprot accession gets that accession's
# taxonomy/GO (GO/taxonomy don't vary by chain or crystallized fragment)
species_by_id = {pid: all_species.get(uniprot_for_id[pid], "") for pid in protein_ids}
go_by_id = {
    pid: all_go.get(uniprot_for_id[pid], {"BP": set(), "MF": set(), "CC": set()})
    for pid in protein_ids
}

write_species_tsv(SPECIES_OUT, protein_ids, species_by_id)
print(f"Wrote species (taxon IDs) for {len(protein_ids)} ids to {SPECIES_OUT}")

write_go_tsv(GO_OUT, protein_ids, go_by_id)
print(f"Wrote GO annotations for {len(protein_ids)} ids to {GO_OUT}")
