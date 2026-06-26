#!/usr/bin/env python3
"""Fetch protein sequences and GO annotations from UniProt for all proteins in a PPI CSV."""

import csv
import ssl
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/stream"
BATCH_SIZE = 100


def _ssl_context():
    try:
        import certifi
        return ssl.create_default_context(cafile=certifi.where())
    except ImportError:
        return ssl.create_default_context()


def get_unique_proteins(ppis_path):
    proteins = set()
    with open(ppis_path) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            proteins.add(row["protein1"].strip())
            proteins.add(row["protein2"].strip())
    return sorted(proteins)


def fetch_batch(accessions, retries=3):
    """Single request returning accession, sequence, and per-category GO IDs as TSV."""
    query = " OR ".join(f"accession:{acc}" for acc in accessions)
    params = urllib.parse.urlencode({
        "query": query,
        "format": "tsv",
        "fields": "accession,sequence,go_p,go_f,go_c",
    })
    req = urllib.request.Request(f"{UNIPROT_URL}?{params}")
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=120, context=_ssl_context()) as resp:
                return resp.read().decode("utf-8")
        except urllib.error.URLError as exc:
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
            else:
                raise RuntimeError(
                    f"Failed to fetch batch after {retries} attempts: {exc}"
                ) from exc


def _parse_go(raw):
    return {t.strip() for t in raw.split(";") if t.strip()}


def parse_batch(text):
    """Return ({accession: sequence}, {accession: {BP/MF/CC: set}}) from UniProt TSV."""
    seqs, go = {}, {}
    lines = text.splitlines()
    for line in lines[1:]:  # skip header row
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        acc = parts[0].strip()
        if acc:
            seqs[acc] = parts[1].strip()
            go[acc] = {
                "BP": _parse_go(parts[2]),
                "MF": _parse_go(parts[3]),
                "CC": _parse_go(parts[4]),
            }
    return seqs, go


def main():
    if len(sys.argv) != 4:
        sys.exit(f"Usage: {sys.argv[0]} ppis.csv sequences.fasta go_annotations.tsv")

    ppis_path, fasta_out, go_out = sys.argv[1], sys.argv[2], sys.argv[3]
    proteins = get_unique_proteins(ppis_path)
    print(
        f"Fetching sequences and GO annotations for {len(proteins)} proteins from UniProt...",
        file=sys.stderr,
    )

    all_seqs, all_go = {}, {}
    for i in range(0, len(proteins), BATCH_SIZE):
        batch = proteins[i : i + BATCH_SIZE]
        seqs, go = parse_batch(fetch_batch(batch))
        all_seqs.update(seqs)
        all_go.update(go)
        print(
            f"  {min(i + BATCH_SIZE, len(proteins))}/{len(proteins)} fetched",
            file=sys.stderr,
        )
        if i + BATCH_SIZE < len(proteins):
            time.sleep(0.5)

    missing = set(proteins) - set(all_seqs)
    if missing:
        print(
            f"Warning: no data found for {len(missing)} proteins: {sorted(missing)}",
            file=sys.stderr,
        )

    with open(fasta_out, "w") as fh:
        for acc in proteins:
            if acc in all_seqs:
                fh.write(f">{acc}\n{all_seqs[acc]}\n")

    with open(go_out, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["protein_id", "go_bp", "go_mf", "go_cc"])
        for acc in proteins:
            cats = all_go.get(acc, {"BP": set(), "MF": set(), "CC": set()})
            writer.writerow([
                acc,
                ";".join(sorted(cats.get("BP", set()))),
                ";".join(sorted(cats.get("MF", set()))),
                ";".join(sorted(cats.get("CC", set()))),
            ])

    print(f"Written {len(all_seqs)} sequences to {fasta_out}", file=sys.stderr)
    print(f"Written GO annotations for {len(proteins)} proteins to {go_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
