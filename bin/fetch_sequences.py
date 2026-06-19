#!/usr/bin/env python3
"""Fetch protein sequences from UniProt for all proteins listed in a PPI CSV."""

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


def fetch_fasta_batch(accessions, retries=3):
    query = " OR ".join(f"accession:{acc}" for acc in accessions)
    params = urllib.parse.urlencode({"query": query, "format": "fasta"})
    req = urllib.request.Request(f"{UNIPROT_URL}?{params}")
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=120, context=_ssl_context()) as resp:
                return resp.read().decode("utf-8")
        except urllib.error.URLError as exc:
            if attempt < retries - 1:
                time.sleep(2**attempt)
            else:
                raise RuntimeError(
                    f"Failed to fetch batch after {retries} attempts: {exc}"
                ) from exc


def parse_uniprot_fasta(text):
    """Return {accession: sequence} parsed from UniProt FASTA text."""
    result = {}
    acc = None
    parts = []
    for line in text.splitlines():
        if line.startswith(">"):
            if acc:
                result[acc] = "".join(parts)
            # Header forms: >sp|ACC|NAME ... or >tr|ACC|NAME ... or >ACC ...
            header = line[1:]
            fields = header.split("|")
            if len(fields) >= 2 and fields[0] in ("sp", "tr"):
                acc = fields[1]
            else:
                acc = header.split()[0]
            parts = []
        elif line:
            parts.append(line)
    if acc and parts:
        result[acc] = "".join(parts)
    return result


def main():
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} ppis.csv sequences.fasta")

    ppis_path, out_path = sys.argv[1], sys.argv[2]
    proteins = get_unique_proteins(ppis_path)
    print(f"Fetching {len(proteins)} unique proteins from UniProt...", file=sys.stderr)

    all_seqs = {}
    for i in range(0, len(proteins), BATCH_SIZE):
        batch = proteins[i : i + BATCH_SIZE]
        text = fetch_fasta_batch(batch)
        all_seqs.update(parse_uniprot_fasta(text))
        print(
            f"  {min(i + BATCH_SIZE, len(proteins))}/{len(proteins)} fetched",
            file=sys.stderr,
        )
        if i + BATCH_SIZE < len(proteins):
            time.sleep(0.5)

    missing = set(proteins) - set(all_seqs)
    if missing:
        print(
            f"Warning: no sequence found for {len(missing)} proteins: {sorted(missing)}",
            file=sys.stderr,
        )

    with open(out_path, "w") as fh:
        for acc in proteins:
            if acc in all_seqs:
                fh.write(f">{acc}\n{all_seqs[acc]}\n")

    print(f"Written {len(all_seqs)} sequences to {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
