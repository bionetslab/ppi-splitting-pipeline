#!/usr/bin/env python3
"""Fetch protein sequences and GO annotations from UniProt for all proteins in a PPI CSV."""

import csv
import ssl
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

UNIPROT_URL      = "https://rest.uniprot.org/uniprotkb/stream"
UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"
BATCH_SIZE = 100


def _ssl_context():
    try:
        import certifi
        return ssl.create_default_context(cafile=certifi.where())
    except ImportError:
        return ssl.create_default_context()


def _canonical(acc):
    """Strip isoform suffix: 'O14836-2' → 'O14836'."""
    return acc.split("-")[0]


def get_unique_proteins(ppis_path):
    """Return (sorted original IDs, {canonical_id: [original_ids]})."""
    proteins = set()
    with open(ppis_path) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            proteins.add(row["protein1"].strip())
            proteins.add(row["protein2"].strip())
    canonical_map = {}
    for acc in proteins:
        canonical_map.setdefault(_canonical(acc), []).append(acc)
    return sorted(proteins), canonical_map


def fetch_isoform_sequence(acc, retries=3):
    """Fetch sequence for one isoform accession via the per-entry FASTA endpoint.

    Returns the sequence string, or None if not found.
    """
    url = UNIPROT_FASTA_URL.format(acc=acc)
    req = urllib.request.Request(url)
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=30, context=_ssl_context()) as resp:
                fasta = resp.read().decode("utf-8").strip()
                lines = fasta.splitlines()
                seq = "".join(l for l in lines if not l.startswith(">"))
                return seq if seq else None
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return None
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
        except urllib.error.URLError:
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
    return None


def fetch_batch(accessions, retries=3):
    """Batch request returning accession, sequence, GO IDs, and taxon ID as TSV.

    accessions should be canonical (no isoform suffix) to ensure UniProt matches them.
    """
    query = " OR ".join(f"accession:{acc}" for acc in accessions)
    params = urllib.parse.urlencode({
        "query": query,
        "format": "tsv",
        "fields": "accession,sequence,go_p,go_f,go_c,organism_id",
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
    """Return ({accession: sequence}, {accession: {BP/MF/CC: set}}, {accession: taxon_id})."""
    seqs, go, species = {}, {}, {}
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
            species[acc] = parts[5].strip() if len(parts) > 5 else ""
    return seqs, go, species


def main():
    if len(sys.argv) != 5:
        sys.exit(f"Usage: {sys.argv[0]} ppis.csv sequences.fasta go_annotations.tsv species.tsv")

    ppis_path, fasta_out, go_out, species_out = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
    proteins, canonical_map = get_unique_proteins(ppis_path)
    canonicals = sorted(canonical_map.keys())
    print(
        f"Fetching data for {len(proteins)} proteins "
        f"({len(canonicals)} canonical accessions) from UniProt...",
        file=sys.stderr,
    )

    # Fetch by canonical accession; UniProt ignores isoform suffixes in queries.
    canon_seqs, canon_go, canon_species = {}, {}, {}
    for i in range(0, len(canonicals), BATCH_SIZE):
        batch = canonicals[i : i + BATCH_SIZE]
        seqs, go, species = parse_batch(fetch_batch(batch))
        canon_seqs.update(seqs)
        canon_go.update(go)
        canon_species.update(species)
        print(
            f"  {min(i + BATCH_SIZE, len(canonicals))}/{len(canonicals)} canonical accessions fetched",
            file=sys.stderr,
        )
        if i + BATCH_SIZE < len(canonicals):
            time.sleep(0.5)

    # Expand results back to all original IDs; GO and species come from canonical entry.
    all_seqs, all_go, all_species = {}, {}, {}
    for canon, originals in canonical_map.items():
        for acc in originals:
            if canon in canon_seqs:
                all_seqs[acc]    = canon_seqs[canon]
                all_go[acc]      = canon_go[canon]
                all_species[acc] = canon_species.get(canon, "")

    # Isoforms may have distinct sequences — fetch each individually.
    isoforms = [acc for acc in proteins if "-" in acc]
    if isoforms:
        print(f"Fetching isoform-specific sequences for {len(isoforms)} isoforms...",
              file=sys.stderr)
        for acc in isoforms:
            seq = fetch_isoform_sequence(acc)
            if seq:
                all_seqs[acc] = seq
            time.sleep(0.1)

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

    with open(species_out, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["protein_id", "taxon_id"])
        for acc in proteins:
            writer.writerow([acc, all_species.get(acc, "")])

    print(f"Written {len(all_seqs)} sequences to {fasta_out}", file=sys.stderr)
    print(f"Written GO annotations for {len(proteins)} proteins to {go_out}", file=sys.stderr)
    print(f"Written species (taxon IDs) for {len(proteins)} proteins to {species_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
