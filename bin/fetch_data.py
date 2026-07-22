#!/usr/bin/env python3
"""Fetch protein sequences and GO annotations from UniProt for all proteins in a PPI CSV."""

import csv
import ssl
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

UNIPROT_URL      = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"
UNIPARC_URL      = "https://rest.uniprot.org/uniparc/search"
BATCH_SIZE = 100
SEARCH_PAGE_SIZE = 500  # UniProt's documented max page size for /uniprotkb/search


class InvalidAccessionBatch(Exception):
    """Raised when UniProt rejects a whole batch (HTTP 400) because it contains
    an accession that isn't valid UniProtKB format, e.g. a UniParc-only ID such
    as a raw EMBL/GenBank protein_id. Not transient, so callers should not retry
    the same batch and should instead fall back to per-accession requests.
    """


def _ssl_context():
    try:
        import certifi
        return ssl.create_default_context(cafile=certifi.where())
    except ImportError:
        return ssl.create_default_context()


def _canonical(acc):
    """Strip isoform suffix: 'O14836-2' → 'O14836'."""
    return acc.split("-")[0]


def read_protein_list(path):
    """Return sorted unique protein IDs from a plain text file, one ID per line."""
    with open(path) as fh:
        proteins = {line.strip() for line in fh if line.strip()}
    return sorted(proteins)


def build_canonical_map(proteins):
    """Return {canonical_id: [original_ids]} for a collection of protein IDs."""
    canonical_map = {}
    for acc in proteins:
        canonical_map.setdefault(_canonical(acc), []).append(acc)
    return canonical_map


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


_STREAM_ERROR_MARKER = "error encountered when streaming data"


def _is_stream_error(text: str) -> bool:
    """True if `text` is UniProt's own error page for a request that failed
    server-side after the response had already started (HTTP 200), so it
    never raises HTTPError/URLError and would otherwise look like a
    successful-but-unparseable TSV response."""
    return _STREAM_ERROR_MARKER in text.lower()


def _next_page_url(link_header):
    """Parse the RFC 5988 `Link` header UniProt returns on /search responses
    and return the rel="next" URL, or None if this was the last page."""
    if not link_header:
        return None
    for part in link_header.split(","):
        segments = [s.strip() for s in part.split(";")]
        if len(segments) >= 2 and segments[1] == 'rel="next"' and segments[0].startswith("<"):
            return segments[0].strip("<>")
    return None


def _fetch_page(url, retries=3):
    """Fetch one page from /uniprotkb/search, returning (body_text, next_page_url_or_None)."""
    req = urllib.request.Request(url)
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=120, context=_ssl_context()) as resp:
                text = resp.read().decode("utf-8")
                next_url = _next_page_url(resp.headers.get("Link"))
        except urllib.error.HTTPError as exc:
            if exc.code == 400:
                # UniProt rejects the whole OR-joined query if one accession is
                # non-UniProtKB format (e.g. a UniParc/EMBL id). Retrying won't help.
                raise InvalidAccessionBatch(str(exc)) from exc
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
                continue
            raise RuntimeError(
                f"Failed to fetch page after {retries} attempts: {exc}"
            ) from exc
        except urllib.error.URLError as exc:
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
                continue
            raise RuntimeError(
                f"Failed to fetch page after {retries} attempts: {exc}"
            ) from exc
        else:
            if _is_stream_error(text):
                if attempt < retries - 1:
                    time.sleep(2 ** attempt)
                    continue
                raise RuntimeError(
                    f"UniProt's search endpoint failed after {retries} attempts: "
                    f"{text.strip()[:200]!r}"
                )
            return text, next_url
    raise RuntimeError(f"fetch_page called with retries={retries} <= 0")


def fetch_batch(accessions, retries=3):
    """Batch request returning accession, sequence, GO IDs, and taxon ID as TSV.

    accessions should be canonical (no isoform suffix) to ensure UniProt matches them.
    Uses UniProt's paginated /uniprotkb/search endpoint (per their guidance for
    fetching large numbers of results, https://www.uniprot.org/help/pagination)
    rather than /stream, which holds one connection open for the whole result
    set and is more prone to being cut off by transient server-side errors.
    Follows the `Link: rel="next"` cursor until UniProt reports no more pages
    (in practice always one page here, since BATCH_SIZE <= SEARCH_PAGE_SIZE).
    """
    query = " OR ".join(f"accession:{acc}" for acc in accessions)
    params = urllib.parse.urlencode({
        "query": query,
        "format": "tsv",
        "fields": "accession,sequence,go_p,go_f,go_c,organism_id",
        "size": SEARCH_PAGE_SIZE,
    })
    url = f"{UNIPROT_URL}?{params}"

    header = None
    rows = []
    while url:
        text, url = _fetch_page(url, retries)
        lines = text.splitlines()
        if not lines:
            break
        if header is None:
            header = lines[0]
        rows.extend(lines[1:])

    return "\n".join([header or "", *rows])


def fetch_uniparc_entry(acc, retries=3):
    """Look up an accession UniProtKB doesn't recognize (e.g. a raw EMBL/GenBank
    protein_id) via UniParc's cross-reference index instead.

    Returns (sequence, taxon_id), or None if UniParc has no record either.
    UniParc is a sequence archive, not a curated database, so it has no GO
    annotations - callers only get sequence + taxonomy back for these IDs.
    """
    params = urllib.parse.urlencode({
        "query": f"dbid:{acc}",
        "format": "tsv",
        "fields": "organism_id,sequence",
    })
    req = urllib.request.Request(f"{UNIPARC_URL}?{params}")
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=30, context=_ssl_context()) as resp:
                lines = resp.read().decode("utf-8").strip().splitlines()
                if len(lines) < 2:
                    return None
                taxon_id, seq = lines[1].split("\t")
                return seq.strip(), taxon_id.strip()
        except urllib.error.HTTPError as exc:
            if exc.code == 400:
                return None
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
        except urllib.error.URLError:
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
    return None


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
        sys.exit(f"Usage: {sys.argv[0]} proteins.txt sequences.fasta go_annotations.tsv species.tsv")

    proteins_path, fasta_out, go_out, species_out = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
    proteins = read_protein_list(proteins_path)
    canonical_map = build_canonical_map(proteins)
    canonicals = sorted(canonical_map.keys())
    print(
        f"Fetching data for {len(proteins)} proteins "
        f"({len(canonicals)} canonical accessions) from UniProt...",
        file=sys.stderr,
    )

    # Fetch by canonical accession; UniProt ignores isoform suffixes in queries.
    canon_seqs, canon_go, canon_species = {}, {}, {}
    uniparc_only = []  # accessions UniProtKB rejects outright; resolved via UniParc below
    for i in range(0, len(canonicals), BATCH_SIZE):
        batch = canonicals[i : i + BATCH_SIZE]
        try:
            seqs, go, species = parse_batch(fetch_batch(batch))
        except InvalidAccessionBatch:
            print(
                f"  Batch rejected (contains an accession UniProt doesn't recognize); "
                f"retrying the {len(batch)} accessions individually to isolate it...",
                file=sys.stderr,
            )
            seqs, go, species = {}, {}, {}
            for acc in batch:
                try:
                    s, g, sp = parse_batch(fetch_batch([acc]))
                    seqs.update(s)
                    go.update(g)
                    species.update(sp)
                except InvalidAccessionBatch:
                    print(f"Failed to fetch {acc}; falling back to UniParc...", file=sys.stderr)
                    uniparc_only.append(acc)
                time.sleep(0.1)
        canon_seqs.update(seqs)
        canon_go.update(go)
        canon_species.update(species)
        print(
            f"  {min(i + BATCH_SIZE, len(canonicals))}/{len(canonicals)} canonical accessions fetched",
            file=sys.stderr,
        )
        if i + BATCH_SIZE < len(canonicals):
            time.sleep(0.5)

    if uniparc_only:
        print(
            f"Falling back to UniParc for {len(uniparc_only)} accession(s) not valid in "
            f"UniProtKB (sequence + taxonomy only, no GO annotations available there): "
            f"{sorted(uniparc_only)}",
            file=sys.stderr,
        )
        for acc in uniparc_only:
            result = fetch_uniparc_entry(acc)
            if result:
                canon_seqs[acc], canon_species[acc] = result
            time.sleep(0.2)

    # Expand results back to all original IDs; GO and species come from canonical entry.
    all_seqs, all_go, all_species = {}, {}, {}
    for canon, originals in canonical_map.items():
        for acc in originals:
            if canon in canon_seqs:
                all_seqs[acc]    = canon_seqs[canon]
                all_go[acc]      = canon_go.get(canon, {"BP": set(), "MF": set(), "CC": set()})
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
