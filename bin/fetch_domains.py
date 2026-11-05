#!/usr/bin/env python3
"""Fetch Pfam domain instances for a list of Pfam families (DDI mode's FETCH_DATA).

Up to --pool-size instances are sampled per family in --families, stratified by --tiers, and written as
sequences.fasta, species.tsv, an empty go_annotations.tsv, instances.tsv (instance ->
family, clan, parent protein) and dropped_families.tsv (families that kept no instance).

Inputs, downloaded unless given locally:
    Pfam-A.regions.tsv.gz  family, parent accession and alignment coordinates per region
    uniprot_*.dat.gz       parent sequence, taxon and review flag; repeatable, Swiss-Prot first
    Pfam-A.clans.tsv.gz    family -> clan; a blank clan means a singleton clan

The domain is cut as sequence[ali_start-1:ali_end]. Every written accession is a UniProt
primary accession; regions on a demoted (secondary) accession are dropped. Instance ids
are {family}_{protein}_{start}_{end}, never '|'-separated; parse back with rsplit("_", 3).
"""

import argparse
import csv
import gzip
import hashlib
import io
import os
import random
import re
import shutil
import sys
import time
import urllib.request
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fetch_data import _ssl_context, write_go_tsv, write_species_tsv
from utils import INSTANCE_COLUMNS, write_fasta

PFAM_BASE = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release"
PFAM_REGIONS_URL = f"{PFAM_BASE}/Pfam-A.regions.tsv.gz"
PFAM_CLANS_URL = f"{PFAM_BASE}/Pfam-A.clans.tsv.gz"
PFAM_DEAD_URL = f"{PFAM_BASE}/Pfam-A.dead.gz"
PFAM_VERSION_URL = f"{PFAM_BASE}/Pfam.version.gz"
SPROT_URL = (
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/" "complete/uniprot_sprot.dat.gz"
)

# Bump whenever the sampling or output format changes, so stale caches are ignored.
FORMAT_VERSION = 4

# Pfam-A.regions.tsv.gz columns, read by name: the file has gained columns before.
REGIONS_FAMILY_COL = "pfamA_acc"
REGIONS_PROTEIN_COL = "pfamseq_acc"
REGIONS_START_COL = "ali_start"
REGIONS_END_COL = "ali_end"
REGIONS_REQUIRED = (REGIONS_FAMILY_COL, REGIONS_PROTEIN_COL, REGIONS_START_COL, REGIONS_END_COL)
# Optional checksum of the sequence Pfam's coordinates refer to.
REGIONS_CRC64_COL = "crc64"
CRC64_RE = re.compile(r"^[0-9A-Fa-f]{16}$")
# Minimum checked regions before an all-mismatch result is treated as a format change.
CRC64_SANITY_MIN = 100

HUMAN_TAXON = "9606"

# At most pool_size // PARENT_CAP_DIVISOR regions of one parent are offered per (family, tier),
# so a repeat-domain family is not filled from a single protein.
PARENT_CAP_DIVISOR = 5

# OX   NCBI_TaxID=9606;
DAT_TAXON_RE = re.compile(r"NCBI_TaxID=(\d+)")

# ID   TP53_HUMAN              Reviewed;         393 AA.
# ID   A0A804HJH4_HUMAN        Unreviewed;       557 AA.
DAT_REVIEW_TOKENS = {"Reviewed;": True, "Unreviewed;": False}

DEAD_AC_RE = re.compile(r"^#=GF\s+AC\s+(PF\d+)")

# Disjoint strata, in fill order. Reviewed outranks human.
TIERS = ("human_reviewed", "other_reviewed", "human_unreviewed", "other_unreviewed")
N_TIERS = len(TIERS)

# Must match TIERS: a wrong set silently makes a human-only run parse all of Swiss-Prot.
HUMAN_TIERS = frozenset({0, 2})
UNREVIEWED_TIERS = frozenset({2, 3})
NONHUMAN_TIERS = frozenset({1, 3})

DROPPED_FAMILIES_COLUMNS = ("family", "reason")
DROP_REASONS = {
    "no_eligible_instances": "have Pfam records but none in a stratum --tiers names",
    "demoted_accession": "have Pfam records, but every one names an accession UniProt has demoted to secondary",
    "sequence_mismatch": "have eligible Pfam records, but every one disagrees with the UniProt sequence "
    "(CRC64 differs or coordinates fall outside it)",
    "dead": "are dead in Pfam",
    "not_in_pfam": "have no regions in Pfam-A.regions and are not listed as dead",
}


# ---------------------------------------------------------------------------
# Downloads
# ---------------------------------------------------------------------------


def _open_url(url, retries=3, timeout=120):
    """Open a URL for streaming, retrying transient failures (OSError) with backoff."""
    req = urllib.request.Request(url, headers={"Accept-Encoding": "identity"})
    for attempt in range(retries):
        try:
            return urllib.request.urlopen(req, timeout=timeout, context=_ssl_context())
        except OSError as exc:
            if attempt < retries - 1:
                time.sleep(2**attempt)
                continue
            raise RuntimeError(f"Failed to fetch {url} after {retries} attempts: {exc}") from exc


def download_text(url, gzipped=False):
    """Return a small download as text, decompressing gzip in memory."""
    with _open_url(url) as resp:
        raw = resp.read()
    if gzipped:
        raw = gzip.decompress(raw)
    return raw.decode("utf-8", "replace")


def download_to_file(url, path):
    """Stream a download to disk."""
    tmp = f"{path}.part"
    with _open_url(url) as resp, open(tmp, "wb") as fh:
        shutil.copyfileobj(resp, fh, length=1 << 20)
    os.replace(tmp, path)
    return path


def iter_gzip_lines(fileobj):
    """Yield raw bytes lines from a gzip stream (HTTP response or local file)."""
    gz = gzip.GzipFile(fileobj=fileobj, mode="rb")
    return io.BufferedReader(gz, buffer_size=1 << 20)


# ---------------------------------------------------------------------------
# Reference tables
# ---------------------------------------------------------------------------


def parse_release(text):
    """Return the Pfam release string, e.g. '38.2', from Pfam.version's contents."""
    for line in text.splitlines():
        if line.lower().startswith("pfam release"):
            return line.split(":", 1)[1].strip()
    raise RuntimeError(f"Could not find a release line in Pfam.version: {text.strip()[:200]!r}")


def parse_clans(text):
    """Return {family: clan}. A blank clan column means the family is its own singleton clan."""
    clans = {}
    for row in csv.reader(io.StringIO(text), delimiter="\t"):
        if not row or not row[0].strip():
            continue
        fam = row[0].strip()
        clan = row[1].strip() if len(row) > 1 else ""
        clans[fam] = clan if clan else f"CL_{fam}"
    return clans


def parse_uniprot_dat(stream, taxa=None, into=None, secondary=None):
    """Fill {primary accession: (taxon_id, sequence, reviewed)} from a UniProt flat-file stream.

    Only the first accession of the first AC line is indexed; the rest go into the
    optional `secondary` set. `taxa` restricts to those taxon ids, and with `into`
    the first writer wins on a duplicate accession.
    """
    universe = {} if into is None else into
    accessions, taxon, seq_parts, in_seq, reviewed = [], "", [], False, None

    def flush():
        if not accessions:
            return
        if reviewed is None:
            raise RuntimeError(
                f"UniProt entry {accessions[0]} has no ID line, so its review status is "
                "unknown -- refusing to guess it"
            )
        if taxon and (taxa is None or taxon in taxa):
            sequence = "".join(seq_parts)
            universe.setdefault(accessions[0], (taxon, sequence, reviewed))
            if secondary is not None:
                secondary.update(accessions[1:])

    for raw in stream:
        line = raw.decode("ascii", "replace") if isinstance(raw, bytes) else raw
        if line.startswith("//"):
            flush()
            accessions, taxon, seq_parts, in_seq, reviewed = [], "", [], False, None
        elif in_seq:
            seq_parts.append(line.strip().replace(" ", ""))
        elif line.startswith("ID   "):
            fields = line.split()
            token = fields[2] if len(fields) > 2 else ""
            if token not in DAT_REVIEW_TOKENS:
                raise RuntimeError(f"UniProt ID line carries neither Reviewed; nor Unreviewed;: {line.rstrip()!r}")
            reviewed = DAT_REVIEW_TOKENS[token]
        elif line.startswith("AC   "):
            accessions.extend(a.strip() for a in line[5:].split(";") if a.strip())
        elif line.startswith("OX   "):
            match = DAT_TAXON_RE.search(line)
            taxon = match.group(1) if match else ""
        elif line.startswith("SQ   "):
            in_seq = True
    flush()
    return universe


def load_universe(paths, taxa=None):
    """Return (universe, secondary) from every --uniprot-dat; the first file to carry an accession wins.

    A file that contributes no entry after the taxon filter is fatal.
    """
    universe, secondary = {}, set()
    for path in paths:
        before, sec_before = len(universe), len(secondary)
        with open(path, "rb") as fh:
            parse_uniprot_dat(iter_gzip_lines(fh), taxa, universe, secondary)
        added = len(universe) - before
        print(
            f"  {path}: {added:,} entries (one primary accession each), "
            f"{len(secondary) - sec_before:,} secondary accessions noted",
            file=sys.stderr,
        )
        if not added:
            raise RuntimeError(
                f"{path} contributed no proteins to the universe"
                + (f" for taxon(s) {', '.join(sorted(taxa))}" if taxa else "")
                + (
                    " -- already covered by an earlier --uniprot-dat?"
                    if before
                    else " -- every region would be dropped as not_in_universe"
                )
            )
    return universe, secondary


def eligible_taxa(eligible):
    """Return {'9606'} when no eligible stratum can hold a non-human protein, else None."""
    return {HUMAN_TAXON} if set(eligible) <= HUMAN_TIERS else None


def check_universe_covers(eligible, universe, sources):
    """Fail when an eligible stratum can never be filled from the parsed universe."""
    where = ", ".join(str(p) for p in sources)
    for tiers, holds, hint in (
        (
            UNREVIEWED_TIERS,
            any(not rec[2] for rec in universe.values()),
            "Pass the matching TrEMBL flat file (--uniprot-dat uniprot_trembl_*.dat.gz).",
        ),
        (
            NONHUMAN_TIERS,
            any(rec[0] != HUMAN_TAXON for rec in universe.values()),
            "Pass a flat file that is not restricted to Homo sapiens.",
        ),
    ):
        wanted = sorted(set(eligible) & tiers)
        if not wanted or holds:
            continue
        adjective = "unreviewed" if tiers is UNREVIEWED_TIERS else "non-human"
        raise RuntimeError(
            f"--tiers includes {', '.join(TIERS[t] for t in wanted)}, but none of the "
            f"{len(universe):,} proteins parsed from {where} are {adjective}. {hint}"
        )


def parse_tiers(value):
    """Parse --tiers (stratum names) into a frozenset of indices into TIERS."""
    names = [name.strip() for name in str(value).split(",") if name.strip()]
    if not names:
        raise argparse.ArgumentTypeError(f"--tiers is empty; name at least one of {', '.join(TIERS)}")
    unknown = [name for name in names if name not in TIERS]
    if unknown:
        raise argparse.ArgumentTypeError(
            f"--tiers has unknown stratum name(s) {', '.join(unknown)}; valid names are {', '.join(TIERS)}"
        )
    return frozenset(TIERS.index(name) for name in names)


def parse_dead(text):
    """Return the set of killed Pfam families."""
    return {match.group(1) for line in text.splitlines() if (match := DEAD_AC_RE.match(line))}


# ---------------------------------------------------------------------------
# Reservoir sampling
# ---------------------------------------------------------------------------


class Reservoir:
    """Algorithm-R reservoir over one (family, tier) stratum.

    Seeded per (seed, family, tier), so a family's sample does not depend on its
    position in the stream.
    """

    __slots__ = ("k", "rng", "buf", "seen")

    def __init__(self, k, seed, family, tier):
        self.k = k
        self.rng = random.Random(f"{seed}:{family}:{tier}")
        self.buf = []
        self.seen = 0

    def offer(self):
        """Reserve a slot for the next record, or return None if it is not sampled."""
        self.seen += 1
        if len(self.buf) < self.k:
            self.buf.append(None)
            return len(self.buf) - 1
        j = self.rng.randrange(self.seen)
        return j if j < self.k else None


_CRC64_POLY = 0xD800000000000000  # ISO 3309, reflected


def _crc64_table():
    table = []
    for i in range(256):
        crc = i
        for _ in range(8):
            crc = (crc >> 1) ^ _CRC64_POLY if crc & 1 else crc >> 1
        table.append(crc)
    return table


_CRC64_TABLE = _crc64_table()


def crc64(sequence):
    """UniProt's CRC64 of a sequence, as the 16-digit uppercase hex on its SQ line."""
    crc = 0
    for byte in sequence.encode("ascii", "replace"):
        crc = _CRC64_TABLE[(crc ^ byte) & 0xFF] ^ (crc >> 8)
    return f"{crc:016X}"


def tier_of(taxon_id, reviewed):
    """Classify one instance into a disjoint tier index (0 = most preferred)."""
    is_human = taxon_id == HUMAN_TAXON
    if reviewed:
        return 0 if is_human else 1
    return 2 if is_human else 3


# ---------------------------------------------------------------------------
# The streaming pass
# ---------------------------------------------------------------------------


def sample_instances(stream, wanted, pool_size, universe, seed, eligible=None, secondary=None):
    """One pass over Pfam-A.regions.tsv, returning ({family: [record]}, stats, seen, family_stats).

    Regions whose parent is outside `universe`, whose tier is not in `eligible`, or whose
    coordinates/crc64 disagree with the UniProt sequence are counted and skipped; the rest
    go through a per-parent cap into one Reservoir per (family, tier). Each family is then
    filled from tier 0 upward up to `pool_size`, preferring distinct parent proteins.
    """
    if eligible is None:
        eligible = frozenset(range(N_TIERS))
    reservoirs = {}
    parent_cap = max(1, pool_size // PARENT_CAP_DIVISOR)
    parent_offers = defaultdict(Counter)
    stats = Counter()
    family_stats = {"demoted": Counter(), "mismatch": Counter(), "eligible_tier": Counter()}
    parent_crc = {}
    seen_families = set()
    last_family = None
    columns = None
    crc_col = None

    for raw in stream:
        line = raw.decode("utf-8", "replace") if isinstance(raw, bytes) else raw
        fields = line.rstrip("\r\n").split("\t")

        if columns is None:
            # Header required: a positional read would silently swap envelope and alignment coordinates.
            columns = {name: i for i, name in enumerate(f.strip() for f in fields)}
            missing = [c for c in REGIONS_REQUIRED if c not in columns]
            if missing:
                raise RuntimeError(f"Pfam-A.regions is missing column(s) {', '.join(missing)}; " f"header was {fields}")
            crc_col = columns.get(REGIONS_CRC64_COL)
            continue

        stats["records"] += 1
        if len(fields) <= max(columns[c] for c in REGIONS_REQUIRED):
            stats["malformed_rows"] += 1
            continue

        family = fields[columns[REGIONS_FAMILY_COL]].strip()
        if family not in wanted:
            continue

        if family != last_family:
            if family in seen_families:
                stats["noncontiguous_blocks"] += 1
            seen_families.add(family)
            last_family = family

        accession = fields[columns[REGIONS_PROTEIN_COL]].strip()
        entry = universe.get(accession)
        if entry is None:
            if secondary and accession in secondary:
                stats["demoted_accession"] += 1
                family_stats["demoted"][family] += 1
            else:
                stats["not_in_universe"] += 1
            continue
        taxon_id, parent, reviewed = entry

        tier = tier_of(taxon_id, reviewed)
        stats[f"tier_{TIERS[tier]}_offered"] += 1
        if tier not in eligible:
            stats["ineligible_tier_records"] += 1
            continue

        start_s, end_s = fields[columns[REGIONS_START_COL]].strip(), fields[columns[REGIONS_END_COL]].strip()
        try:
            start_i, end_i = int(start_s), int(end_s)
        except ValueError:
            stats["malformed_rows"] += 1
            continue
        # 1-based inclusive.
        if start_i < 1 or end_i > len(parent) or end_i < start_i:
            stats["coords_out_of_range"] += 1
            family_stats["mismatch"][family] += 1
            continue
        if crc_col is not None and crc_col < len(fields) and CRC64_RE.match(want := fields[crc_col].strip()):
            stats["crc64_checked"] += 1
            have = parent_crc.get(accession)
            if have is None:
                have = parent_crc[accession] = crc64(parent)
            if have != want.upper():
                stats["crc64_mismatch"] += 1
                family_stats["mismatch"][family] += 1
                continue
        sequence = parent[start_i - 1 : end_i]
        family_stats["eligible_tier"][family] += 1

        key = (family, tier)
        # Cap rows offered per parent (see PARENT_CAP_DIVISOR).
        offered = parent_offers[key]
        if offered[accession] >= parent_cap:
            stats["parent_cap_skipped"] += 1
            continue
        offered[accession] += 1

        res = reservoirs.get(key)
        if res is None:
            res = reservoirs[key] = Reservoir(pool_size, seed, family, TIERS[tier])
        slot = res.offer()
        if slot is None:
            continue

        res.buf[slot] = {
            "instance_id": f"{family}_{accession}_{start_s}_{end_s}",
            "family": family,
            "protein_id": accession,
            "taxon_id": taxon_id,
            "start": start_s,
            "end": end_s,
            "tier": tier,
            "reviewed": reviewed,
            "sequence": sequence,
        }

    if columns is None:
        raise RuntimeError("Pfam-A.regions stream was empty -- not even a header row")
    if stats["crc64_checked"] >= CRC64_SANITY_MIN and stats["crc64_mismatch"] == stats["crc64_checked"]:
        # Every checksum disagreeing means a format change, not a release skew.
        raise RuntimeError(
            f"all {stats['crc64_checked']:,} Pfam crc64 values disagree with the UniProt sequences' CRC64 "
            "-- has the regions file's crc64 column changed format?"
        )

    # Fill each family from tier 0 upward until pool_size; a non-empty tier does not stop
    # lower tiers from filling the remaining room. Oversized tiers are sampled, not truncated.
    by_family = {}
    for family in sorted(seen_families):
        rng = random.Random(f"{seed}:{family}:pick")
        chosen = []
        for tier in range(N_TIERS):
            room = pool_size - len(chosen)
            if room <= 0:
                break
            res = reservoirs.get((family, tier))
            if res is None:
                continue
            buf = res.buf
            picked = buf if len(buf) <= room else _spread_over_parents(rng, buf, room)
            chosen.extend(picked)
        if chosen:
            by_family[family] = sorted(chosen, key=lambda r: r["instance_id"])
            for rec in by_family[family]:
                stats[f"tier_{TIERS[rec['tier']]}_kept"] += 1

    return by_family, stats, seen_families, family_stats


def _spread_over_parents(rng, records, room):
    """Pick `room` records, round-robin over shuffled per-parent queues to maximise distinct parents."""
    by_parent = defaultdict(list)
    for rec in records:
        by_parent[rec["protein_id"]].append(rec)
    queues = [by_parent[p] for p in sorted(by_parent)]
    for queue in queues:
        rng.shuffle(queue)
    rng.shuffle(queues)

    picked = []
    while len(picked) < room:
        for queue in queues:
            if not queue:
                continue
            picked.append(queue.pop())
            if len(picked) == room:
                break
    return picked


# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------


def write_instances_tsv(path, by_family, clans):
    """Write instances.tsv (utils.INSTANCE_COLUMNS)."""
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(INSTANCE_COLUMNS)
        for family in sorted(by_family):
            for rec in by_family[family]:
                writer.writerow(
                    [
                        rec["instance_id"],
                        family,
                        clans.get(family, f"CL_{family}"),
                        rec["protein_id"],
                        rec["start"],
                        rec["end"],
                        rec["taxon_id"],
                        "reviewed" if rec["reviewed"] else "unreviewed",
                    ]
                )


def classify_dropped(families, by_family, seen_families, dead, family_stats):
    """Return [(family, reason)] for every requested family that kept no instance (see DROP_REASONS).

    demoted_accession and sequence_mismatch are tested before no_eligible_instances, because
    seen_families is filled before the universe lookup.
    """
    demoted, eligible = family_stats["demoted"], family_stats["eligible_tier"]
    mismatch = family_stats["mismatch"]
    dropped = []
    for family in sorted(families):
        if family in by_family:
            continue
        if demoted.get(family) and not eligible.get(family):
            dropped.append((family, "demoted_accession"))
        elif mismatch.get(family) and not eligible.get(family):
            dropped.append((family, "sequence_mismatch"))
        elif family in seen_families:
            dropped.append((family, "no_eligible_instances"))
        elif family in dead:
            dropped.append((family, "dead"))
        else:
            dropped.append((family, "not_in_pfam"))
    return dropped


def write_dropped_families(path, dropped):
    """Write the drop report; header-only when nothing dropped (it is a required output)."""
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(DROPPED_FAMILIES_COLUMNS)
        writer.writerows(dropped)


def read_dropped_families(path):
    """Read back what write_dropped_families wrote."""
    with open(path) as fh:
        return [(row["family"], row["reason"]) for row in csv.DictReader(fh, delimiter="\t")]


def read_instances_tsv(path):
    """Read back what write_instances_tsv wrote, as {family: [record]}."""
    by_family = defaultdict(list)
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            by_family[row["family"]].append(row)
    return dict(by_family)


def build_species(by_family):
    """Return (ordered_ids, {id: taxon_id}) for instances and families; a family takes its modal taxon."""
    ids, species = [], {}
    for family in sorted(by_family):
        counts = Counter()
        for rec in by_family[family]:
            taxon = rec["taxon_id"]
            ids.append(rec["instance_id"])
            species[rec["instance_id"]] = taxon
            counts[taxon] += 1
        ids.append(family)
        species[family] = counts.most_common(1)[0][0]
    return ids, species


# ---------------------------------------------------------------------------
# Cache
# ---------------------------------------------------------------------------


class Cache:
    """Release-keyed on-disk cache. A miss falls through to the full stream with identical output."""

    def __init__(self, root, release):
        self.dir = os.path.join(root, f"pfam-{release}") if root else None
        if self.dir:
            os.makedirs(self.dir, exist_ok=True)

    def path(self, name):
        return os.path.join(self.dir, name) if self.dir else None

    def text(self, name, url, gzipped=False):
        """Return a small reference table, reading from cache when present."""
        path = self.path(name)
        if path and os.path.exists(path):
            with open(path) as fh:
                return fh.read()
        text = download_text(url, gzipped=gzipped)
        if path:
            tmp = f"{path}.part"
            with open(tmp, "w") as fh:
                fh.write(text)
            os.replace(tmp, path)
        return text


def content_digest(items):
    """Return a stable digest of a set of ids, for the cache key."""
    return hashlib.sha256("\n".join(sorted(items)).encode()).hexdigest()


def derived_key(release, families, pool_size, seed, universe_digest, tiers, clans_digest):
    """Hash every input that can change the cached instance set or its columns.

    `universe_digest` covers the UniProt flat files (sequence, taxon, review flag), and
    `clans_digest` the requested families' clans, since both are baked into the cache.
    """
    payload = "|".join(
        [
            str(FORMAT_VERSION),
            release,
            str(pool_size),
            str(seed),
            tiers,
            universe_digest,
            clans_digest,
            ",".join(sorted(families)),
        ]
    )
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _sample(items, limit=20):
    """Render a list for a log line, truncated to `limit` items."""
    items = sorted(items)
    if len(items) <= limit:
        return ", ".join(items)
    return f"{', '.join(items[:limit])}, ... and {len(items) - limit} more"


def read_families(path):
    """Return the sorted set of Pfam families listed one per line."""
    with open(path) as fh:
        return sorted({line.strip() for line in fh if line.strip()})


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--families", required=True, help="plain text file, one Pfam family accession per line")
    parser.add_argument("--pool-size", type=int, default=5, help="instances sampled per family (M); default 5")
    parser.add_argument("--seed", type=int, default=42, help="sampling seed; default 42")
    parser.add_argument(
        "--tiers",
        type=parse_tiers,
        required=True,
        metavar="NAME[,NAME...]",
        help="comma-separated strata to sample, out of " + ", ".join(TIERS) + "; unreviewed ones need a "
        "TrEMBL --uniprot-dat",
    )
    parser.add_argument("--clans", help="precomputed Pfam-A.clans.tsv(.gz), skipping that download")
    parser.add_argument("--pfam-regions", help="local Pfam-A.regions.tsv.gz, skipping the 4.7 GB download")
    parser.add_argument(
        "--uniprot-dat",
        action="append",
        default=[],
        metavar="PATH",
        help="local UniProt flat file (.dat.gz); repeatable, Swiss-Prot first. "
        "Swiss-Prot is downloaded when none is given",
    )
    parser.add_argument("--interpro-cache", help="directory for cached downloads and sampled instances")
    parser.add_argument("--pfam-release", help="skip the Pfam.version lookup and use this release string")
    return parser.parse_args()


def main():
    args = parse_args()
    fasta_out = "sequences.fasta"
    species_out = "species.tsv"
    go_out = "go_annotations.tsv"
    instances_out = "instances.tsv"
    dropped_out = "dropped_families.tsv"

    families = read_families(args.families)
    eligible = args.tiers
    tier_names = [TIERS[t] for t in sorted(eligible)]
    print(
        f"Requested {len(families):,} Pfam families (--tiers {','.join(tier_names)})",
        file=sys.stderr,
    )

    if args.pfam_regions and not args.pfam_release:
        print(
            "Warning: --pfam-regions given without --pfam-release; the release is taken from EBI's "
            "current_release, which may not match that copy.",
            file=sys.stderr,
        )
    release = args.pfam_release or parse_release(download_text(PFAM_VERSION_URL, gzipped=True))
    print(f"Pfam release {release}", file=sys.stderr)
    cache = Cache(args.interpro_cache, release)

    if args.clans:
        opener = gzip.open if args.clans.endswith(".gz") else open
        with opener(args.clans, "rt") as fh:
            clans = parse_clans(fh.read())
    else:
        clans = parse_clans(cache.text("Pfam-A.clans.tsv", PFAM_CLANS_URL, gzipped=True))
    # Parsed before the cache key, which digests its accession set.
    taxa_filter = eligible_taxa(eligible)
    dat_paths = list(args.uniprot_dat)
    if not dat_paths:
        cached_dat = cache.path("uniprot_sprot.dat.gz")
        if not cached_dat:
            raise RuntimeError("--uniprot-dat is required when --interpro-cache is not set")
        if not os.path.exists(cached_dat):
            print(f"Caching {SPROT_URL} to {cached_dat} (~600 MB)...", file=sys.stderr)
            download_to_file(SPROT_URL, cached_dat)
        dat_paths = [cached_dat]
        print(
            "No --uniprot-dat given: falling back to Swiss-Prot, which can only fill " f"{TIERS[0]} and {TIERS[1]}",
            file=sys.stderr,
        )
    print(f"Reading the protein universe from {len(dat_paths)} flat file(s)...", file=sys.stderr)
    universe, secondary = load_universe(dat_paths, taxa_filter)
    check_universe_covers(eligible, universe, dat_paths)
    print(
        f"Reference tables: {len(clans):,} families, {len(universe):,} proteins"
        + (f" (taxon {', '.join(sorted(taxa_filter))})" if taxa_filter else " (all taxa)"),
        file=sys.stderr,
    )

    key = derived_key(
        release,
        families,
        args.pool_size,
        args.seed,
        content_digest(universe),
        ",".join(sorted(tier_names)),
        content_digest(f"{f}\t{clans.get(f, '')}" for f in families),
    )
    cached_instances = cache.path(f"instances-{key}.tsv")
    cached_sequences = cache.path(f"sequences-{key}.fasta")
    cached_dropped = cache.path(f"dropped-{key}.tsv")
    stats = Counter()

    # All three or none: the drop report cannot be rebuilt from instances.tsv.
    if cached_instances and all(os.path.exists(p) for p in (cached_instances, cached_sequences, cached_dropped)):
        print(f"Reusing cached instances for key {key}", file=sys.stderr)
        by_family = read_instances_tsv(cached_instances)
        dropped = read_dropped_families(cached_dropped)
        shutil.copyfile(cached_instances, instances_out)
        shutil.copyfile(cached_sequences, fasta_out)
        shutil.copyfile(cached_dropped, dropped_out)
    else:
        wanted = set(families)
        regions_path = args.pfam_regions
        if regions_path:
            print(f"Streaming {regions_path}...", file=sys.stderr)
            with open(regions_path, "rb") as fh:
                by_family, stats, seen, family_stats = sample_instances(
                    iter_gzip_lines(fh), wanted, args.pool_size, universe, args.seed, eligible, secondary
                )
        else:
            print(f"Streaming {PFAM_REGIONS_URL} (~4.7 GB gz)...", file=sys.stderr)
            with _open_url(PFAM_REGIONS_URL, timeout=600) as resp:
                by_family, stats, seen, family_stats = sample_instances(
                    iter_gzip_lines(resp), wanted, args.pool_size, universe, args.seed, eligible, secondary
                )

        write_instances_tsv(instances_out, by_family, clans)
        seqs = {rec["instance_id"]: rec["sequence"] for recs in by_family.values() for rec in recs}
        write_fasta(seqs, seqs.keys(), fasta_out)

        # The dead table is only downloaded when something is absent.
        absent = [f for f in families if f not in by_family and f not in seen]
        dead = parse_dead(cache.text("Pfam-A.dead", PFAM_DEAD_URL, gzipped=True)) if absent else set()
        dropped = classify_dropped(families, by_family, seen, dead, family_stats)
        write_dropped_families(dropped_out, dropped)

        if cached_instances:
            shutil.copyfile(instances_out, cached_instances)
            shutil.copyfile(fasta_out, cached_sequences)
            shutil.copyfile(dropped_out, cached_dropped)

    ids, species = build_species(by_family)
    write_species_tsv(species_out, ids, species)
    # Empty, but consumers still require the file.
    write_go_tsv(go_out, [], {})

    n_instances = sum(len(v) for v in by_family.values())
    print(
        f"Written {n_instances:,} instances across {len(by_family):,} families to {instances_out}",
        file=sys.stderr,
    )
    print(f"Written {n_instances:,} domain sequences to {fasta_out}", file=sys.stderr)
    print(f"Written {len(ids):,} species rows (instances + families) to {species_out}", file=sys.stderr)
    print(f"Written empty GO annotations to {go_out} (dropped in DDI mode)", file=sys.stderr)
    print(f"Written {len(dropped):,} dropped families to {dropped_out}", file=sys.stderr)

    if stats:
        print(f"Scanned {stats['records']:,} Pfam-A region rows", file=sys.stderr)
        for index, tier in enumerate(TIERS):
            offered, kept = stats[f"tier_{tier}_offered"], stats[f"tier_{tier}_kept"]
            gate = "" if index in eligible else "  [not named in --tiers]"
            print(f"  tier {tier}: {offered:,} offered, {kept:,} kept{gate}", file=sys.stderr)
        if stats["parent_cap_skipped"]:
            print(
                f"  {stats['parent_cap_skipped']:,} region rows skipped: parent already at "
                f"its per-family cap of {max(1, args.pool_size // PARENT_CAP_DIVISOR)}",
                file=sys.stderr,
            )
        if stats["not_in_universe"]:
            print(
                f"  {stats['not_in_universe']:,} region rows skipped: parent protein outside "
                "the universe (expected -- Pfam covers every reference proteome)",
                file=sys.stderr,
            )
        if stats["demoted_accession"]:
            print(
                f"  {stats['demoted_accession']:,} region rows skipped: parent is a secondary " "(demoted) accession",
                file=sys.stderr,
            )
        for label, message in (
            ("malformed_rows", "region rows with missing or non-numeric fields"),
            ("coords_out_of_range", "regions whose coordinates fall outside the parent sequence"),
            ("crc64_mismatch", "regions whose Pfam crc64 differs from the UniProt sequence (revised sequence)"),
            ("noncontiguous_blocks", "families whose block was re-entered (file not contiguous)"),
        ):
            if stats[label]:
                print(f"Warning: {stats[label]:,} {message}", file=sys.stderr)

    by_reason = defaultdict(list)
    for family, reason in dropped:
        by_reason[reason].append(family)
    for reason, message in DROP_REASONS.items():
        hit = by_reason.get(reason)
        if hit:
            print(f"Warning: {len(hit)} requested families {message}: {_sample(hit)}", file=sys.stderr)


if __name__ == "__main__":
    main()
