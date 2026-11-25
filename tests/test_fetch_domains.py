#!/usr/bin/env python3
"""Tests for bin/fetch_domains.py, on local fixtures with no network.

Run directly (python3 tests/test_fetch_domains.py) or via pytest.
"""

import csv
import gzip
import io
import os
import subprocess
import sys
import tempfile
from collections import Counter

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(REPO, "bin")
sys.path.insert(0, BIN)

FETCH = os.path.join(BIN, "fetch_domains.py")

from fetch_domains import (  # noqa: E402
    INSTANCE_COLUMNS,
    TIERS,
    check_universe_covers,
    crc64,
    eligible_taxa,
    load_universe,
    sample_instances,
    parse_tiers,
    parse_uniprot_dat,
    tier_of,
)

ENV = dict(os.environ, PYTHONPATH=BIN + os.pathsep + os.environ.get("PYTHONPATH", ""))

ALL_TIERS = ",".join(TIERS)

REGIONS_HEADER = [
    "pfamseq_acc",
    "seq_version",
    "crc64",
    "md5",
    "pfamA_acc",
    "seq_start",
    "seq_end",
    "ali_start",
    "ali_end",
]

# PF00001: two parents, one with two regions. PF00002: human + mouse. PF00003: mouse only.
# PF00006: mouse-reviewed vs human-unreviewed. PF00007: human unreviewed. PF00008: non-human
# unreviewed. PF00009: parent not in the DAT. PF00010: only a secondary accession.
REGIONS = (
    [
        ("P11111", "PF00001", 10, 20, 11, 19),
        ("P11111", "PF00001", 50, 62, 51, 60),
        ("P22222", "PF00001", 5, 15, 6, 14),
        ("P33333", "PF00002", 1, 12, 2, 11),
        ("Q44444", "PF00002", 3, 14, 4, 13),
        ("Q55555", "PF00003", 7, 17, 8, 16),
        ("Q66666", "PF00006", 1, 12, 1, 12),
        ("U77777", "PF00006", 2, 14, 2, 13),
        ("U77777", "PF00007", 3, 15, 3, 14),
        ("U88888", "PF00008", 1, 12, 1, 11),
        ("P99999", "PF00009", 1, 10, 1, 10),  # P99999 absent from the DAT
        # Secondary accessions of P11111; PF00001 must be unaffected.
        ("S11111", "PF00001", 30, 40, 31, 39),
        ("S11112", "PF00001", 60, 70, 61, 69),
        ("S33333", "PF00010", 1, 12, 2, 11),
        ("P66666", "PF00004", 900, 999, 901, 990),  # coordinates past the sequence end
        # PF00011: coordinates fit, but the crc64 does not (revised sequence).
        ("P12121", "PF00011", 1, 10, 1, 10),
    ]
    + [
        # PF00005: repeat domain, ten copies on P77777 and one on P88888 (per-parent cap).
        ("P77777", "PF00005", i * 3 + 1, i * 3 + 3, i * 3 + 1, i * 3 + 3)
        for i in range(10)
    ]
    + [
        ("P88888", "PF00005", 1, 6, 1, 6),
    ]
)

# accession -> (taxon, reviewed, sequence); every stratum has a member.
PROTEINS = {
    "P11111": ("9606", True, "".join(f"{i % 10}" for i in range(80)).replace("0", "A")),
    "P22222": ("9606", True, "M" * 40),
    "P33333": ("9606", True, "K" * 40),
    "Q44444": ("10090", True, "L" * 40),
    "Q55555": ("10090", True, "V" * 40),
    "Q66666": ("10090", True, "C" * 40),
    "U77777": ("9606", False, "D" * 40),
    "U88888": ("10090", False, "E" * 40),
    "P66666": ("9606", True, "W" * 40),
    "P77777": ("9606", True, "R" * 40),
    "P88888": ("9606", True, "S" * 40),
    "P12121": ("9606", True, "T" * 40),
}

# Parents whose Pfam crc64 is written wrong.
CRC_MISMATCH = {"P12121"}

# accession -> its secondary accessions.
SECONDARY = {
    "P11111": ["S11111", "S11112"],
    "P33333": ["S33333"],
}

FAMILIES = [
    "PF00001",
    "PF00002",
    "PF00003",
    "PF00004",
    "PF00005",
    "PF00006",
    "PF00007",
    "PF00008",
    "PF00009",
    "PF00010",
    "PF00011",
]

HUMAN_REVIEWED = {acc for acc, (taxon, reviewed, _s) in PROTEINS.items() if taxon == "9606" and reviewed}


def dat_entry(acc, taxon, reviewed, seq, secondaries=()):
    out = io.StringIO()
    status = "Reviewed;" if reviewed else "Unreviewed;"
    out.write(f"ID   {acc}_TEST   {status}   {len(seq)} AA.\n")
    # First secondary on the primary's AC line, the rest on a later one.
    first = "; ".join([acc] + list(secondaries[:1]))
    out.write(f"AC   {first};\n")
    if len(secondaries) > 1:
        out.write("AC   " + "; ".join(secondaries[1:]) + ";\n")
    out.write(f"OX   NCBI_TaxID={taxon};\n")
    out.write(f"SQ   SEQUENCE   {len(seq)} AA;  0 MW;  0 CRC64;\n")
    for i in range(0, len(seq), 60):
        out.write(f"     {seq[i:i + 60]}\n")
    out.write("//\n")
    return out.getvalue()


def write_dat(path, entries, secondary=None):
    secondary = SECONDARY if secondary is None else secondary
    with gzip.open(path, "wt") as fh:
        for acc, (taxon, reviewed, seq) in entries.items():
            fh.write(dat_entry(acc, taxon, reviewed, seq, secondary.get(acc, ())))
    return path


def write_fixtures(tmp):
    regions = os.path.join(tmp, "Pfam-A.regions.tsv.gz")
    with gzip.open(regions, "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(REGIONS_HEADER)
        for acc, fam, s_start, s_end, a_start, a_end in REGIONS:
            crc = crc64(PROTEINS[acc][2]) if acc in PROTEINS else "CRC"
            if acc in CRC_MISMATCH:
                crc = "0000000000000000"
            w.writerow([acc, "1", crc, "MD5", fam, s_start, s_end, a_start, a_end])

    dat = write_dat(os.path.join(tmp, "uniprot.dat.gz"), PROTEINS)

    clans = os.path.join(tmp, "Pfam-A.clans.tsv.gz")
    with gzip.open(clans, "wt") as fh:
        fh.write("PF00001\tCL0001\tname\tshort\tdesc\n")  # PF00002+ have no clan row

    families = os.path.join(tmp, "families.txt")
    with open(families, "w") as fh:
        fh.write("\n".join(FAMILIES) + "\n")
    return regions, dat, clans, families


def slug(tiers):
    return tiers.replace(",", "+")


def run(tmp, tiers, pool_size=25, dats=None, check=True, out_name=None):
    regions, dat, clans, families = write_fixtures(tmp)
    out = os.path.join(tmp, out_name or slug(tiers))
    os.makedirs(out, exist_ok=True)
    cmd = [
        sys.executable,
        FETCH,
        "--families",
        families,
        "--pool-size",
        str(pool_size),
        "--seed",
        "42",
        "--tiers",
        tiers,
        "--clans",
        clans,
        "--pfam-regions",
        regions,
        "--pfam-release",
        "38.2",
    ]
    for path in (dats if dats is not None else [dat]):
        cmd += ["--uniprot-dat", path]
    proc = subprocess.run(cmd, cwd=out, env=ENV, capture_output=True, text=True)
    if check and proc.returncode != 0:
        raise AssertionError(f"fetch_domains.py failed:\n{proc.stderr}")
    if check:
        return out, proc.stderr
    return out, proc


def read_instances(out):
    with open(os.path.join(out, "instances.tsv"), newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_dropped(out):
    with open(os.path.join(out, "dropped_families.tsv"), newline="") as fh:
        return {r["family"]: r["reason"] for r in csv.DictReader(fh, delimiter="\t")}


def read_fasta(out):
    seqs, key = {}, None
    with open(os.path.join(out, "sequences.fasta")) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                key = line[1:]
            elif key:
                seqs[key] = seqs.get(key, "") + line
    return seqs


def test_unit_helpers():
    assert eligible_taxa(parse_tiers("human_reviewed")) == {"9606"}
    assert eligible_taxa(parse_tiers("human_reviewed,human_unreviewed")) == {"9606"}
    assert eligible_taxa(parse_tiers("human_reviewed,other_reviewed")) is None
    assert eligible_taxa(parse_tiers(ALL_TIERS)) is None

    # UniProt's own SQ-line checksum for Q6GZX4 (the first Swiss-Prot entry).
    assert (
        crc64(
            "MAFSAEDVLKEYDRRRRMEALLLSLYYPNDRKLLDYKEWSPPRVQVECPKAPVEWNNPPSEKGLIVGHFSGIKYKGEKAQASEVDVNKMCCWVSKFKDAMRRY"
            "QGIQTCKIPGKVLSDLDAKIKAYNLTVEGVEGFVRYSRVTKQHVAAFLKELRHSKQYENVNLIHYILTDKRVDIQHLEKDLVKDFKALVESAHRMRQGHMINV"
            "KYILYQLLKKHGHGPDGPDILTVKTGSKGVLYDDSFRKIYTDLGWKFTPL"
        )
        == "B4840739BF7D4121"
    )

    assert parse_tiers("human_reviewed") == frozenset({0})
    assert parse_tiers("other_unreviewed, human_reviewed") == frozenset({0, 3})
    for bad in ("", "human_only", "any", "0,2"):
        try:
            parse_tiers(bad)
        except Exception as exc:
            assert "valid names are" in str(exc) or "name at least one of" in str(exc), str(exc)
        else:
            raise AssertionError(f"--tiers {bad!r} should have been rejected")

    # Reviewed outranks human.
    assert tier_of("9606", True) == 0
    assert tier_of("10090", True) == 1
    assert tier_of("9606", False) == 2
    assert tier_of("10090", False) == 3

    with tempfile.TemporaryDirectory() as tmp:
        _regions, dat, _clans, _families = write_fixtures(tmp)
        with gzip.open(dat, "rb") as fh:
            everything = parse_uniprot_dat(fh)
        assert set(everything) == set(PROTEINS), sorted(everything)
        assert everything["P11111"] == ("9606", PROTEINS["P11111"][2], True)
        assert everything["U77777"][2] is False, everything["U77777"]
        # Primary accessions only; secondaries come back as a separate set.
        seen_secondary = set()
        with gzip.open(dat, "rb") as fh:
            primaries = parse_uniprot_dat(fh, None, None, seen_secondary)
        assert set(primaries) == set(PROTEINS), sorted(primaries)
        assert seen_secondary == {"S11111", "S11112", "S33333"}, sorted(seen_secondary)
        assert not (seen_secondary & set(primaries)), sorted(seen_secondary & set(primaries))

        # The taxon filter applies to the secondaries too.
        mouse_only = set()
        with gzip.open(dat, "rb") as fh:
            parse_uniprot_dat(fh, {"10090"}, None, mouse_only)
        assert mouse_only == set(), sorted(mouse_only)

        with gzip.open(dat, "rb") as fh:
            human = parse_uniprot_dat(fh, {"9606"})
        assert set(human) == {a for a, (t, _r, _s) in PROTEINS.items() if t == "9606"}, sorted(human)
        assert "Q44444" not in human and "Q55555" not in human, sorted(human)

        # An ID line with neither token is fatal.
        broken = "ID   X99999_TEST   393 AA.\nAC   X99999;\nOX   NCBI_TaxID=9606;\nSQ   SEQUENCE\n     MMMM\n//\n"
        try:
            parse_uniprot_dat(io.StringIO(broken))
        except RuntimeError as exc:
            assert "Reviewed;" in str(exc), str(exc)
        else:
            raise AssertionError("an ID line without a review token should be fatal")

        # First writer wins.
        merged = {}
        with gzip.open(dat, "rb") as fh:
            parse_uniprot_dat(fh, None, merged)
        promoted = write_dat(os.path.join(tmp, "promoted.dat.gz"), {"U77777": ("9606", True, "D" * 40)})
        with gzip.open(promoted, "rb") as fh:
            parse_uniprot_dat(fh, None, merged)
        assert merged["U77777"][2] is False, merged["U77777"]

        # A file that contributes nothing is fatal.
        try:
            load_universe([dat, dat])
        except RuntimeError as exc:
            assert "contributed no proteins" in str(exc), str(exc)
        else:
            raise AssertionError("a zero-contribution --uniprot-dat should be fatal")

        # An eligible stratum the universe cannot hold is fatal.
        reviewed_only = {a: (t, s, r) for a, (t, r, s) in PROTEINS.items() if r}
        check_universe_covers(parse_tiers("human_reviewed,other_reviewed"), reviewed_only, ["x.dat.gz"])
        for tiers, word in (("human_reviewed,human_unreviewed", "unreviewed"), ("other_reviewed", "non-human")):
            universe = (
                reviewed_only
                if word == "unreviewed"
                else {a: (t, s, r) for a, (t, r, s) in PROTEINS.items() if t == "9606"}
            )
            try:
                check_universe_covers(parse_tiers(tiers), universe, ["x.dat.gz"])
            except RuntimeError as exc:
                assert word in str(exc), str(exc)
            else:
                raise AssertionError(f"--tiers {tiers} should not survive that universe")


def test_formats_and_coordinates():
    with tempfile.TemporaryDirectory() as tmp:
        out, log = run(tmp, ALL_TIERS)

        rows = read_instances(out)
        with open(os.path.join(out, "instances.tsv"), newline="") as fh:
            header = fh.readline().rstrip("\r\n").split("\t")
        assert header == list(INSTANCE_COLUMNS), header

        by_id = {r["instance_id"]: r for r in rows}
        # instance_id is {family}_{accession}_{ali_start}_{ali_end} ...
        assert "PF00001_P11111_11_19" in by_id, sorted(by_id)
        # ... and the sequence is the 1-based inclusive slice of the parent.
        parent = PROTEINS["P11111"][2]
        seqs = read_fasta(out)
        assert seqs["PF00001_P11111_11_19"] == parent[10:19]
        assert len(seqs["PF00001_P11111_11_19"]) == 9

        rec = by_id["PF00001_P11111_11_19"]
        assert (rec["family"], rec["protein_id"], rec["start"], rec["end"]) == ("PF00001", "P11111", "11", "19")
        assert rec["taxon_id"] == "9606"
        assert rec["source_db"] == "reviewed"
        assert rec["clan"] == "CL0001"
        # A family with no clan row gets the synthetic singleton clan.
        assert by_id["PF00002_P33333_2_11"]["clan"] == "CL_PF00002"

        # source_db is the record's own review flag.
        assert by_id["PF00007_U77777_3_14"]["source_db"] == "unreviewed"
        assert {r["source_db"] for r in rows} == {"reviewed", "unreviewed"}

        # The mouse-only family survives the full tier set, via other_reviewed.
        assert any(r["family"] == "PF00003" for r in rows), sorted(by_id)
        # ... and the non-human unreviewed one via the bottom stratum.
        assert any(r["family"] == "PF00008" for r in rows), sorted(by_id)

        # Every FASTA key is an instance id and every instance has a sequence.
        assert set(seqs) == set(by_id), sorted(set(seqs) ^ set(by_id))

        # species.tsv: instance rows plus one row per family, whose taxon is the
        # modal taxon of its instances.
        with open(os.path.join(out, "species.tsv"), newline="") as fh:
            species = {r["protein_id"]: r["taxon_id"] for r in csv.DictReader(fh, delimiter="\t")}
        assert species["PF00001"] == "9606"
        assert species["PF00003"] == "10090"
        assert species["PF00001_P11111_11_19"] == "9606"

        # PF00004: out of range, PF00011: wrong CRC64, PF00009: parent outside the universe.
        dropped = read_dropped(out)
        assert set(dropped) == {"PF00004", "PF00009", "PF00010", "PF00011"}, dropped
        assert dropped["PF00004"] == "sequence_mismatch", dropped
        assert dropped["PF00011"] == "sequence_mismatch", dropped
        assert dropped["PF00009"] == "no_eligible_instances", dropped
        assert dropped["PF00010"] == "demoted_accession", dropped

        # No secondary accession is ever written.
        secondaries = {s for accs in SECONDARY.values() for s in accs}
        assert not ({r["protein_id"] for r in rows} & secondaries), rows
        assert not any(s in i for i in by_id for s in secondaries), sorted(by_id)
        assert {r["protein_id"] for r in rows if r["family"] == "PF00001"} == {"P11111", "P22222"}

        # go_annotations.tsv is header-only.
        with open(os.path.join(out, "go_annotations.tsv")) as fh:
            go = fh.read().splitlines()
        assert len(go) == 1 and go[0].split("\t")[0] == "protein_id", go

        assert "region rows skipped" in log, log
        assert "coordinates fall outside" in log, log
        assert "revised sequence" in log, log
        assert "secondary" in log and "demoted" in log, log
        assert "primary accession each" in log, log


def test_reviewed_outranks_human_when_filling():
    with tempfile.TemporaryDirectory() as tmp:
        # PF00006 with room for one: the mouse-reviewed candidate beats the human-unreviewed one.
        out, _ = run(tmp, ALL_TIERS, pool_size=1)
        picked = [r for r in read_instances(out) if r["family"] == "PF00006"]
        assert [r["protein_id"] for r in picked] == ["Q66666"], picked
        assert picked[0]["source_db"] == "reviewed"

        # With room for two, tier 2 still contributes.
        out, _ = run(tmp, ALL_TIERS, pool_size=2, out_name="pool2")
        parents = sorted(r["protein_id"] for r in read_instances(out) if r["family"] == "PF00006")
        assert parents == ["Q66666", "U77777"], parents


def test_human_reviewed_is_a_subset_and_keeps_the_same_instances():
    with tempfile.TemporaryDirectory() as tmp:
        all_out, _ = run(tmp, ALL_TIERS)
        human_out, human_log = run(tmp, "human_reviewed")

        all_rows = read_instances(all_out)
        all_ids = {r["instance_id"] for r in all_rows}
        human_rows = read_instances(human_out)
        human_ids = {r["instance_id"] for r in human_rows}

        assert human_ids < all_ids, sorted(human_ids - all_ids)
        # The human-reviewed instances are identical between the two tier sets.
        assert human_ids == {
            r["instance_id"] for r in all_rows if r["taxon_id"] == "9606" and r["source_db"] == "reviewed"
        }
        assert {r["protein_id"] for r in human_rows} <= HUMAN_REVIEWED
        assert all(r["source_db"] == "reviewed" for r in human_rows)

        # A family with nothing in a named stratum is dropped, not filled from another.
        dropped = read_dropped(human_out)
        for family in ("PF00003", "PF00006", "PF00007", "PF00008"):
            assert dropped[family] == "no_eligible_instances", dropped
        assert dropped["PF00010"] == "demoted_accession", dropped
        assert "9606" in human_log, human_log


def test_several_flat_files_merge_first_writer_first():
    with tempfile.TemporaryDirectory() as tmp:
        _regions, dat, _clans, _families = write_fixtures(tmp)
        # U77777 promoted to reviewed: listed first it wins; a file adding nothing is fatal.
        promoted = write_dat(
            os.path.join(tmp, "promoted.dat.gz"),
            {"U77777": ("9606", True, "D" * 40), "X12345": ("9606", True, "Y" * 40)},
        )

        out, _ = run(tmp, "human_reviewed", dats=[promoted, dat], out_name="promoted_first")
        rows = {r["family"]: r for r in read_instances(out)}
        assert "PF00007" in rows, sorted(rows)
        assert rows["PF00007"]["source_db"] == "reviewed", rows["PF00007"]

        _out, proc = run(tmp, "human_reviewed", dats=[dat, dat], check=False, out_name="dupe")
        assert proc.returncode != 0
        assert "contributed no proteins" in proc.stderr, proc.stderr


def test_a_tier_set_the_universe_cannot_fill_is_fatal():
    with tempfile.TemporaryDirectory() as tmp:
        reviewed_only = write_dat(
            os.path.join(tmp, "reviewed_only.dat.gz"),
            {a: v for a, v in PROTEINS.items() if v[1]},
        )
        human_only = write_dat(
            os.path.join(tmp, "human_only.dat.gz"),
            {a: v for a, v in PROTEINS.items() if v[0] == "9606"},
        )

        _out, proc = run(tmp, "human_reviewed,human_unreviewed", dats=[reviewed_only], check=False, out_name="g1")
        assert proc.returncode != 0
        assert "human_unreviewed" in proc.stderr and "are unreviewed" in proc.stderr, proc.stderr

        _out, proc = run(tmp, "human_reviewed,other_reviewed", dats=[human_only], check=False, out_name="g2")
        assert proc.returncode != 0
        assert "other_reviewed" in proc.stderr and "are non-human" in proc.stderr, proc.stderr

        _out, proc = run(tmp, "human_only", check=False, out_name="g3")
        assert proc.returncode != 0
        assert "valid names are" in proc.stderr and "other_unreviewed" in proc.stderr, proc.stderr


def test_filling_prefers_distinct_parent_proteins():
    with tempfile.TemporaryDirectory() as tmp:
        # Room for two of PF00001's three regions: one from each parent.
        out, _ = run(tmp, ALL_TIERS, pool_size=2)
        parents = [r["protein_id"] for r in read_instances(out) if r["family"] == "PF00001"]
        assert len(parents) == 2, parents
        assert len(set(parents)) == 2, parents


def test_a_repeat_domain_family_is_not_filled_from_one_protein():
    with tempfile.TemporaryDirectory() as tmp:
        # pool_size 25 -> cap 5. P77777 carries ten PF00005 regions, P88888 one.
        out, log = run(tmp, ALL_TIERS, pool_size=25)
        rows = [r for r in read_instances(out) if r["family"] == "PF00005"]
        parents = Counter(r["protein_id"] for r in rows)
        assert parents["P77777"] == 5, parents
        assert parents["P88888"] == 1, parents
        assert "per-family cap of 5" in log, log

        # pool_size 5 -> cap 1.
        out, _ = run(tmp, ALL_TIERS, pool_size=5, out_name="cap5")
        rows = [r for r in read_instances(out) if r["family"] == "PF00005"]
        parents = Counter(r["protein_id"] for r in rows)
        assert parents == {"P77777": 1, "P88888": 1}, parents


def test_a_rejected_region_never_evicts_a_sampled_one():
    # Pool size 1: a good region, then one whose coordinates run past its parent.
    stream = [
        "\t".join(["pfamA_acc", "pfamseq_acc", "ali_start", "ali_end"]),
        "\t".join(["PF09999", "A00001", "1", "10"]),
        "\t".join(["PF09999", "B00002", "5", "500"]),
    ]
    universe = {"A00001": ("9606", "M" * 40, True), "B00002": ("9606", "K" * 40, True)}
    for seed in range(40):
        by_family, stats, _seen, _fs = sample_instances(iter(stream), {"PF09999"}, 1, universe, seed)
        assert [r["protein_id"] for r in by_family.get("PF09999", [])] == ["A00001"], (seed, by_family)
        assert stats["coords_out_of_range"] == 1


if __name__ == "__main__":
    test_unit_helpers()
    test_a_rejected_region_never_evicts_a_sampled_one()
    test_formats_and_coordinates()
    test_reviewed_outranks_human_when_filling()
    test_human_reviewed_is_a_subset_and_keeps_the_same_instances()
    test_several_flat_files_merge_first_writer_first()
    test_a_tier_set_the_universe_cannot_fill_is_fatal()
    test_filling_prefers_distinct_parent_proteins()
    test_a_repeat_domain_family_is_not_filled_from_one_protein()
    print("ok: regions-based sampling satisfies the output contract")
