#!/usr/bin/env python3
"""
Apply the run-wide dataset order to the collected MultiQC files.

Runs inside MULTIQC, the only task that sees every dataset, so the order does not
enter the hash of the tasks that write the *_mqc.tsv / *_mqc.html files.

Rewrites (zero-padded NN):
  Sample      "<label>_<k>_<split>[_<tail>]" -> "<NN>_<label>_<split>[_<tail>]"
              "<label>"                      -> "<NN>_<label>"
  section id  "<kind>_<label>"               -> "<kind>_<NN>_<label>" (per-dataset kinds only)
and writes a MultiQC config whose report_section_order pins the section order.
"""

import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from utils import SPLIT_ORDER  # noqa: E402

# Section-id prefixes followed by a dataset label.
DATASET_KINDS = [
    "split_bar",
    "bias_scatter",
    "similarity_heatmap",
]

# Report section order: shared summaries first, then each per-dataset kind.
KIND_ORDER = [
    # shared, one section merging every dataset's rows
    "neg_generalstats",
    "pos_neg_bar",
    "classifier_metrics",
    "ddi_examples_bar",
    "ddi_examples_stats",
    "ddi_examples_ilp",
    "neg_ilp_diagnostics",
    "neg_ilp_residuals",
    # per-dataset, one section each
    "split_bar",
    "bias_scatter",
    "similarity_heatmap",
]

# Longest first: "test" would otherwise swallow "test_balanced".
_SPLITS_LONGEST_FIRST = sorted(SPLIT_ORDER, key=len, reverse=True)

_TSV_SECTION_ID = re.compile(r"^# id:\s*(['\"])(?P<id>.*?)\1\s*$")
_HTML_SECTION_ID = re.compile(r"^id:\s*(['\"])(?P<id>.*?)\1\s*$")


class Renamer:
    """Holds the resolved order and every rewrite derived from it."""

    def __init__(self, order):
        self.order = order
        self.n_datasets = len(order)
        self.dataset_width = len(str(self.n_datasets))
        self.split_width = len(str(self.n_datasets * len(SPLIT_ORDER)))
        # Longest first so "minimal_leakage" cannot shadow "minimal_leakage_hcni".
        self.labels_longest_first = sorted(order, key=len, reverse=True)
        self.section_ids = []

    def dataset_index(self, label):
        return f"{self.order.index(label) + 1:0{self.dataset_width}d}"

    def split_index(self, label, split):
        idx = self.order.index(label) * len(SPLIT_ORDER) + SPLIT_ORDER.index(split) + 1
        return f"{idx:0{self.split_width}d}"

    def _split_labels(self, sample):
        """Yield every (label, rest) reading of `sample`, longest label first; several can match."""
        for label in self.labels_longest_first:
            if sample == label:
                yield label, None
            elif sample.startswith(label + "_"):
                yield label, sample[len(label) + 1 :]

    def _rewrite(self, label, rest):
        """Return the renumbered cell for one reading, or None if `rest` is not a split suffix."""
        if rest is None:
            # Dataset-level, from utils.mqc_dataset().
            return f"{self.dataset_index(label)}_{label}"
        # Split-level, from utils.mqc_sample(); "_<tail>" for per-protein rows.
        m = re.match(r"^(\d+)_(?P<tail>.+)$", rest)
        if not m:
            return None
        tail = m.group("tail")
        for split in _SPLITS_LONGEST_FIRST:
            if tail == split:
                return f"{self.split_index(label, split)}_{label}_{split}"
            if tail.startswith(split + "_"):
                return f"{self.split_index(label, split)}_{label}_{tail}"
        return None

    def sample(self, value):
        """Rewrite one Sample cell; utils.mqc_category() values carry no label and are returned unchanged."""
        for label, rest in self._split_labels(value):
            out = self._rewrite(label, rest)
            if out is not None:
                return out
        return value

    def section_id(self, sid):
        """Insert the dataset index into a per-dataset section id."""
        kinds = sorted((k for k in DATASET_KINDS if sid.startswith(k + "_")), key=len, reverse=True)
        for kind in kinds:
            label = sid[len(kind) + 1 :]
            return f"{kind}_{self.dataset_index(label)}_{label}"
        return sid


def kind_of(sid, dataset_width):
    """The section id with its "_<NN>_<label>" tail removed, for KIND_ORDER."""
    for kind in sorted(DATASET_KINDS, key=len, reverse=True):
        prefix = kind + "_"
        if sid.startswith(prefix):
            rest = sid[len(prefix) :]
            if re.match(rf"^\d{{{dataset_width}}}_", rest):
                return kind
    for kind in sorted(KIND_ORDER, key=len, reverse=True):
        if sid == kind or sid.startswith(kind + "_"):
            return kind


def rewrite_tsv(lines, renamer):
    """Rewrite the section id and the Sample column of a custom-content TSV."""
    out = []
    seen_header_row = False
    for line in lines:
        if line.startswith("#"):
            m = _TSV_SECTION_ID.match(line.rstrip("\n"))
            if m:
                sid = m.group("id")
                new = renamer.section_id(sid)
                renamer.section_ids.append(new)
                out.append(f"# id: '{new}'\n")
                continue
            out.append(line)
            continue
        if not seen_header_row:
            seen_header_row = True
            out.append(line)
            continue
        fields = line.rstrip("\n").split("\t")
        fields[0] = renamer.sample(fields[0])
        out.append("\t".join(fields) + "\n")
    return out


def rewrite_html(lines, renamer):
    """Rewrite the leading HTML-comment header of a custom-content HTML file."""
    out = []
    in_header = False
    done = False
    for line in lines:
        stripped = line.strip()
        if not done:
            if stripped.startswith("<!--"):
                in_header = True
            if in_header:
                m = _HTML_SECTION_ID.match(stripped)
                if m:
                    sid = m.group("id")
                    new = renamer.section_id(sid)
                    renamer.section_ids.append(new)
                    out.append(f"id: '{new}'\n")
                    continue
                if stripped.endswith("-->"):
                    in_header = False
                    done = True
        out.append(line)
    return out


def write_config(path, section_ids, dataset_width):
    """Write report_section_order listing every section id; higher order renders first."""
    by_kind = {}
    for sid in sorted(set(section_ids)):
        by_kind.setdefault(kind_of(sid, dataset_width), []).append(sid)

    ranked = []
    for kind in KIND_ORDER:
        # Ids carry the zero-padded dataset index, so a plain sort is the resolved order.
        ranked.extend(sorted(by_kind.get(kind, [])))

    # Base well above MultiQC's defaults (10, 20, ...) so no unlisted section ties.
    step = 10
    base = 100000
    with open(path, "w") as fh:
        fh.write("# Generated per run by bin/relabel_mqc.py -- do not edit.\n")
        fh.write("# Higher order renders nearer the top. Shared cross-dataset summaries\n")
        fh.write("# first, then each per-dataset kind with the datasets in the resolved\n")
        fh.write("# --mqc-order.\n")
        fh.write("report_section_order:\n")
        fh.write(f"  general_stats:\n    order: {base + step}\n")
        for i, sid in enumerate(ranked):
            fh.write(f"  {sid}:\n    order: {base - i * step}\n")
    return ranked


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("source", help="directory holding the collected *_mqc.tsv / *_mqc.html files")
    ap.add_argument("--out", required=True, help="directory to write the rewritten tree into")
    ap.add_argument(
        "--mqc-order",
        default="",
        help="comma-joined display labels, in report order (main.nf's resolved order)",
    )
    ap.add_argument("--config", required=True, help="MultiQC config to write (report_section_order)")
    args = ap.parse_args()

    order = [s.strip() for s in args.mqc_order.split(",") if s.strip()]
    renamer = Renamer(order)

    out_root = os.path.abspath(args.out)
    src_root = os.path.abspath(args.source)
    os.makedirs(out_root, exist_ok=True)

    n_files = 0
    for dirpath, dirnames, filenames in os.walk(src_root):
        # Never descend into our own output, and skip Nextflow's task files.
        dirnames[:] = [d for d in dirnames if os.path.join(dirpath, d) != out_root and not d.startswith(".")]
        for name in sorted(filenames):
            if not (name.endswith("_mqc.tsv") or name.endswith("_mqc.html")):
                continue
            src = os.path.join(dirpath, name)
            rel = os.path.relpath(src, src_root)
            dst = os.path.join(out_root, rel)
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            with open(src) as fh:
                lines = fh.readlines()
            if name.endswith("_mqc.tsv"):
                lines = rewrite_tsv(lines, renamer)
            else:
                lines = rewrite_html(lines, renamer)
            with open(dst, "w") as fh:
                fh.writelines(lines)
            n_files += 1

    ranked = write_config(args.config, renamer.section_ids, renamer.dataset_width)

    print(
        f"Relabelled {n_files} MultiQC file(s), {len(ranked)} section(s), "
        f"{len(order)} dataset(s): {', '.join(order)}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
