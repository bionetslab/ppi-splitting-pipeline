#!/usr/bin/env python3
"""Retrain the RF baseline on the ESM-2 robustness splits with other embeddings (e.g. ProtT5).

The ESM-2 robustness runs keep every (split, sampling) seed combination's splits. This reuses one combination
(default split seed 1234, sampling seed 42, training seed 42, the pipeline defaults) so an embedding comparison
differs only in the embeddings. Training and evaluation are train_classifier.py's own.

Writes <outdir>/<id>/classifier_metrics_test_{balanced,realistic}_mqc.tsv (as the pipeline does) and
<outdir>/multiqc_classifier_metrics_test_{balanced,realistic}_table.txt with all IDs (as MultiQC does).

Run it for both embeddings in the same environment, so only the embeddings differ (scikit-learn isn't pinned, so the
pipeline's own ESM-2 numbers can differ slightly):
    retrain_on_esm2_splits.py --embeddings data/all_embeddings_esm2.npz data/structural_embeddings_esm2.npz \\
        --results results_robustness results_robustness_hcn --outdir results_esm2_same_splits
    retrain_on_esm2_splits.py --embeddings data/embeddings_prot_t5_esm2_proteins.npz \\
        --results results_robustness results_robustness_hcn --outdir results_prott5_same_splits
data/embeddings_prot_t5_esm2_proteins.npz holds only the proteins that also have an ESM-2 embedding, so the same
pairs are dropped as with ESM-2 (build_X skips pairs with a missing embedding).
"""

import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from train_classifier import train_and_evaluate, write_mqc  # noqa: E402
from utils import load_embeddings  # noqa: E402

COLS = ["auroc", "auprc", "f1", "mcc", "precision", "recall", "accuracy"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--embeddings", required=True, nargs="+",
                    help=".npz file(s) to train with; several are merged (e.g. the two ESM-2 files)")
    ap.add_argument("--results", required=True, nargs="+", type=Path, help="Robustness --outdir(s)")
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--split_seed", type=int, default=1234)
    ap.add_argument("--sampling_seed", type=int, default=42)
    ap.add_argument("--training_seed", type=int, default=42)
    ap.add_argument("--ids", nargs="*", help="Only these IDs (default: every ID with that seed combination)")
    args = ap.parse_args()

    runs = sorted(
        (d.parts[-4], d)
        for root in args.results
        for d in root.glob(f"*/robustness/split{args.split_seed}/sampling{args.sampling_seed}")
        if args.ids is None or d.parts[-4] in args.ids
    )
    if not runs:
        sys.exit("No splits found for that seed combination.")

    embeddings = {}
    for f in args.embeddings:
        print(f"Loading embeddings from {f} ...", file=sys.stderr)
        embeddings.update(load_embeddings(f))

    rows = {"test_balanced": [], "test_realistic": []}
    for id_, d in runs:
        print(f"=== {id_} ===", file=sys.stderr)
        csvs = [d / f"{s}.csv" for s in ("train", "val", "test_balanced", "test_realistic")]
        # build_X drops pairs with a missing embedding, as in the pipeline; report them since they change the data
        missing = {p for f in csvs for line in f.read_text().splitlines()[1:] for p in line.split(",")[:2]} - embeddings.keys()
        if missing:
            print(f"  {len(missing)} proteins have no embedding, their pairs are dropped: {sorted(missing)[:5]}",
                  file=sys.stderr)
        results = train_and_evaluate(*csvs, embeddings, args.training_seed)

        out = args.outdir / id_
        out.mkdir(parents=True, exist_ok=True)
        cwd = os.getcwd()
        os.chdir(out)  # write_mqc writes into the working directory, like in the pipeline
        try:
            write_mqc(results, id_)
        finally:
            os.chdir(cwd)
        for name, m in results:
            rows[name].append("\t".join([id_, id_] + [f"{m[c]:.4f}" for c in COLS]))

    for name, lines in rows.items():
        with open(args.outdir / f"multiqc_classifier_metrics_{name}_table.txt", "w") as fh:
            fh.write("Sample\tID\tAUROC\tAUPRC\tF1\tMCC\tPrecision\tRecall\tAccuracy\n")
            fh.write("\n".join(lines) + "\n")
    print(f"Wrote {len(runs)} IDs to {args.outdir}", file=sys.stderr)


if __name__ == "__main__":
    main()
