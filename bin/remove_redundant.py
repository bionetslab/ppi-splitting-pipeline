#!/usr/bin/env python3
"""
Remove cross-partition redundant sequences using CD-HIT-2D output.

cd-hit-2d writes to its output file the sequences from db2 (i2) that are NOT
similar to any sequence in db1 (i).  So the IDs present in each output file
are exactly the sequences to KEEP in the smaller partition.

Removal logic:
  val_nr   = val  ∩ {not similar to train}   (from sim_train_val.out)
  test_nr  = test ∩ {not similar to train}   (from sim_train_test.out)
  train_nr = train (unchanged; it is the reference)

PPIs are filtered so both partners must still be present.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import mqc_category, read_fasta, read_ppis, write_fasta, write_ppi_csv


def fasta_ids(path):
    """Return set of protein IDs present in a FASTA file."""
    ids = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.add(line[1:].rstrip().split()[0])
    return ids


def filter_ppis(rows, keep):
    return [row for row in rows if row["protein1"] in keep and row["protein2"] in keep]


def write_mqc(split_results, id_, n_ppis_discarded):
    """The sole contributor to the per-dataset "split_bar_{id_}" PPI
    Partitioning bar chart for kahip/ilp datasets (sort_ppis.py/solve_ilp.py
    contribute nothing themselves -- see their modules): train/val/test show
    their final (post-CD-HIT) Kept size, and the single "discarded" bar is
    stacked by discard reason -- cross-partition (KaHIP/ILP, computed here
    from the original pre-split ppis.csv) vs CD-HIT-2D redundancy removal
    (summed across val+test)."""
    total_cdhit_removed = sum(r["n_ppis_removed"] for r in split_results)

    with open("remove_redundant_bar_mqc.tsv", "w") as fh:
        fh.write(
            f"# id: 'split_bar_{id_}'\n"
            f"# section_name: 'PPI Partitioning: {id_}'\n"
            "# description: 'PPI counts per split. The discarded bar is coloured by why "
            "a PPI never made it into a split: cross-partition (KaHIP/ILP) or removed by "
            "CD-HIT-2D redundancy filtering.'\n"
            "# plot_type: 'bargraph'\n"
            "# pconfig:\n"
            f"#     id: 'split_bar_plot_{id_}'\n"
            f"#     title: 'PPI Partitioning: edges per split ({id_})'\n"
            "#     ylab: '# PPIs'\n"
            "Sample\tKept\tDiscarded (KaHIP/ILP)\tDiscarded (CD-HIT-2D)\n"
        )
        for r in split_results:
            fh.write(f"{mqc_category(r['name'])}\t{r['n_ppis_nr']}\t0\t0\n")
        fh.write(f"{mqc_category('discarded')}\t0\t{n_ppis_discarded}\t{total_cdhit_removed}\n")



def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ppis", required=True,
                    help="Original, pre-split PPI CSV -- used only to compute how many "
                         "PPIs were discarded by KaHIP/ILP for the PPI Partitioning chart")
    ap.add_argument("--train_ppis", required=True)
    ap.add_argument("--val_ppis", required=True)
    ap.add_argument("--test_ppis", required=True)
    ap.add_argument("--train_fasta", required=True)
    ap.add_argument("--val_fasta", required=True)
    ap.add_argument("--test_fasta", required=True)
    ap.add_argument("--sim_train_val", required=True)
    ap.add_argument("--sim_train_test", required=True)
    ap.add_argument("--id", required=True, help="Dataset ID, for MultiQC tagging")
    args = ap.parse_args()

    n_input = len(read_ppis(args.ppis))

    train_seqs = read_fasta(args.train_fasta)
    val_seqs   = read_fasta(args.val_fasta)
    test_seqs  = read_fasta(args.test_fasta)

    train_ppis = read_ppis(args.train_ppis)
    val_ppis   = read_ppis(args.val_ppis)
    test_ppis  = read_ppis(args.test_ppis)

    # IDs to keep (sequences in CD-HIT-2D output = NOT similar to the reference)
    val_keep  = fasta_ids(args.sim_train_val)
    test_keep = fasta_ids(args.sim_train_test)

    train_prot_nr = set(train_seqs)
    val_prot_nr   = set(val_seqs)  & val_keep
    test_prot_nr  = set(test_seqs) & test_keep

    train_ppis_nr = filter_ppis(train_ppis, train_prot_nr)
    val_ppis_nr   = filter_ppis(val_ppis,   val_prot_nr)
    test_ppis_nr  = filter_ppis(test_ppis,  test_prot_nr)

    write_ppi_csv(train_ppis_nr, "train_nr.csv")
    write_ppi_csv(val_ppis_nr,   "val_nr.csv")
    write_ppi_csv(test_ppis_nr,  "test_nr.csv")
    write_fasta(train_seqs, train_prot_nr, "train_nr.fasta")
    write_fasta(val_seqs,   val_prot_nr,   "val_nr.fasta")
    write_fasta(test_seqs,  test_prot_nr,  "test_nr.fasta")

    for name, ppis, prot in [
        ("train_nr", train_ppis_nr, train_prot_nr),
        ("val_nr",   val_ppis_nr,   val_prot_nr),
        ("test_nr",  test_ppis_nr,  test_prot_nr),
    ]:
        print(f"{name}: {len(ppis)} PPIs, {len(prot)} proteins", file=sys.stderr)

    # kahip/ilp assign every PPI to one split or discard it, so comparing input
    # total to the pre-CD-HIT split totals gives the discard count directly,
    # without threading it in from SORT_PPIS/SOLVE_ILP.
    n_ppis_discarded = n_input - (len(train_ppis) + len(val_ppis) + len(test_ppis))

    write_mqc([
        {
            "name": "train",
            "n_ppis_nr": len(train_ppis_nr),
            "n_ppis_removed": 0,
        },
        {
            "name": "val",
            "n_ppis_nr": len(val_ppis_nr),
            "n_ppis_removed": len(val_ppis) - len(val_ppis_nr),
        },
        {
            "name": "test",
            "n_ppis_nr": len(test_ppis_nr),
            "n_ppis_removed": len(test_ppis) - len(test_ppis_nr),
        },
    ], args.id, n_ppis_discarded)


if __name__ == "__main__":
    main()
