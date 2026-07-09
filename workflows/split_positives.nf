include { SORT_PPIS; SOLVE_ILP; CDHIT2D; REMOVE_REDUNDANT } from '../modules/splitting'

// Assigns PPIs to train/val/test (via KaHIP-based sorting or the ILP
// splitter), then removes cross-split redundancy with CD-HIT-2D.
workflow SPLIT_POSITIVES {
    take:
    ppis_ch
    sequences_ch
    partition_ch
    node_mapping_ch

    main:
    if (params.split_method == "ilp") {
        gurobi_license_ch = params.gurobi_license
            ? channel.value(file(params.gurobi_license, checkIfExists: true))
            : channel.value([])
        sorted = SOLVE_ILP(ppis_ch, sequences_ch, partition_ch, node_mapping_ch, gurobi_license_ch)
    } else {
        sorted = SORT_PPIS(ppis_ch, partition_ch, sequences_ch, node_mapping_ch)
    }

    // One channel, one process: each (label, fasta1, fasta2) item becomes its
    // own task, so Nextflow runs both CD-HIT-2D comparisons in parallel.
    cdhit_inputs_ch = sorted.train_fasta.combine(sorted.val_fasta).map { t, v -> tuple("train_val", t, v) }
        .mix(sorted.train_fasta.combine(sorted.test_fasta).map { t, te -> tuple("train_test", t, te) })

    cdhit_out = CDHIT2D(cdhit_inputs_ch)

    cdhit_branched = cdhit_out.sim.branch {
        label, f ->
            train_val:  label == "train_val"
            train_test: label == "train_test"
    }
    sim_tv = cdhit_branched.train_val.map { label, f -> f }
    sim_tt = cdhit_branched.train_test.map { label, f -> f }

    nr = REMOVE_REDUNDANT(
        sorted.train_ppis,
        sorted.val_ppis,
        sorted.test_ppis,
        sorted.train_fasta,
        sorted.val_fasta,
        sorted.test_fasta,
        sim_tv,
        sim_tt
    )

    emit:
    train_ppis  = nr.train_ppis
    val_ppis    = nr.val_ppis
    test_ppis   = nr.test_ppis
    train_fasta = nr.train_fasta
    val_fasta   = nr.val_fasta
    test_fasta  = nr.test_fasta
    sorted_mqc  = sorted.mqc
    nr_mqc      = nr.mqc
}