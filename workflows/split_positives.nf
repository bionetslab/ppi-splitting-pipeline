include { SORT_PPIS; SOLVE_ILP; CDHIT2D; REMOVE_REDUNDANT } from '../processes/splitting'

// Assigns PPIs to train/val/test (via KaHIP-based sorting or the ILP
// splitter, chosen per-dataset via meta.split_method), then removes
// cross-split redundancy with CD-HIT-2D.
workflow SPLIT_POSITIVES {
    take:
    ppis_ch          // tuple(meta, ppis)
    sequences_ch     // tuple(meta, fasta)
    partition_ch     // tuple(meta, partition)
    node_mapping_ch  // tuple(meta, node_mapping)

    main:
    joined = ppis_ch.join(sequences_ch).join(partition_ch).join(node_mapping_ch)
    // tuple(meta, ppis, fasta, partition, node_mapping)

    branched = joined.branch { meta, ppis, fasta, partition, node_mapping ->
        ilp:   meta.split_method == "ilp"
        kahip: true
    }

    gurobi_license_ch = params.gurobi_license
        ? channel.value(file(params.gurobi_license, checkIfExists: true))
        : channel.value([])

    ilp_out   = SOLVE_ILP(branched.ilp, gurobi_license_ch)
    kahip_out = SORT_PPIS(branched.kahip)

    train_ppis  = ilp_out.train_ppis.mix(kahip_out.train_ppis)
    val_ppis    = ilp_out.val_ppis.mix(kahip_out.val_ppis)
    test_ppis   = ilp_out.test_ppis.mix(kahip_out.test_ppis)
    train_fasta = ilp_out.train_fasta.mix(kahip_out.train_fasta)
    val_fasta   = ilp_out.val_fasta.mix(kahip_out.val_fasta)
    test_fasta  = ilp_out.test_fasta.mix(kahip_out.test_fasta)
    sorted_mqc  = ilp_out.mqc.mix(kahip_out.mqc)

    // One channel, one process: each (meta, label, fasta1, fasta2) item
    // becomes its own task, so Nextflow runs both CD-HIT-2D comparisons
    // for every dataset in parallel.
    cdhit_inputs_ch = train_fasta.join(val_fasta).map { meta, t, v -> tuple(meta, "train_val", t, v) }
        .mix(train_fasta.join(test_fasta).map { meta, t, te -> tuple(meta, "train_test", t, te) })

    cdhit_out = CDHIT2D(cdhit_inputs_ch)

    cdhit_branched = cdhit_out.sim.branch {
        meta, label, f ->
            train_val:  label == "train_val"
            train_test: label == "train_test"
    }
    sim_tv = cdhit_branched.train_val.map { meta, label, f -> tuple(meta, f) }
    sim_tt = cdhit_branched.train_test.map { meta, label, f -> tuple(meta, f) }

    nr_inputs = train_ppis.join(val_ppis).join(test_ppis)
        .join(train_fasta).join(val_fasta).join(test_fasta)
        .join(sim_tv).join(sim_tt)

    nr = REMOVE_REDUNDANT(nr_inputs)

    emit:
    train_ppis  = nr.train_ppis
    val_ppis    = nr.val_ppis
    test_ppis   = nr.test_ppis
    train_fasta = nr.train_fasta
    val_fasta   = nr.val_fasta
    test_fasta  = nr.test_fasta
    sorted_mqc  = sorted_mqc
    nr_mqc      = nr.mqc
}
