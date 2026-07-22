include { SORT_PPIS; SOLVE_ILP; SPLIT_RANDOM; CDHIT2D; REMOVE_REDUNDANT } from '../processes/splitting'

// Assigns PPIs to train/val/test (KaHIP, ILP, or random shuffle, per
// meta.split_method), then runs CD-HIT-2D redundancy removal for the two
// leakage-aware methods. "random" skips CD-HIT -- see SPLIT_RANDOM.
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
        ilp:    meta.split_method == "ilp"
        random: meta.split_method == "random"
        kahip:  true
    }

    gurobi_license_ch = params.gurobi_license
        ? channel.value(file(params.gurobi_license, checkIfExists: true))
        : channel.value([])

    ilp_out    = SOLVE_ILP(branched.ilp, gurobi_license_ch)
    kahip_out  = SORT_PPIS(branched.kahip)
    random_out = SPLIT_RANDOM(branched.random.map { meta, ppis, fasta, partition, node_mapping -> tuple(meta, ppis, fasta) })

    // ilp/kahip go through CD-HIT redundancy removal below; random doesn't
    // -- it feeds straight into the final emit channels further down.
    homology_train_ppis  = ilp_out.train_ppis.mix(kahip_out.train_ppis)
    homology_val_ppis    = ilp_out.val_ppis.mix(kahip_out.val_ppis)
    homology_test_ppis   = ilp_out.test_ppis.mix(kahip_out.test_ppis)
    homology_train_fasta = ilp_out.train_fasta.mix(kahip_out.train_fasta)
    homology_val_fasta   = ilp_out.val_fasta.mix(kahip_out.val_fasta)
    homology_test_fasta  = ilp_out.test_fasta.mix(kahip_out.test_fasta)
    sorted_mqc = ilp_out.mqc.mix(kahip_out.mqc).mix(random_out.mqc)

    // One channel, one process: each (meta, label, fasta1, fasta2) item
    // becomes its own task, so Nextflow runs both CD-HIT-2D comparisons
    // for every dataset in parallel.
    cdhit_inputs_ch = homology_train_fasta.join(homology_val_fasta).map { meta, t, v -> tuple(meta, "train_val", t, v) }
        .mix(homology_train_fasta.join(homology_test_fasta).map { meta, t, te -> tuple(meta, "train_test", t, te) })

    cdhit_out = CDHIT2D(cdhit_inputs_ch)

    cdhit_branched = cdhit_out.sim.branch {
        meta, label, f ->
            train_val:  label == "train_val"
            train_test: label == "train_test"
    }
    sim_tv = cdhit_branched.train_val.map { meta, label, f -> tuple(meta, f) }
    sim_tt = cdhit_branched.train_test.map { meta, label, f -> tuple(meta, f) }

    // ppis_ch (pre-split) is threaded in so REMOVE_REDUNDANT can compute the
    // KaHIP/ILP discard count itself for the PPI Partitioning chart.
    nr_inputs = ppis_ch.join(homology_train_ppis).join(homology_val_ppis).join(homology_test_ppis)
        .join(homology_train_fasta).join(homology_val_fasta).join(homology_test_fasta)
        .join(sim_tv).join(sim_tt)

    nr = REMOVE_REDUNDANT(nr_inputs)

    emit:
    train_ppis  = nr.train_ppis.mix(random_out.train_ppis)
    val_ppis    = nr.val_ppis.mix(random_out.val_ppis)
    test_ppis   = nr.test_ppis.mix(random_out.test_ppis)
    train_fasta = nr.train_fasta.mix(random_out.train_fasta)
    val_fasta   = nr.val_fasta.mix(random_out.val_fasta)
    test_fasta  = nr.test_fasta.mix(random_out.test_fasta)
    sorted_mqc  = sorted_mqc
    nr_mqc      = nr.mqc
}
