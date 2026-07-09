// Aliased because the "random" sampler process is itself named SAMPLE_NEGATIVES,
// which would otherwise collide with this subworkflow's own name.
include { SAMPLE_NEGATIVES_DEGREE; SAMPLE_NEGATIVES_ILP } from '../modules/negative_sampling'

// Samples negative PPIs for each split, using either the default random
// sampler or the bias-aware ILP sampler, selected via
// params.negative_sampling_method.
workflow SAMPLE_NEGATIVES {
    take:
    train_ppis
    val_ppis
    test_ppis
    species_ch
    go_annotations_ch

    main:
    if (params.negative_sampling_method == "ilp") {
        candidate_network_ch = params.candidate_network
            ? channel.value(file(params.candidate_network, checkIfExists: true))
            : channel.value([])
        neg_gurobi_license_ch = params.gurobi_license
            ? channel.value(file(params.gurobi_license, checkIfExists: true))
            : channel.value([])

        // One channel, one process: each item (label, positives, neg_ratio)
        // becomes its own task, so Nextflow runs all four splits in parallel.
        neg_splits_ch = train_ppis.map { f -> tuple("train", f, 1.0) }
            .mix(val_ppis.map  { f -> tuple("val", f, 1.0) })
            .mix(test_ppis.map { f -> tuple("test_balanced", f, 1.0) })
            .mix(test_ppis.map { f -> tuple("test_realistic", f, 10.0) })

        neg_out = SAMPLE_NEGATIVES_ILP(
            neg_splits_ch, species_ch, go_annotations_ch, candidate_network_ch, neg_gurobi_license_ch
        )
    } else {
        // One channel, one process: each item (label, positives, ratio, uniform)
        // becomes its own task, so Nextflow runs all four splits in parallel.
        neg_splits_ch = train_ppis.map { f -> tuple("train", f, 1.0, false) }
            .mix(val_ppis.map  { f -> tuple("val", f, 1.0, false) })
            .mix(test_ppis.map { f -> tuple("test_balanced", f, 1.0, false) })
            .mix(test_ppis.map { f -> tuple("test_realistic", f, 10.0, true) })

        neg_out = SAMPLE_NEGATIVES_DEGREE(neg_splits_ch)
    }

    neg_branched = neg_out.labelled.branch {
        label, f ->
            train:          label == "train"
            val:            label == "val"
            test_balanced:  label == "test_balanced"
            test_realistic: label == "test_realistic"
    }

    emit:
    train          = neg_branched.train.map { label, f -> f }
    val            = neg_branched.val.map { label, f -> f }
    test_balanced  = neg_branched.test_balanced.map { label, f -> f }
    test_realistic = neg_branched.test_realistic.map { label, f -> f }
    mqc            = neg_out.mqc
}