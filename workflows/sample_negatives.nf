// Aliased because the "random" sampler process is itself named SAMPLE_NEGATIVES,
// which would otherwise collide with this subworkflow's own name.
include { SAMPLE_NEGATIVES_DEGREE; SAMPLE_NEGATIVES_ILP } from '../modules/negative_sampling'

// Samples negative PPIs for each split, using either the default random
// sampler or the bias-aware ILP sampler, selected per-dataset via
// meta.negative_sampling_method.
workflow SAMPLE_NEGATIVES {
    take:
    train_ppis            // tuple(meta, path)
    val_ppis               // tuple(meta, path)
    test_ppis               // tuple(meta, path)
    species_ch              // tuple(meta, path)
    go_annotations_ch        // tuple(meta, path)
    candidate_network_ch      // tuple(meta, candidate_network_or_[])

    main:
    // One channel, one process: each (meta, label) item becomes its own
    // task, so Nextflow runs all four splits for every dataset in parallel.
    splits_ch = train_ppis.map { meta, f -> tuple(meta, "train", f, 1.0) }
        .mix(val_ppis.map  { meta, f -> tuple(meta, "val", f, 1.0) })
        .mix(test_ppis.map { meta, f -> tuple(meta, "test_balanced", f, 1.0) })
        .mix(test_ppis.map { meta, f -> tuple(meta, "test_realistic", f, 10.0) })

    branched = splits_ch.branch { meta, label, f, ratio ->
        ilp:     meta.negative_sampling_method == "ilp"
        default: true
    }

    neg_gurobi_license_ch = params.gurobi_license
        ? channel.value(file(params.gurobi_license, checkIfExists: true))
        : channel.value([])

    // species/go_annotations/candidate_network are one-per-dataset;
    // combine(by: 0) broadcasts each dataset's single file to every one
    // of that dataset's (up to 4) splits, rather than a full cross-join.
    ilp_inputs = branched.ilp
        .combine(species_ch, by: 0)
        .combine(go_annotations_ch, by: 0)
        .combine(candidate_network_ch, by: 0)
    ilp_out = SAMPLE_NEGATIVES_ILP(ilp_inputs, neg_gurobi_license_ch)

    default_inputs = branched.default.map { meta, label, f, ratio ->
        tuple(meta, label, f, ratio, label == "test_realistic")
    }
    default_out = SAMPLE_NEGATIVES_DEGREE(default_inputs)

    neg_labelled = ilp_out.labelled.mix(default_out.labelled)
    neg_mqc      = ilp_out.mqc.mix(default_out.mqc)

    neg_branched = neg_labelled.branch {
        meta, label, f ->
            train:          label == "train"
            val:            label == "val"
            test_balanced:  label == "test_balanced"
            test_realistic: label == "test_realistic"
    }

    emit:
    train          = neg_branched.train.map          { meta, label, f -> tuple(meta, f) }
    val            = neg_branched.val.map            { meta, label, f -> tuple(meta, f) }
    test_balanced  = neg_branched.test_balanced.map  { meta, label, f -> tuple(meta, f) }
    test_realistic = neg_branched.test_realistic.map { meta, label, f -> tuple(meta, f) }
    mqc            = neg_mqc
}
