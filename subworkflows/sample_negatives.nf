include { SAMPLE_NEGATIVES_DEGREE; SAMPLE_NEGATIVES_ILP } from '../processes/negative_sampling'

// Output filename suffix: "" for a single negative set, "_<negset>" for several.
def negsuffix(negsets, negset) {
    negsets.size() > 1 ? "_${negset}" : ""
}

// Re-emit the one test_realistic item under every negative set of the row. Every
// set samples it uniformly with the same seed, so it is sampled once, unsuffixed.
def broadcastRealistic(ch, negsets_ch) {
    def branched = ch.branch { meta, negset, label, suffix, f ->
        shared: label == "test_realistic"
        other:  true
    }
    return branched.other.map { meta, negset, label, suffix, f -> tuple(meta, negset, label, f) }
        .mix(branched.shared.combine(negsets_ch, by: 0).flatMap { meta, negset, label, suffix, f, negsets ->
            negsets.collect { ns -> tuple(meta, ns, label, f) }
        })
}

// Samples negative PPIs for each split, using the degree-preserving
// sampler, a fully-uniform variant, or the bias-aware ILP sampler --
// selected per entry of meta.negative_sampling_method.
workflow SAMPLE_NEGATIVES {
    take:
    train_ppis            // tuple(meta, path)
    val_ppis               // tuple(meta, path)
    test_ppis               // tuple(meta, path)
    species_ch              // tuple(meta, path)
    go_annotations_ch        // tuple(meta, path)
    candidate_network_ch      // tuple(meta, candidate_network_or_[])
    negsets_ch                   // tuple(meta, [negset, ...])

    main:
    // One channel, one process: each (meta, negset, label) item becomes its own
    // task. negset is a tuple field, not a meta key, so joins on meta still match.
    fan_splits = train_ppis.map { meta, f -> tuple(meta, "train", f, 1.0) }
        .mix(val_ppis.map  { meta, f -> tuple(meta, "val", f, 1.0) })
        .mix(test_ppis.map { meta, f -> tuple(meta, "test_balanced", f, 1.0) })
        .combine(negsets_ch, by: 0)
        .flatMap { meta, label, f, ratio, negsets ->
            negsets.collect { ns -> tuple(meta, ns, label, negsuffix(negsets, ns), f, ratio) }
        }

    // test_realistic is sampled once per row (see broadcastRealistic()), routed
    // through the first negative set's branch.
    realistic_split = test_ppis.combine(negsets_ch, by: 0)
        .map { meta, f, negsets ->
            if (negsets.size() > 1) {
                log.info "${meta.id}: sampling test_realistic once for negative sets ${negsets.join(', ')}, published unsuffixed"
            }
            tuple(meta, negsets[0], "test_realistic", "", f, 10.0)
        }

    splits_ch = fan_splits.mix(realistic_split)

    branched = splits_ch.branch { meta, negset, label, suffix, f, ratio ->
        // test_realistic always falls through to "default" (uniform), even
        // under "ilp" -- it simulates an uncontrolled random screen, so
        // bias-matching it would defeat the point.
        ilp:     negset in ["ilp", "ilp_candidates"] && label != "test_realistic"
        // Deliberately naive baseline: uniform negatives for every split,
        // paired with split_method=random, to showcase the degree shortcut
        // instead of avoiding it -- see bin/bias_analysis.py's
        // topology_shortcut attribute.
        uniform: negset == "uniform"
        default: true
    }

    neg_gurobi_license_ch = params.gurobi_license
        ? channel.value(file(params.gurobi_license, checkIfExists: true))
        : channel.value([])

    // species/go_annotations/candidate_network are one-per-dataset;
    // combine(by: 0) broadcasts each dataset's single file to every one
    // of that dataset's splits, rather than a full cross-join.
    ilp_inputs = branched.ilp
        .combine(species_ch, by: 0)
        .combine(go_annotations_ch, by: 0)
        .combine(candidate_network_ch, by: 0)
        // Only the "ilp_candidates" negative set is restricted to the candidate network.
        .map { meta, negset, label, suffix, f, ratio, species, go, cand ->
            tuple(meta, negset, label, suffix, f, ratio, species, go, negset == "ilp_candidates" ? cand : [])
        }
    ilp_out = SAMPLE_NEGATIVES_ILP(ilp_inputs, neg_gurobi_license_ch)

    uniform_inputs = branched.uniform.map { meta, negset, label, suffix, f, ratio ->
        tuple(meta, negset, label, suffix, f, ratio, true)
    }
    // label == "test_realistic" here also catches the ILP negative sets'
    // test_realistic item, since the branch above routes it here instead of
    // to SAMPLE_NEGATIVES_ILP.
    default_inputs = branched.default.map { meta, negset, label, suffix, f, ratio ->
        tuple(meta, negset, label, suffix, f, ratio, label == "test_realistic")
    }
    degree_out = SAMPLE_NEGATIVES_DEGREE(uniform_inputs.mix(default_inputs))

    fam_labelled = ilp_out.labelled.mix(degree_out.labelled)  // tuple(meta, negset, label, suffix, csv)
    neg_mqc      = ilp_out.mqc.mix(degree_out.mqc)
    neg_labelled = broadcastRealistic(fam_labelled, negsets_ch)

    neg_branched = neg_labelled.branch {
        meta, negset, label, f ->
            train:          label == "train"
            val:            label == "val"
            test_balanced:  label == "test_balanced"
            test_realistic: label == "test_realistic"
    }

    emit:
    train          = neg_branched.train.map          { meta, negset, label, f -> tuple(meta, negset, f) }
    val            = neg_branched.val.map            { meta, negset, label, f -> tuple(meta, negset, f) }
    test_balanced  = neg_branched.test_balanced.map  { meta, negset, label, f -> tuple(meta, negset, f) }
    test_realistic = neg_branched.test_realistic.map { meta, negset, label, f -> tuple(meta, negset, f) }
    mqc            = neg_mqc
    labelled       = neg_labelled  // tuple(meta, negset, label, csv)
}
