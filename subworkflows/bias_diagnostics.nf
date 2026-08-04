include { BIAS_ANALYSIS } from '../processes/qc'

// Runs BIAS_ANALYSIS for every attribute applicable to each dataset. Shared
// by the full pipeline's QC subworkflow (subworkflows/qc.nf) and the
// --bias_only shortcut (subworkflows/bias_only.nf), so the attribute list
// and the same_species rule live in exactly one place.
workflow BIAS_DIAGNOSTICS {
    take:
    train_ppis            // tuple(meta, path)
    val_ppis
    test_balanced_ppis
    test_realistic_ppis
    blast_out
    embeddings
    go_annotations_ch
    species_ch

    main:
    // Whether to include "same_species" depends on each dataset's own
    // species.tsv, so it's computed per-dataset here rather than with a
    // single run-wide collect().
    attrs_ch = species_ch.map { meta, sp ->
        def taxa = sp.splitCsv(header: true, sep: '\t').collect { it.taxon_id }.unique()
        def attrs = ["sequence_similarity", "embedding_similarity",
                     "functional_relatedness_BP", "functional_relatedness_MF",
                     "functional_relatedness_CC", "self_interactions",
                     "topology_shortcut"]
        if (taxa.size() > 1) attrs << "same_species"
        tuple(meta, attrs)
    }.flatMap { meta, attrs -> attrs.collect { a -> tuple(meta, a) } }

    // train/val/test/blast/embeddings/go/species are one-per-dataset;
    // combine(by: 0) broadcasts each dataset's single set of files to
    // every one of that dataset's attributes, rather than a full cross-join.
    bias_inputs = attrs_ch
        .combine(train_ppis,          by: 0)
        .combine(val_ppis,            by: 0)
        .combine(test_balanced_ppis,  by: 0)
        .combine(test_realistic_ppis, by: 0)
        .combine(blast_out,           by: 0)
        .combine(embeddings,          by: 0)
        .combine(go_annotations_ch,   by: 0)
        .combine(species_ch,          by: 0)

    bias = BIAS_ANALYSIS(bias_inputs)

    emit:
    mqc = bias.mqc
}