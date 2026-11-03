include { BIAS_ANALYSIS; COLLECT_BIAS } from '../processes/qc'

// Runs BIAS_ANALYSIS for every attribute applicable to each dataset and
// negative set, and collects the results into one scatter per pair. Shared by
// the full pipeline and --bias_only, so the attribute list and the
// same_species rule live in exactly one place.
workflow BIAS_DIAGNOSTICS {
    take:
    train_ppis            // tuple(meta, negset, path)
    val_ppis
    test_balanced_ppis
    test_realistic_ppis
    blast_out             // tuple(meta, path)
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

    // join(by: [0, 1]): a join on meta alone would pair one negative set's
    // train CSV with another's val CSV.
    negset_splits = train_ppis.join(val_ppis, by: [0, 1])
        .join(test_balanced_ppis,  by: [0, 1])
        .join(test_realistic_ppis, by: [0, 1])
    // tuple(meta, negset, train, val, test_balanced, test_realistic)

    // blast/embeddings/go/species are one-per-dataset; combine(by: 0) broadcasts
    // them to every (attribute, negative set) pair rather than a full cross-join.
    bias_inputs = negset_splits
        .combine(attrs_ch, by: 0)
        .map { meta, negset, train, val, tb, tr, attr -> tuple(meta, attr, negset, train, val, tb, tr) }
        .combine(blast_out,         by: 0)
        .combine(embeddings,        by: 0)
        .combine(go_annotations_ch, by: 0)
        .combine(species_ch,        by: 0)

    bias = BIAS_ANALYSIS(bias_inputs)

    // Group by the negset-qualified label so a row's sets get separate scatters.
    bias_by_negset = bias.mqc.flatMap { meta, negset, f ->
        def files = (f instanceof List) ? f : [f]
        def multi = meta.negative_sampling_method.toString().tokenize(',').size() > 1
        files.collect { ff -> tuple((multi ? "${meta.id}_${negset}" : "${meta.id}").toString(), ff) }
    }
    scatter = COLLECT_BIAS(bias_by_negset.groupTuple())

    emit:
    scatter = scatter.mqc   // tuple(label, bias_scatter_mqc.html)
}
