include { BIAS_ANALYSIS; COLLECT_BIAS; SIMILARITY_HEATMAP; MULTIQC } from '../modules/qc'

// Runs the per-attribute bias analyses, collects them into a scatter plot,
// builds the train/val/test similarity heatmap, and assembles the final
// MultiQC report from every stage's diagnostics.
workflow QC {
    take:
    train_ppis
    val_ppis
    test_balanced_ppis
    test_realistic_ppis
    blast_out
    embeddings
    go_annotations_ch
    species_ch
    train_fasta
    val_fasta
    test_fasta
    sorted_mqc
    nr_mqc
    neg_mqc
    clf_mqc

    main:
    same_species_ch = species_ch
        .splitCsv(header: true, sep: '\t')
        .map    { row -> row.taxon_id }
        .collect()
        .map    { ids -> ids.unique() }
        .filter { ids -> ids.size() > 1 }
        .map    { "same_species" }

    bias = BIAS_ANALYSIS(
        channel.of("sequence_similarity", "embedding_similarity",
                   "functional_relatedness_BP", "functional_relatedness_MF",
                   "functional_relatedness_CC", "self_interactions")
               .mix(same_species_ch).collect(),
        train_ppis,
        val_ppis,
        test_balanced_ppis,
        test_realistic_ppis,
        blast_out,
        embeddings,
        go_annotations_ch,
        species_ch
    )
    scatter = COLLECT_BIAS(bias.mqc.collect())
    heatmap = SIMILARITY_HEATMAP(train_fasta, val_fasta, test_fasta, blast_out)

    mqc_files = sorted_mqc
        .mix(nr_mqc)
        .mix(neg_mqc)
        .mix(clf_mqc)
        .mix(bias.mqc)
        .mix(scatter.mqc)
        .mix(heatmap)
        .collect()

    MULTIQC(mqc_files)
}