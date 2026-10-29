include { SIMILARITY_HEATMAP; MULTIQC } from '../processes/qc'

// Some mqc-emitting processes glob-match more than one file per task, which
// Nextflow packs into a List. Flatten to one (id, file) pair per file.
def flattenMqc(ch) {
    ch.flatMap { meta, f ->
        def files = (f instanceof List) ? f : [f]
        files.collect { ff -> tuple(meta.id, ff) }
    }
}

// Builds the train/val/test similarity heatmap and assembles one combined
// MultiQC report for the whole run from every dataset's diagnostics.
workflow QC {
    take:
    bias_scatter          // tuple(label, bias_scatter_mqc.html)
    blast_out
    train_fasta
    val_fasta
    test_fasta
    sorted_mqc
    nr_mqc
    neg_mqc
    clf_mqc

    main:
    heatmap_inputs = train_fasta.join(val_fasta).join(test_fasta).join(blast_out)
        .map { meta, t, v, te, b -> tuple(meta.id, t, v, te, b) }
    heatmap = SIMILARITY_HEATMAP(heatmap_inputs)

    // Bias tables are published by BIAS_ANALYSIS; only the scatter goes into the report.
    mqc_files = flattenMqc(sorted_mqc)
        .mix(flattenMqc(nr_mqc))
        .mix(flattenMqc(neg_mqc))
        .mix(flattenMqc(clf_mqc))
        .mix(bias_scatter)
        .mix(heatmap)
        .map { id, f -> f }
        .collect()

    multiqc = MULTIQC(mqc_files)

    emit:
    multiqc_report = multiqc.report
}
