include { COLLECT_BIAS; SIMILARITY_HEATMAP; MULTIQC } from '../processes/qc'

// Some mqc-emitting processes glob-match more than one file per task, which
// Nextflow packs into a List. Flatten to one (id, file) pair per file so
// groupTuple() below doesn't nest a List inside the grouped list.
def flattenMqc(ch) {
    ch.flatMap { meta, f ->
        def files = (f instanceof List) ? f : [f]
        files.collect { ff -> tuple(meta.id, ff) }
    }
}

// Collects the bias diagnostics (computed upstream by BIAS_DIAGNOSTICS,
// subworkflows/bias_diagnostics.nf) into a scatter plot, builds the
// train/val/test similarity heatmap, and assembles one combined MultiQC
// report for the whole run from every dataset's diagnostics.
workflow QC {
    take:
    bias_mqc               // tuple(meta, path) -- from BIAS_DIAGNOSTICS
    blast_out
    train_fasta
    val_fasta
    test_fasta
    sorted_mqc
    nr_mqc
    neg_mqc
    clf_mqc

    main:
    scatter = COLLECT_BIAS(flattenMqc(bias_mqc).groupTuple())

    heatmap_inputs = train_fasta.join(val_fasta).join(test_fasta).join(blast_out)
        .map { meta, t, v, te, b -> tuple(meta.id, t, v, te, b) }
    heatmap = SIMILARITY_HEATMAP(heatmap_inputs)

    mqc_files = flattenMqc(sorted_mqc)
        .mix(flattenMqc(nr_mqc))
        .mix(flattenMqc(neg_mqc))
        .mix(flattenMqc(clf_mqc))
        .mix(scatter.mqc)
        .mix(heatmap)
        .map { id, f -> f }
        .collect()

    MULTIQC(mqc_files)
}
