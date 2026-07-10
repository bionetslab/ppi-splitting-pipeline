process BIAS_ANALYSIS {
    tag "${meta.id}_${attribute}"

    input:
    tuple val(meta), val(attribute),
          path(train_csv), path(val_csv), path(test_balanced_csv), path(test_realistic_csv),
          path(blast_tsv), path(embeddings), path(go_annotations), path(species)

    output:
    tuple val(meta), path("*_bias_mqc.tsv"), emit: mqc, optional: true

    script:
    """
    bias_analysis.py \\
        --attribute       ${attribute} \\
        --train           ${train_csv} \\
        --val             ${val_csv} \\
        --test_balanced   ${test_balanced_csv} \\
        --test_realistic  ${test_realistic_csv} \\
        --blast           ${blast_tsv} \\
        --embeddings      ${embeddings} \\
        --go_annotations  ${go_annotations} \\
        --species         ${species} \\
        --seed            ${params.seed}
    """
}

process COLLECT_BIAS {
    tag "${id}_bias_scatter"

    input:
    tuple val(id), path(tsvs)

    output:
    tuple val(id), path("bias_scatter_mqc.html"), emit: mqc

    script:
    """
    collect_bias.py ${tsvs} --id ${id}
    """
}

process SIMILARITY_HEATMAP {
    publishDir(path: { "${params.outdir}/${id}/multiqc" }, mode: 'copy')
    tag "${id}_heatmap"

    input:
    tuple val(id), path(train_fasta), path(val_fasta), path(test_fasta), path(blast_tsv)

    output:
    tuple val(id), path("similarity_heatmap_mqc.html")

    script:
    """
    plot_similarity_heatmap.py \\
        --train_fasta ${train_fasta} \\
        --val_fasta   ${val_fasta} \\
        --test_fasta  ${test_fasta} \\
        --blast       ${blast_tsv} \\
        --max_per_split ${params.heatmap_max_per_split} \\
        --seed        ${params.seed} \\
        --id          ${id}
    """
}

// One combined report for the whole run. stageAs: "?/*" stages every
// dataset's contributions into its own numbered subdirectory, since
// different datasets independently write same-named files (e.g. every
// dataset has its own classifier_metrics_mqc.tsv) that would otherwise
// collide in one task's work directory -- MultiQC scans subdirectories
// recursively, so this is transparent to it.
process MULTIQC {
    publishDir(path: { "${params.outdir}/multiqc" }, mode: 'copy')
    tag "multiqc"

    input:
    path multiqc_files, stageAs: "?/*"

    output:
    path "multiqc_report.html"
    path "multiqc_report_data"

    script:
    """
    multiqc . --title "PPI Splitting Pipeline" --filename multiqc_report.html
    """
}
