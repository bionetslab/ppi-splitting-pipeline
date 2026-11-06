process BIAS_ANALYSIS {
    publishDir(path: { "${params.outdir}/${meta.id}/multiqc" }, mode: 'copy')
    tag "${meta.id}_${attribute}"
    label 'error_retry'

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
        --seed            ${params.bias_seed}
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
        --seed        ${params.bias_seed} \\
        --id          ${id}
    """
}

// --robustness_test: aggregates one dataset's classifier metrics over every
// (split, sampling, training) seed combination. seeds[i] holds the
// [split, sampling, training] seeds of metrics[i]; both lists come from the
// same groupTuple(), so their order matches.
process ROBUSTNESS_SUMMARY {
    publishDir(path: { "${params.outdir}/${id}/robustness" }, mode: 'copy')
    tag "${id}_robustness"

    input:
    tuple val(id), val(seeds), path(metrics, stageAs: "?/*")

    output:
    tuple val(id), path("robustness_summary_mqc.tsv"), emit: mqc
    path "robustness_runs.tsv"

    script:
    def files = metrics instanceof List ? metrics : [metrics]
    def rows = [files, seeds].transpose().collect { f, s -> ([f] + s).join(' ') }.join(' ')
    """
    printf '%s\\t%s\\t%s\\t%s\\n' file split_seed sampling_seed training_seed ${rows} > manifest.tsv
    robustness_summary.py --manifest manifest.tsv --id ${id}
    """
}

// One combined report for the whole run. stageAs: "?/*" stages each
// dataset's contributions into its own numbered subdirectory, since
// same-named files (e.g. classifier_metrics_mqc.tsv) would otherwise
// collide -- MultiQC scans subdirectories recursively, so this is transparent.
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
