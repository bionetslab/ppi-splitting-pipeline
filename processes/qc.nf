// One task per (dataset, attribute, negative set).
process BIAS_ANALYSIS {
    // Per negative set: the filename carries no set name.
    publishDir(path: { "${params.outdir}/${meta.id}/bias/${negset}" }, mode: 'copy')
    tag "${meta.id}_${negset}_${attribute}"
    label 'error_retry'

    input:
    tuple val(meta), val(attribute), val(negset),
          path(train_csv), path(val_csv), path(test_balanced_csv), path(test_realistic_csv),
          path(blast_tsv), path(embeddings), path(go_annotations), path(species)

    output:
    tuple val(meta), val(negset), path("*_bias_mqc.tsv"), emit: mqc, optional: true

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

// id is the MultiQC display label for one (dataset, negative set), not meta.id.
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

// `id` (meta.id) is the publishDir component, `label` the MultiQC display label;
// they differ for a row with several negative sets.
process SIMILARITY_HEATMAP {
    publishDir(path: { "${params.outdir}/${id}/multiqc" }, mode: 'copy')
    tag "${label}_heatmap"

    input:
    tuple val(id), val(label), path(train_fasta), path(val_fasta), path(test_fasta), path(blast_tsv)

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
        --id          ${label}
    """
}

// One combined report for the whole run. stageAs: "?/*" stages each
// dataset's contributions into its own numbered subdirectory, since
// same-named files (e.g. classifier_metrics_mqc.tsv) would otherwise
// collide -- MultiQC scans subdirectories recursively, so this is transparent.
// relabel_mqc.py applies the run-wide dataset order here, so it stays out of the
// emitting tasks' hashes, and pins the section order MultiQC would otherwise take
// from channel arrival.
process MULTIQC {
    publishDir(path: { "${params.outdir}/multiqc" }, mode: 'copy')
    tag "multiqc"

    input:
    path multiqc_files, stageAs: "?/*"
    val  mqc_order  // comma-joined display labels, in report order

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_report_data", emit: data

    script:
    """
    relabel_mqc.py . \\
        --out       relabelled \\
        --mqc-order '${mqc_order}' \\
        --config    multiqc_config.yml

    multiqc relabelled \\
        -c multiqc_config.yml \\
        --title "PPI Splitting Pipeline" \\
        --filename multiqc_report.html
    """
}
