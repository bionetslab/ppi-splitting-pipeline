process BIAS_ANALYSIS {
    tag "${attribute}"

    input:
    each attribute
    path train_csv
    path val_csv
    path test_balanced_csv
    path test_realistic_csv
    path blast_tsv
    path embeddings
    path go_annotations
    path species

    output:
    path "*_bias_mqc.tsv", emit: mqc, optional: true

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
    tag "bias_scatter"

    input:
    path tsvs

    output:
    path "bias_scatter_mqc.html", emit: mqc

    script:
    """
    collect_bias.py ${tsvs}
    """
}

process SIMILARITY_HEATMAP {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    tag "heatmap"

    input:
    path train_fasta
    path val_fasta
    path test_fasta
    path blast_tsv

    output:
    path "similarity_heatmap_mqc.html"

    script:
    """
    plot_similarity_heatmap.py \\
        --train_fasta ${train_fasta} \\
        --val_fasta   ${val_fasta} \\
        --test_fasta  ${test_fasta} \\
        --blast       ${blast_tsv} \\
        --max_per_split ${params.heatmap_max_per_split} \\
        --seed        ${params.seed}
    """
}

process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    tag "multiqc"

    input:
    path mqc_files

    output:
    path "multiqc_report.html"
    path "multiqc_report_data"

    script:
    """
    multiqc . --title "PPI Splitting Pipeline" --filename multiqc_report.html
    """
}
