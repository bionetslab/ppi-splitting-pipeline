process RUN_BLAST {
    publishDir(path: { "${params.outdir}/${meta.id}/similarities" }, mode: 'copy')
    tag "${meta.id}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("all_vs_all.tsv")

    script:
    """
    makeblastdb -dbtype prot -in ${fasta} -out blastdb -parse_seqids
    blastp \\
        -query ${fasta} \\
        -db blastdb \\
        -outfmt "6 qseqid sseqid evalue bitscore pident" \\
        -max_hsps 1 \\
        -num_threads ${task.cpus} \\
        -out all_vs_all.tsv
    """
}

process MAKE_METIS {
    publishDir(path: { "${params.outdir}/${meta.id}/similarities" }, mode: 'copy')
    tag "${meta.id}"

    input:
    tuple val(meta), path(blast_results), path(lengths)

    output:
    tuple val(meta), path("similarity.graph"), emit: graph
    tuple val(meta), path("node_mapping.tsv"), emit: node_mapping

    script:
    """
    make_metis.py \\
        ${blast_results} \\
        ${lengths} \\
        similarity.graph \\
        node_mapping.tsv \\
        --edge_weight ${meta.edge_weight}
    """
}

process RUN_KAHIP {
    publishDir(path: { "${params.outdir}/${meta.id}/similarities" }, mode: 'copy')
    tag "${meta.id}: k=${k}"
    label 'error_retry'

    input:
    tuple val(meta), path(graph), val(k)

    output:
    tuple val(meta), path("partitioned_proteome.txt")

    script:
    """
    kaffpa \\
        ${graph} \\
        --seed=${params.kahip_seed} \\
        --output_filename=partitioned_proteome.txt \\
        --k=${k} \\
        --preconfiguration=${params.kahip_preconfiguration}
    """
}
