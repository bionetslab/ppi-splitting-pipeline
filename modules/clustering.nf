process RUN_BLAST {
    publishDir "${params.outdir}/similarities", mode: 'copy'
    tag "blast"

    input:
    path fasta

    output:
    path "all_vs_all.tsv"

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
    publishDir "${params.outdir}/similarities", mode: 'copy'
    tag "metis"

    input:
    path blast_results
    path lengths

    output:
    path "similarity.graph", emit: graph
    path "node_mapping.tsv", emit: node_mapping

    script:
    """
    make_metis.py \\
        ${blast_results} \\
        ${lengths} \\
        similarity.graph \\
        node_mapping.tsv \\
        --edge_weight ${params.edge_weight}
    """
}

process RUN_KAHIP {
    publishDir "${params.outdir}/similarities", mode: 'copy'
    tag "kahip: k=${k}"

    input:
    path graph
    val k

    output:
    path "partitioned_proteome.txt"

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
