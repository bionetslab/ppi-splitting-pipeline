process FETCH_DATA {
    publishDir "${params.outdir}/data", mode: 'copy'
    tag "fetch"

    input:
    path ppis

    output:
    path "sequences.fasta",    emit: sequences
    path "go_annotations.tsv", emit: go_annotations

    script:
    """
    fetch_data.py ${ppis} sequences.fasta go_annotations.tsv
    """
}

process GET_LENGTHS {
    tag "lengths"

    input:
    path fasta

    output:
    path "lengths.tsv"

    script:
    """
    { printf 'protein_id\\tlength\\n'; \
      awk '/^>/{if(acc) print acc"\\t"len; acc=substr(\$1,2); len=0; next} {len+=length(\$0)} END{if(acc) print acc"\\t"len}' ${fasta} \
          | sort; \
    } > lengths.tsv
    """
}

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
    tag "kahip"

    input:
    path graph

    output:
    path "partitioned_proteome.txt"

    script:
    """
    kaffpa \\
        ${graph} \\
        --seed=${params.kahip_seed} \\
        --output_filename=partitioned_proteome.txt \\
        --k=${params.kahip_k} \\
        --preconfiguration=${params.kahip_preconfiguration}
    """
}

process SORT_PPIS {
    tag "sort"

    input:
    path ppis
    path partition
    path fasta
    path node_mapping

    output:
    path "train.csv",          emit: train_ppis
    path "val.csv",            emit: val_ppis
    path "test.csv",           emit: test_ppis
    path "train.fasta",        emit: train_fasta
    path "val.fasta",          emit: val_fasta
    path "test.fasta",         emit: test_fasta
    path "*_mqc.tsv",          emit: mqc

    script:
    """
    sort_ppis.py \\
        --ppis         ${ppis} \\
        --partition    ${partition} \\
        --fasta        ${fasta} \\
        --node_mapping ${node_mapping}
    """
}

process CDHIT {
    input:
    path db1_fasta
    path db2_fasta

    output:
    path "cdhit.out"

    script:
    """
    cd-hit-2d \\
        -i  ${db1_fasta} \\
        -i2 ${db2_fasta} \\
        -o  cdhit.out \\
        -c  ${params.cdhit_identity} \\
        -n  ${params.cdhit_wordsize} \\
        -T  ${task.cpus} \\
        -M  4000
    """
}

process REMOVE_REDUNDANT {
    tag "remove_redundant"

    input:
    path train_ppis
    path val_ppis
    path test_ppis
    path train_fasta
    path val_fasta
    path test_fasta
    path sim_train_val,  stageAs: 'sim_train_val.out'
    path sim_train_test, stageAs: 'sim_train_test.out'

    output:
    path "train_nr.csv",              emit: train_ppis
    path "val_nr.csv",                emit: val_ppis
    path "test_nr.csv",               emit: test_ppis
    path "train_nr.fasta",            emit: train_fasta
    path "val_nr.fasta",              emit: val_fasta
    path "test_nr.fasta",             emit: test_fasta
    path "*_mqc.tsv",                 emit: mqc

    script:
    """
    remove_redundant.py \\
        --train_ppis     ${train_ppis} \\
        --val_ppis       ${val_ppis} \\
        --test_ppis      ${test_ppis} \\
        --train_fasta    ${train_fasta} \\
        --val_fasta      ${val_fasta} \\
        --test_fasta     ${test_fasta} \\
        --sim_train_val  ${sim_train_val} \\
        --sim_train_test ${sim_train_test}
    """
}

process SAMPLE_NEGATIVES {
    publishDir "${params.outdir}", mode: 'copy', saveAs: { f -> f.endsWith('_mqc.tsv') ? null : f }
    tag "negatives"

    input:
    path train_ppis
    path val_ppis
    path test_ppis

    output:
    path "train.csv",          emit: train
    path "val.csv",            emit: val
    path "test_balanced.csv",  emit: test_balanced
    path "test_realistic.csv", emit: test_realistic
    path "*_mqc.tsv",          emit: mqc

    script:
    """
    sample_negatives.py \\
        --train ${train_ppis} \\
        --val   ${val_ppis} \\
        --test  ${test_ppis}
    """
}

process EMBED_SEQUENCES {
    publishDir "${params.outdir}/data"
    tag "embed"

    input:
    path train_fasta
    path val_fasta
    path test_fasta

    output:
    path "embeddings.npz"

    script:
    """
    embed_sequences.py \\
        --fasta ${train_fasta} ${val_fasta} ${test_fasta} \\
        --model ${params.embedding_model}
    """
}

process TRAIN_CLASSIFIER {
    tag "classifier"

    input:
    path train_csv
    path val_csv
    path test_balanced_csv
    path test_realistic_csv
    path embeddings

    output:
    path "classifier_metrics_mqc.tsv", emit: mqc

    script:
    """
    train_classifier.py \\
        --train          ${train_csv} \\
        --val            ${val_csv} \\
        --test_balanced  ${test_balanced_csv} \\
        --test_realistic ${test_realistic_csv} \\
        --embeddings     ${embeddings} \\
        --seed           ${params.seed}
    """
}

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
