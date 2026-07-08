process FETCH_DATA {
    publishDir "${params.outdir}/data", mode: 'copy'
    tag "fetch"

    input:
    path ppis

    output:
    path "sequences.fasta",    emit: sequences
    path "go_annotations.tsv", emit: go_annotations
    path "species.tsv",        emit: species

    script:
    """
    fetch_data.py ${ppis} sequences.fasta go_annotations.tsv species.tsv
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

process CDHIT2D {
    tag "cdhit2d_${label}"

    input:
    tuple val(label), path(db1_fasta), path(db2_fasta)  // label: "train_val" | "train_test"

    output:
    tuple val(label), path("cdhit.out"), emit: sim

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

process SOLVE_ILP {
    tag "solve_ilp"

    input:
    path ppis
    path fasta
    path partition
    path node_mapping
    path gurobi_license

    output:
    path "train.csv",   emit: train_ppis
    path "val.csv",     emit: val_ppis
    path "test.csv",    emit: test_ppis
    path "train.fasta", emit: train_fasta
    path "val.fasta",   emit: val_fasta
    path "test.fasta",  emit: test_fasta
    path "*_mqc.tsv",   emit: mqc

    script:
    def license_export = gurobi_license ? "export GRB_LICENSE_FILE=\$PWD/${gurobi_license}" : ""
    """
    ${license_export}
    solve_ilp.py \\
        --ppis          ${ppis} \\
        --fasta         ${fasta} \\
        --partition     ${partition} \\
        --node_mapping  ${node_mapping} \\
        --train-split ${params.train_split} \\
        --val-split   ${params.val_split} \\
        --test-split  ${params.test_split} \\
        --epsilon     ${params.ilp_epsilon} \\
        --max-sec     ${params.ilp_max_sec} \\
        ${params.ilp_solver ? "--solver ${params.ilp_solver}" : ""}
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
    tag "negatives_${label}"

    input:
    tuple val(label), path(positives), val(ratio), val(uniform)  // label: "train" | "val" | "test_balanced" | "test_realistic"

    output:
    tuple val(label), path("${label}.csv"), emit: labelled
    path "${label}*_mqc.tsv",               emit: mqc

    script:
    def uniform_flag = uniform ? '--uniform' : ''
    """
    sample_negatives.py \\
        --positives  ${positives} \\
        --output     ${label}.csv \\
        --split-name ${label} \\
        --ratio      ${ratio} \\
        ${uniform_flag} \\
        --seed       ${params.seed}
    """
}

process SAMPLE_NEGATIVES_ILP {
    publishDir "${params.outdir}", mode: 'copy', saveAs: { f -> f.endsWith('_mqc.tsv') ? null : f }
    tag "negatives_ilp_${label}"

    input:
    tuple val(label), path(positives), val(neg_ratio)  // label: "train" | "val" | "test_balanced" | "test_realistic"
    path species
    path go_annotations
    path candidate_network  // optional; [] if params.candidate_network is unset
    path gurobi_license     // optional; [] if params.gurobi_license is unset

    output:
    tuple val(label), path("${label}.csv"), emit: labelled
    path "${label}*_mqc.tsv",               emit: mqc

    script:
    def cand_arg = candidate_network ? "--candidate-network ${candidate_network}" : ''
    def lic_arg  = gurobi_license    ? "--gurobi-license ${gurobi_license}"        : ''
    """
    sample_negatives_ilp.py \\
        --positives          ${positives} \\
        --output             ${label}.csv \\
        --split-name         ${label} \\
        --neg-ratio          ${neg_ratio} \\
        --species            ${species} \\
        --go-annotations     ${go_annotations} \\
        ${cand_arg} \\
        ${lic_arg} \\
        --alpha-confidence   ${params.neg_ilp_alpha_confidence} \\
        --alpha-bias         ${params.neg_ilp_alpha_bias} \\
        --lambda-degree      ${params.neg_ilp_lambda_degree} \\
        --lambda-taxon-pair  ${params.neg_ilp_lambda_taxon_pair} \\
        --lambda-self-loop   ${params.neg_ilp_lambda_self_loop} \\
        --lambda-jaccard     ${params.neg_ilp_lambda_jaccard} \\
        --degree-bias-mode   ${params.neg_ilp_degree_bias_mode} \\
        --solver             ${params.neg_ilp_solver} \\
        --time-limit         ${params.neg_ilp_time_limit} \\
        --mip-gap            ${params.neg_ilp_mip_gap} \\
        --threads            ${task.cpus} \\
        --seed               ${params.seed} \\
        --diagnostics-out    ${label}_mqc.tsv \\
        --residuals-out      ${label}_residuals_mqc.tsv \\
        --max-candidates    70000000 \\
        --verbose
    """
}

process EMBED_SEQUENCES {
    publishDir "${params.outdir}/data"
    tag "embed"
    label "process_gpu"

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
