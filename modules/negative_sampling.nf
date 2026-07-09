process SAMPLE_NEGATIVES_DEGREE {
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
        --max-candidates    200000 \\
        --verbose
    """
}
