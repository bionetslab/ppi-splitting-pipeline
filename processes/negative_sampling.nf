process SAMPLE_NEGATIVES_DEGREE {
    publishDir(path: { "${params.outdir}/${meta.id}" }, mode: 'copy', saveAs: { f -> f.endsWith('_mqc.tsv') ? null : f })
    tag "${meta.id}_${label}"
    label 'error_retry'

    input:
    // Staged under a fixed name distinct from "${label}.csv" -- otherwise for
    // label in {train, val} the output would overwrite the staged input
    // symlink, corrupting the upstream task's cached output on -resume.
    tuple val(meta), val(label), path(positives, stageAs: 'positives_in.csv'), val(ratio), val(uniform)  // label: "train" | "val" | "test_balanced" | "test_realistic"

    output:
    tuple val(meta), val(label), path("${label}.csv"), emit: labelled
    tuple val(meta), path("${label}*_mqc.tsv"),         emit: mqc

    script:
    def uniform_flag = uniform ? '--uniform' : ''
    """
    sample_negatives.py \\
        --positives  ${positives} \\
        --output     ${label}.csv \\
        --split-name ${label} \\
        --ratio      ${ratio} \\
        ${uniform_flag} \\
        --seed       ${params.seed} \\
        --id         ${meta.id}
    """
}

process SAMPLE_NEGATIVES_ILP {
    publishDir(path: { "${params.outdir}/${meta.id}" }, mode: 'copy', saveAs: { f -> f.endsWith('_mqc.tsv') ? null : f })
    tag "${meta.id}_${label}"
    label 'error_retry'
    label 'gurobi'

    input:
    // positives is staged under a fixed name distinct from any "${label}.csv"
    // output -- see SAMPLE_NEGATIVES_DEGREE above for why that matters.
    tuple val(meta), val(label), path(positives, stageAs: 'positives_in.csv'), val(neg_ratio), path(species), path(go_annotations), path(candidate_network)  // label: "train" | "val" | "test_balanced" | "test_realistic"; candidate_network optional, [] if unset
    path gurobi_license  // optional; [] if params.gurobi_license is unset

    output:
    tuple val(meta), val(label), path("${label}.csv"), emit: labelled
    tuple val(meta), path("${label}*_mqc.tsv"),         emit: mqc

    script:
    def cand_arg = candidate_network ? "--candidate-network ${candidate_network}" : ''
    def lic_arg  = gurobi_license    ? "--gurobi-license ${gurobi_license}"        : ''
    """
    n_positives=\$(( \$(wc -l < ${positives}) - 1 ))  # -1 for the header row
    max_candidates=\$(( 4 * n_positives * ${task.attempt} ))

    sample_negatives_ilp.py \\
        --positives          ${positives} \\
        --output             ${label}.csv \\
        --split-name         ${label} \\
        --neg-ratio          ${neg_ratio} \\
        --species            ${species} \\
        --go-annotations     ${go_annotations} \\
        ${cand_arg} \\
        ${lic_arg} \\
        --alpha-confidence   ${meta.neg_ilp_alpha_confidence} \\
        --alpha-bias         ${meta.neg_ilp_alpha_bias} \\
        --lambda-degree      ${meta.neg_ilp_lambda_degree} \\
        --lambda-taxon-pair  ${meta.neg_ilp_lambda_taxon_pair} \\
        --lambda-self-loop   ${meta.neg_ilp_lambda_self_loop} \\
        --lambda-jaccard     ${meta.neg_ilp_lambda_jaccard} \\
        --degree-bias-mode   ${params.neg_ilp_degree_bias_mode} \\
        --solver             ${params.neg_ilp_solver} \\
        --time-limit         ${meta.neg_ilp_time_limit} \\
        --mip-gap            ${params.neg_ilp_mip_gap} \\
        --threads            ${task.cpus} \\
        --seed               ${params.seed} \\
        --diagnostics-out    ${label}_mqc.tsv \\
        --residuals-out      ${label}_residuals_mqc.tsv \\
        --max-candidates     \$max_candidates \\
        --verbose \\
        --id                 ${meta.id}
    """
}
