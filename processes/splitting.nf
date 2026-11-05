include { mqcRowLabel } from '../helpers/mqc_labels'

// --id takes the row's MultiQC display label; publishDir paths still use meta.id.
process SORT_PPIS {
    tag "${meta.id}"

    input:
    // instances ([] in PPI mode) must stay last: the channel appends it with .join().
    tuple val(meta), path(ppis), path(fasta), path(partition), path(node_mapping), path(instances)

    output:
    tuple val(meta), path("train.csv"),   emit: train_ppis
    tuple val(meta), path("val.csv"),     emit: val_ppis
    tuple val(meta), path("test.csv"),    emit: test_ppis
    tuple val(meta), path("train.fasta"), emit: train_fasta
    tuple val(meta), path("val.fasta"),   emit: val_fasta
    tuple val(meta), path("test.fasta"),  emit: test_fasta
    tuple val(meta), path("*_mqc.tsv"),   emit: mqc, optional: true

    script:
    def inst_arg = instances ? "--instances ${instances}" : ''
    """
    sort_ppis.py \\
        --ppis         ${ppis} \\
        --partition    ${partition} \\
        --fasta        ${fasta} \\
        --node_mapping ${node_mapping} \\
        ${inst_arg}
    """
}

// Deliberately naive baseline: shuffles PPIs randomly instead of using a
// KaHIP partition, so the same protein can land in more than one split --
// see bin/bias_analysis.py's "topology_shortcut" attribute. No CD-HIT runs
// downstream, since it would just strip the shared proteins back out.
process SPLIT_RANDOM {
    tag "${meta.id}"

    input:
    tuple val(meta), path(ppis), path(fasta), path(instances)  // instances last; [] in PPI mode

    output:
    tuple val(meta), path("train.csv"),   emit: train_ppis
    tuple val(meta), path("val.csv"),     emit: val_ppis
    tuple val(meta), path("test.csv"),    emit: test_ppis
    tuple val(meta), path("train.fasta"), emit: train_fasta
    tuple val(meta), path("val.fasta"),   emit: val_fasta
    tuple val(meta), path("test.fasta"),  emit: test_fasta
    tuple val(meta), path("*_mqc.tsv"),   emit: mqc

    script:
    def inst_arg = instances ? "--instances ${instances}" : ''
    """
    sort_ppis_random.py \\
        --ppis        ${ppis} \\
        --fasta       ${fasta} \\
        --train-split ${meta.train_split} \\
        --val-split   ${meta.val_split} \\
        --test-split  ${meta.test_split} \\
        --seed        ${params.seed} \\
        --id          ${mqcRowLabel(meta)} \\
        ${inst_arg}
    """
}

process CDHIT2D {
    tag "${meta.id}_${label}"

    input:
    tuple val(meta), val(label), path(db1_fasta), path(db2_fasta)  // label: "train_val" | "train_test"

    output:
    tuple val(meta), val(label), path("cdhit.out"), emit: sim

    script:
    // cd-hit-2d fails on an empty input (a split fraction of 0); an empty
    // cdhit.out means "no redundancy" and keeps the dataset in the join.
    """
    if [ -s ${db1_fasta} ] && [ -s ${db2_fasta} ]; then
        cd-hit-2d \\
            -i  ${db1_fasta} \\
            -i2 ${db2_fasta} \\
            -o  cdhit.out \\
            -c  ${meta.cdhit_identity} \\
            -n  ${meta.cdhit_wordsize} \\
            -T  ${task.cpus} \\
            -M  4000
    else
        touch cdhit.out
    fi
    """
}

process SOLVE_ILP {
    tag "${meta.id}"
    label 'error_retry'
    label 'gurobi'

    input:
    tuple val(meta), path(ppis), path(fasta), path(partition), path(node_mapping), path(instances)  // instances last
    path gurobi_license

    output:
    tuple val(meta), path("train.csv"),   emit: train_ppis
    tuple val(meta), path("val.csv"),     emit: val_ppis
    tuple val(meta), path("test.csv"),    emit: test_ppis
    tuple val(meta), path("train.fasta"), emit: train_fasta
    tuple val(meta), path("val.fasta"),   emit: val_fasta
    tuple val(meta), path("test.fasta"),  emit: test_fasta
    tuple val(meta), path("*_mqc.tsv"),   emit: mqc, optional: true

    script:
    def license_export = gurobi_license ? "export GRB_LICENSE_FILE=\$PWD/${gurobi_license}" : ""
    def inst_arg       = instances ? "--instances ${instances}" : ''
    """
    ${license_export}
    solve_ilp.py \\
        --ppis          ${ppis} \\
        --fasta         ${fasta} \\
        --partition     ${partition} \\
        --node_mapping  ${node_mapping} \\
        ${inst_arg} \\
        --train-split ${meta.train_split} \\
        --val-split   ${meta.val_split} \\
        --test-split  ${meta.test_split} \\
        --epsilon     ${meta.ilp_epsilon} \\
        --max-sec     ${meta.ilp_max_sec} \\
        --seed        ${params.seed} \\
        ${params.ilp_solver ? "--solver ${params.ilp_solver}" : ""}
    """
}

// DDI mode only. Turns each family pair into up to N domain-instance pairs and
// claims every parent protein for at most one split. Re-emits the family CSVs
// without DDIs left with zero examples.
process SELECT_EXAMPLES {
    publishDir(path: { "${params.outdir}/${meta.id}/examples" }, mode: 'copy', saveAs: { f -> f.endsWith('_mqc.tsv') ? null : f })
    tag "${meta.id}"
    label 'error_retry'
    label 'gurobi'

    input:
    // instances and candidate_network are .join()ed on last, in that order.
    tuple val(meta),
          path(train_ppis), path(val_ppis), path(test_ppis),
          path(train_fasta), path(val_fasta), path(test_fasta),
          path(instances), path(candidate_network)
    path gurobi_license

    output:
    tuple val(meta), path("train_sel.csv"), emit: train_ppis
    tuple val(meta), path("val_sel.csv"),   emit: val_ppis
    tuple val(meta), path("test_sel.csv"),  emit: test_ppis
    tuple val(meta), path("train_examples.csv"), emit: train_examples
    tuple val(meta), path("val_examples.csv"),   emit: val_examples
    tuple val(meta), path("test_examples.csv"),  emit: test_examples
    tuple val(meta), path("train_universe.txt"), emit: train_universe
    tuple val(meta), path("val_universe.txt"),   emit: val_universe
    tuple val(meta), path("test_universe.txt"),  emit: test_universe
    tuple val(meta), path("train_candidate_examples.csv"), emit: train_candidates
    tuple val(meta), path("val_candidate_examples.csv"),   emit: val_candidates
    tuple val(meta), path("test_candidate_examples.csv"),  emit: test_candidates
    // Parents no candidate reached.
    tuple val(meta), path("unclaimed.txt"), emit: unclaimed
    tuple val(meta), path("*_mqc.tsv"),     emit: mqc

    script:
    def license_export = gurobi_license ? "export GRB_LICENSE_FILE=\$PWD/${gurobi_license}" : ""
    def cand_arg       = candidate_network ? "--candidate-network ${candidate_network}" : ''
    // The random baseline must keep its leak, so parents may be shared there.
    def shared_arg     = meta.split_method == "random" ? "--allow-shared-parents" : ''
    def verbose_arg    = params.ddi_select_verbose ? '--verbose' : ''
    """
    ${license_export}
    select_examples.py \\
        --train_ppis  ${train_ppis} \\
        --val_ppis    ${val_ppis} \\
        --test_ppis   ${test_ppis} \\
        --train_fasta ${train_fasta} \\
        --val_fasta   ${val_fasta} \\
        --test_fasta  ${test_fasta} \\
        --instances   ${instances} \\
        ${cand_arg} \\
        ${shared_arg} \\
        --examples-target   ${params.ddi_examples_target} \\
        --shortlist-factor  ${params.ddi_shortlist_factor} \\
        --candidate-factor  ${params.ddi_candidate_factor} \\
        --lambda-diversity  ${params.ddi_lambda_diversity} \\
        --max-sec              ${params.ddi_select_max_sec} \\
        --max-ilp-candidates   ${params.ddi_max_ilp_candidates} \\
        --seed              ${params.seed} \\
        --id                ${mqcRowLabel(meta)} \\
        ${verbose_arg} \\
        ${params.ilp_solver ? "--solver ${params.ilp_solver}" : ""}
    """
}

process REMOVE_REDUNDANT {
    tag "${meta.id}"

    input:
    tuple val(meta),
          path(orig_ppis),
          path(train_ppis), path(val_ppis), path(test_ppis),
          path(train_fasta), path(val_fasta), path(test_fasta),
          path(sim_train_val,  stageAs: 'sim_train_val.out'),
          path(sim_train_test, stageAs: 'sim_train_test.out'),
          path(instances)  // [] in PPI mode; last, matching the .join() order

    output:
    tuple val(meta), path("train_nr.csv"),   emit: train_ppis
    tuple val(meta), path("val_nr.csv"),     emit: val_ppis
    tuple val(meta), path("test_nr.csv"),    emit: test_ppis
    tuple val(meta), path("train_nr.fasta"), emit: train_fasta
    tuple val(meta), path("val_nr.fasta"),   emit: val_fasta
    tuple val(meta), path("test_nr.fasta"),  emit: test_fasta
    tuple val(meta), path("*_mqc.tsv"),      emit: mqc

    script:
    def inst_arg = instances ? "--instances ${instances}" : ''
    """
    remove_redundant.py \\
        --ppis           ${orig_ppis} \\
        --train_ppis     ${train_ppis} \\
        --val_ppis       ${val_ppis} \\
        --test_ppis      ${test_ppis} \\
        --train_fasta    ${train_fasta} \\
        --val_fasta      ${val_fasta} \\
        --test_fasta     ${test_fasta} \\
        --sim_train_val  ${sim_train_val} \\
        --sim_train_test ${sim_train_test} \\
        --id             ${mqcRowLabel(meta)} \\
        ${inst_arg}
    """
}
