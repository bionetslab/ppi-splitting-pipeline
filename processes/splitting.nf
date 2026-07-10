process SORT_PPIS {
    tag "${meta.id}"

    input:
    tuple val(meta), path(ppis), path(fasta), path(partition), path(node_mapping)

    output:
    tuple val(meta), path("train.csv"),   emit: train_ppis
    tuple val(meta), path("val.csv"),     emit: val_ppis
    tuple val(meta), path("test.csv"),    emit: test_ppis
    tuple val(meta), path("train.fasta"), emit: train_fasta
    tuple val(meta), path("val.fasta"),   emit: val_fasta
    tuple val(meta), path("test.fasta"),  emit: test_fasta
    tuple val(meta), path("*_mqc.tsv"),   emit: mqc, optional: true

    script:
    """
    sort_ppis.py \\
        --ppis         ${ppis} \\
        --partition    ${partition} \\
        --fasta        ${fasta} \\
        --node_mapping ${node_mapping}
    """
}

// Deliberately naive baseline: shuffles PPIs randomly instead of using a
// KaHIP partition, so the same protein can (and typically does) land in
// more than one split -- see bin/sort_ppis_random.py and
// bin/bias_analysis.py's "topology_shortcut" attribute. No redundancy
// removal runs downstream of this in SPLIT_POSITIVES: CD-HIT would treat a
// protein shared between train and test as trivially self-similar and
// strip it back out, defeating the point of this baseline.
process SPLIT_RANDOM {
    tag "${meta.id}"

    input:
    tuple val(meta), path(ppis), path(fasta)

    output:
    tuple val(meta), path("train.csv"),   emit: train_ppis
    tuple val(meta), path("val.csv"),     emit: val_ppis
    tuple val(meta), path("test.csv"),    emit: test_ppis
    tuple val(meta), path("train.fasta"), emit: train_fasta
    tuple val(meta), path("val.fasta"),   emit: val_fasta
    tuple val(meta), path("test.fasta"),  emit: test_fasta
    tuple val(meta), path("*_mqc.tsv"),   emit: mqc

    script:
    """
    sort_ppis_random.py \\
        --ppis        ${ppis} \\
        --fasta       ${fasta} \\
        --train-split ${meta.train_split} \\
        --val-split   ${meta.val_split} \\
        --test-split  ${meta.test_split} \\
        --seed        ${params.seed} \\
        --id          ${meta.id}
    """
}

process CDHIT2D {
    tag "${meta.id}_${label}"

    input:
    tuple val(meta), val(label), path(db1_fasta), path(db2_fasta)  // label: "train_val" | "train_test"

    output:
    tuple val(meta), val(label), path("cdhit.out"), emit: sim

    script:
    """
    cd-hit-2d \\
        -i  ${db1_fasta} \\
        -i2 ${db2_fasta} \\
        -o  cdhit.out \\
        -c  ${meta.cdhit_identity} \\
        -n  ${meta.cdhit_wordsize} \\
        -T  ${task.cpus} \\
        -M  4000
    """
}

process SOLVE_ILP {
    tag "${meta.id}"

    input:
    tuple val(meta), path(ppis), path(fasta), path(partition), path(node_mapping)
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
    """
    ${license_export}
    solve_ilp.py \\
        --ppis          ${ppis} \\
        --fasta         ${fasta} \\
        --partition     ${partition} \\
        --node_mapping  ${node_mapping} \\
        --train-split ${meta.train_split} \\
        --val-split   ${meta.val_split} \\
        --test-split  ${meta.test_split} \\
        --epsilon     ${meta.ilp_epsilon} \\
        --max-sec     ${params.ilp_max_sec} \\
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
          path(sim_train_test, stageAs: 'sim_train_test.out')

    output:
    tuple val(meta), path("train_nr.csv"),   emit: train_ppis
    tuple val(meta), path("val_nr.csv"),     emit: val_ppis
    tuple val(meta), path("test_nr.csv"),    emit: test_ppis
    tuple val(meta), path("train_nr.fasta"), emit: train_fasta
    tuple val(meta), path("val_nr.fasta"),   emit: val_fasta
    tuple val(meta), path("test_nr.fasta"),  emit: test_fasta
    tuple val(meta), path("*_mqc.tsv"),      emit: mqc

    script:
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
        --id             ${meta.id}
    """
}
