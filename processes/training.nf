process EMBED_SEQUENCES {
    publishDir(path: { "${params.outdir}/${meta.id}/data" })
    tag "${meta.id}"
    label "process_gpu"

    input:
    tuple val(meta), path(train_fasta), path(val_fasta), path(test_fasta)

    output:
    tuple val(meta), path("embeddings.npz")

    script:
    """
    embed_sequences.py \\
        --fasta ${train_fasta} ${val_fasta} ${test_fasta} \\
        --model ${meta.embedding_model}
    """
}

process TRAIN_CLASSIFIER {
    tag "${meta.id}"

    input:
    tuple val(meta), path(train_csv), path(val_csv), path(test_balanced_csv), path(test_realistic_csv), path(embeddings)

    output:
    tuple val(meta), path("classifier_metrics_mqc.tsv"), emit: mqc

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
