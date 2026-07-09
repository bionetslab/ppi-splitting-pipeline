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
