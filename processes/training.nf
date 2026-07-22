// Datasets requesting the same embedding_model share one embedding call over
// the union of their train/val/test sequences, avoiding recomputation for
// proteins in more than one dataset. stageAs auto-numbers fasta_files since
// every dataset's split fasta shares the same name (train_nr.fasta, etc.)
// and would otherwise collide; embed_sequences.py merges/dedupes them.
process EMBED_SEQUENCES {
    publishDir(path: { "${params.outdir}/_shared/embeddings" }, mode: 'copy', saveAs: { f -> "embeddings_${embedding_model}.npz" })
    tag "embed_${embedding_model}"
    label "process_gpu"

    input:
    tuple val(embedding_model), path(fasta_files, stageAs: "input_*")

    output:
    tuple val(embedding_model), path("embeddings.npz")

    script:
    """
    embed_sequences.py \\
        --fasta ${fasta_files} \\
        --model ${embedding_model}
    """
}

process TRAIN_CLASSIFIER {
    tag "${meta.id}"
    label 'error_retry'

    input:
    tuple val(meta), path(train_csv), path(val_csv), path(test_balanced_csv), path(test_realistic_csv), path(embeddings)

    output:
    tuple val(meta), path("classifier_metrics_*_mqc.tsv"), emit: mqc

    script:
    """
    train_classifier.py \\
        --train          ${train_csv} \\
        --val            ${val_csv} \\
        --test_balanced  ${test_balanced_csv} \\
        --test_realistic ${test_realistic_csv} \\
        --embeddings     ${embeddings} \\
        --seed           ${params.seed} \\
        --id             ${meta.id}
    """
}
