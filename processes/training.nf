// Datasets that request the same embedding_model share one embedding call
// over the union of their train/val/test sequences (ESM2/ProtT5 embedding
// is the most expensive per-protein step, so this avoids recomputing it for
// proteins that appear in more than one dataset). fasta_files is the flat
// list of every dataset-in-the-group's train/val/test fasta; stageAs
// auto-numbers them since every dataset's split fasta is literally named
// train_nr.fasta/val_nr.fasta/test_nr.fasta and would otherwise collide
// when staged together. embed_sequences.py already merges/dedupes across
// however many fasta files it's given.
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
