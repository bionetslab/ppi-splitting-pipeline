include { EMBED_SEQUENCES; TRAIN_CLASSIFIER } from '../modules/training'

// Embeds train/val/test sequences (unless a precomputed .npz is supplied via
// params.embedding_model) and trains the baseline classifier on top of them.
workflow TRAIN_BASELINE {
    take:
    train_fasta
    val_fasta
    test_fasta
    train_csv
    val_csv
    test_balanced_csv
    test_realistic_csv

    main:
    if (params.embedding_model in ["none", "esm2", "prot_t5"]) {
        embeddings = EMBED_SEQUENCES(train_fasta, val_fasta, test_fasta)
    } else {
        embeddings = channel.value(file(params.embedding_model, checkIfExists: true))
    }

    clf = TRAIN_CLASSIFIER(train_csv, val_csv, test_balanced_csv, test_realistic_csv, embeddings)

    emit:
    embeddings = embeddings
    mqc        = clf.mqc
}