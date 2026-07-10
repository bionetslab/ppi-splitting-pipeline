include { EMBED_SEQUENCES; TRAIN_CLASSIFIER } from '../processes/training'

// Embeds train/val/test sequences (unless a precomputed .npz is supplied
// per-dataset via meta.embedding_model) and trains the baseline classifier
// on top of them.
workflow TRAIN_BASELINE {
    take:
    train_fasta        // tuple(meta, path)
    val_fasta            // tuple(meta, path)
    test_fasta            // tuple(meta, path)
    train_csv              // tuple(meta, path)
    val_csv                  // tuple(meta, path)
    test_balanced_csv          // tuple(meta, path)
    test_realistic_csv          // tuple(meta, path)

    main:
    branched = train_fasta.join(val_fasta).join(test_fasta).branch { meta, train, val, test ->
        compute:     meta.embedding_model in ["none", "esm2", "prot_t5"]
        precomputed: true
    }

    computed    = EMBED_SEQUENCES(branched.compute)
    precomputed = branched.precomputed.map { meta, train, val, test ->
        tuple(meta, file(meta.embedding_model, checkIfExists: true))
    }
    embeddings = computed.mix(precomputed)

    clf_inputs = train_csv.join(val_csv).join(test_balanced_csv).join(test_realistic_csv).join(embeddings)
    clf = TRAIN_CLASSIFIER(clf_inputs)

    emit:
    embeddings = embeddings
    mqc        = clf.mqc
}
