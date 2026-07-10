include { EMBED_SEQUENCES; TRAIN_CLASSIFIER } from '../processes/training'

// Embeds train/val/test sequences (unless a precomputed .npz is supplied
// per-dataset via meta.embedding_model) and trains the baseline classifier
// on top of them. Datasets sharing the same embedding_model are embedded
// together in one EMBED_SEQUENCES call -- see that process for why.
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

    grouped_by_model = branched.compute
        .map { meta, train, val, test -> tuple(meta.embedding_model, meta, [train, val, test]) }
        .groupTuple()
    // tuple(embedding_model, [meta, meta, ...], [[train,val,test], [train,val,test], ...])

    embedded = EMBED_SEQUENCES(grouped_by_model.map { model, metas, triples -> tuple(model, triples.flatten()) })

    // Broadcast each model-group's one shared embeddings.npz back out to
    // every dataset that requested that model.
    computed = grouped_by_model
        .map { model, metas, triples -> tuple(model, metas) }
        .join(embedded)
        .flatMap { model, metas, emb -> metas.collect { meta -> tuple(meta, emb) } }

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
