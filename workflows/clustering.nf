include { RUN_BLAST; MAKE_METIS; RUN_KAHIP } from '../modules/clustering'

// Builds the protein similarity graph (BLAST all-vs-all -> METIS graph) and
// partitions it with KaHIP, ready for SPLIT_POSITIVES to assign train/val/test.
// Each dataset independently decides whether to skip BLAST (pre-supplied via
// the samplesheet's blast_results column) and how many KaHIP partitions to use.
workflow CLUSTERING {
    take:
    sequences_ch      // tuple(meta, fasta)
    lengths_ch        // tuple(meta, lengths)
    blast_results_ch  // tuple(meta, blast_results_or_[])

    main:
    branched = sequences_ch.join(blast_results_ch).branch { meta, fasta, blast_results ->
        precomputed: blast_results
            return tuple(meta, blast_results)
        needs_blast: true
            return tuple(meta, fasta)
    }

    blast_out = RUN_BLAST(branched.needs_blast).mix(branched.precomputed)

    metis_out = MAKE_METIS(blast_out.join(lengths_ch))

    // The ILP splitter clusters proteins into many small KaHIP partitions
    // first, whereas the default splitter partitions straight into train/val/test.
    kahip_inputs = metis_out.graph.map { meta, graph ->
        def k = (meta.split_method == "ilp") ? meta.ilp_kahip_k : meta.kahip_k
        tuple(meta, graph, k)
    }
    partition = RUN_KAHIP(kahip_inputs)

    emit:
    blast_out    = blast_out
    node_mapping = metis_out.node_mapping
    partition    = partition
}
