include { RUN_BLAST; MAKE_METIS; RUN_KAHIP } from '../modules/clustering'

// Builds the protein similarity graph (BLAST all-vs-all -> METIS graph) and
// partitions it with KaHIP, ready for SPLIT_POSITIVES to assign train/val/test.
workflow CLUSTERING {
    take:
    sequences_ch
    lengths_ch

    main:
    if (params.blast_results) {
        blast_out = channel.value(file(params.blast_results, checkIfExists: true))
    } else {
        blast_out = RUN_BLAST(sequences_ch)
    }

    metis_out = MAKE_METIS(blast_out, lengths_ch)

    // The ILP splitter clusters proteins into many small KaHIP partitions
    // first, whereas the default splitter partitions straight into train/val/test.
    kahip_k   = (params.split_method == "ilp") ? params.ilp_kahip_k : params.kahip_k
    partition = RUN_KAHIP(metis_out.graph, kahip_k)

    emit:
    blast_out    = blast_out
    node_mapping = metis_out.node_mapping
    partition    = partition
}