include { RUN_BLAST; MAKE_METIS; RUN_KAHIP } from '../processes/clustering'

// Builds the protein similarity graph (BLAST all-vs-all -> METIS graph) and
// partitions it with KaHIP, ready for SPLIT_POSITIVES. Each dataset can
// skip BLAST via the samplesheet's blast_results column. BLAST still runs
// for "random"-split datasets (bias diagnostics/QC need it), but METIS/KaHIP
// don't -- SPLIT_POSITIVES's random branch never looks at the partition.
workflow CLUSTERING {
    take:
    sequences_ch      // tuple(meta, fasta)
    lengths_ch        // tuple(meta, lengths)
    blast_results_ch  // tuple(meta, blast_results_or_[])
    instances_ch      // tuple(meta, instances_or_[])

    main:
    branched = sequences_ch.join(blast_results_ch).branch { meta, fasta, blast_results ->
        precomputed: blast_results
            return tuple(meta, blast_results)
        needs_blast: true
            return tuple(meta, fasta)
    }

    blast_out = RUN_BLAST(branched.needs_blast).mix(branched.precomputed)

    // In DDI mode instances.tsv contracts the graph to one node per Pfam clan.
    // "random" datasets skip METIS/KaHIP and get (meta, []) placeholders below
    // so SPLIT_POSITIVES's partition/node_mapping joins still pair them.
    metis_branched = blast_out.join(lengths_ch).join(instances_ch).branch { meta, results, lengths, instances ->
        skip_partition: meta.split_method == "random"
            return meta
        needs_partition: true
            return tuple(meta, results, lengths, instances)
    }

    metis_out = MAKE_METIS(metis_branched.needs_partition)

    // The ILP splitter clusters proteins into many small KaHIP partitions
    // first, whereas the default splitter partitions straight into train/val/test.
    // k is clamped to the node count: a larger k leaves empty KaHIP blocks and
    // makes the ILP infeasible. The clan graph in DDI mode is often that small.
    kahip_inputs = metis_out.graph.map { meta, graph ->
        def k      = (meta.split_method == "ilp") ? meta.ilp_kahip_k : meta.kahip_k
        def nodes  = graph.withReader { reader -> reader.readLine() }.trim().split(/\s+/)[0] as int
        if (k > nodes) {
            log.warn "${meta.id}: similarity.graph has only ${nodes} nodes -- clamping k from ${k} to ${nodes}"
            k = nodes
        }
        tuple(meta, graph, k)
    }
    kahip_out = RUN_KAHIP(kahip_inputs)

    skipped = metis_branched.skip_partition.map { meta -> tuple(meta, []) }

    emit:
    blast_out    = blast_out
    node_mapping = metis_out.node_mapping.mix(skipped)
    partition    = kahip_out.mix(skipped)
}
