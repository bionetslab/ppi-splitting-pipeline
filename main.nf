#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { samplesheetToList } from 'plugin/nf-schema'

include { DATA_PREP }        from './subworkflows/data_prep'
include { CLUSTERING }       from './subworkflows/clustering'
include { SPLIT_POSITIVES }  from './subworkflows/split_positives'
include { SAMPLE_NEGATIVES } from './subworkflows/sample_negatives'
include { TRAIN_BASELINE }   from './subworkflows/train_baseline'
include { QC }               from './subworkflows/qc'

// samplesheetToList() represents a blank optional cell as [] (empty list),
// not null, regardless of the field's declared type -- so a plain `!= null`
// check isn't enough to detect "not given" for numeric fields where 0 is a
// legitimate override value.
def isGiven(v) {
    !(v == null || v == [])
}

// One row per PPI dataset. Anything left blank in the samplesheet falls
// back to the corresponding default in nextflow.config, so a single run
// can process several datasets in parallel, each with its own overrides.
def buildDatasetsChannel() {
    // samplesheetToList() returns each row as a plain positional list (not
    // a map) unless schema properties are marked "meta" -- this must match
    // assets/schema_input.json's `properties` order exactly.
    def fields = [
        "id", "ppis", "sequences", "go_annotations", "species", "blast_results", "candidate_network",
        "partition", "node_mapping",
        "embedding_model", "cdhit_identity", "cdhit_wordsize", "split_method", "edge_weight",
        "kahip_k", "ilp_kahip_k", "train_split", "val_split", "test_split", "ilp_epsilon", "ilp_max_sec",
        "negative_sampling_method",
        "neg_ilp_alpha_confidence", "neg_ilp_time_limit", "neg_ilp_alpha_bias", "neg_ilp_lambda_degree",
        "neg_ilp_lambda_taxon_pair", "neg_ilp_lambda_self_loop", "neg_ilp_lambda_jaccard",
    ]
    def rows = samplesheetToList(params.samplesheet, "${projectDir}/assets/schema_input.json")

    return channel.fromList(rows).map { rowList ->
        def row = [fields, rowList].transpose().collectEntries { k, v -> [(k): v] }

        // split_only skips FETCH_DATA/CLUSTERING/TRAIN_BASELINE/QC entirely,
        // so every one of those steps' precomputed-input escape hatches
        // becomes mandatory, and the split/negative-sampling method choice
        // is no longer per-dataset -- it's always the ILP path.
        if (params.split_only) {
            if (!(row.sequences && row.go_annotations && row.species && row.partition && row.node_mapping)) {
                error("--split_only requires every samplesheet row to supply sequences, go_annotations, species, partition, and node_mapping (row '${row.id}' is missing at least one).")
            }
        }

        def meta = [
            id                       : row.id,
            embedding_model          : isGiven(row.embedding_model)          ? row.embedding_model          : params.embedding_model,
            cdhit_identity           : isGiven(row.cdhit_identity)           ? row.cdhit_identity           : params.cdhit_identity,
            cdhit_wordsize           : isGiven(row.cdhit_wordsize)           ? row.cdhit_wordsize           : params.cdhit_wordsize,
            split_method             : params.split_only ? "ilp" : (isGiven(row.split_method)             ? row.split_method             : params.split_method),
            edge_weight              : isGiven(row.edge_weight)              ? row.edge_weight              : params.edge_weight,
            kahip_k                  : isGiven(row.kahip_k)                  ? row.kahip_k                  : params.kahip_k,
            ilp_kahip_k              : isGiven(row.ilp_kahip_k)              ? row.ilp_kahip_k              : params.ilp_kahip_k,
            train_split              : isGiven(row.train_split)              ? row.train_split              : params.train_split,
            val_split                : isGiven(row.val_split)                ? row.val_split                : params.val_split,
            test_split               : isGiven(row.test_split)               ? row.test_split               : params.test_split,
            ilp_epsilon              : isGiven(row.ilp_epsilon)              ? row.ilp_epsilon              : params.ilp_epsilon,
            ilp_max_sec              : isGiven(row.ilp_max_sec)              ? row.ilp_max_sec              : params.ilp_max_sec,
            negative_sampling_method : params.split_only ? "ilp" : (isGiven(row.negative_sampling_method) ? row.negative_sampling_method : params.negative_sampling_method),
            neg_ilp_time_limit       : isGiven(row.neg_ilp_time_limit)       ? row.neg_ilp_time_limit       : params.neg_ilp_time_limit,
            neg_ilp_alpha_confidence : isGiven(row.neg_ilp_alpha_confidence)  ? row.neg_ilp_alpha_confidence  : params.neg_ilp_alpha_confidence,
            neg_ilp_alpha_bias       : isGiven(row.neg_ilp_alpha_bias)        ? row.neg_ilp_alpha_bias        : params.neg_ilp_alpha_bias,
            neg_ilp_lambda_degree    : isGiven(row.neg_ilp_lambda_degree)     ? row.neg_ilp_lambda_degree     : params.neg_ilp_lambda_degree,
            neg_ilp_lambda_taxon_pair: isGiven(row.neg_ilp_lambda_taxon_pair) ? row.neg_ilp_lambda_taxon_pair : params.neg_ilp_lambda_taxon_pair,
            neg_ilp_lambda_self_loop : isGiven(row.neg_ilp_lambda_self_loop)  ? row.neg_ilp_lambda_self_loop  : params.neg_ilp_lambda_self_loop,
            neg_ilp_lambda_jaccard   : isGiven(row.neg_ilp_lambda_jaccard)    ? row.neg_ilp_lambda_jaccard    : params.neg_ilp_lambda_jaccard,
        ]
        tuple(meta,
            file(row.ppis, checkIfExists: true),
            row.sequences         ? file(row.sequences,         checkIfExists: true) : [],
            row.go_annotations    ? file(row.go_annotations,    checkIfExists: true) : [],
            row.species           ? file(row.species,           checkIfExists: true) : [],
            row.blast_results     ? file(row.blast_results,     checkIfExists: true) : [],
            row.candidate_network ? file(row.candidate_network, checkIfExists: true) : [],
            row.partition         ? file(row.partition,         checkIfExists: true) : [],
            row.node_mapping      ? file(row.node_mapping,      checkIfExists: true) : [],
        )
    }
}

workflow {
    datasets_ch = buildDatasetsChannel()

    ppis_ch = datasets_ch.map { meta, ppis, sequences, go_annotations, species, blast_results, candidate_network, partition, node_mapping -> tuple(meta, ppis) }

    data = DATA_PREP(
        datasets_ch.map { meta, ppis, sequences, go_annotations, species, blast_results, candidate_network, partition, node_mapping ->
            tuple(meta, ppis, sequences, go_annotations, species, blast_results, candidate_network)
        }
    )

    if (params.split_only) {
        // --split_only: partition/node_mapping are precomputed and required
        // (validated in buildDatasetsChannel), so CLUSTERING (FETCH_DATA/
        // RUN_BLAST/MAKE_METIS/RUN_KAHIP) never needs to run at all.
        partition_ch    = datasets_ch.map { meta, ppis, sequences, go_annotations, species, blast_results, candidate_network, partition, node_mapping -> tuple(meta, partition) }
        node_mapping_ch = datasets_ch.map { meta, ppis, sequences, go_annotations, species, blast_results, candidate_network, partition, node_mapping -> tuple(meta, node_mapping) }
    } else {
        clustered = CLUSTERING(
            data.sequences, data.lengths,
            datasets_ch.map { meta, ppis, sequences, go_annotations, species, blast_results, candidate_network, partition, node_mapping -> tuple(meta, blast_results) }
        )
        partition_ch    = clustered.partition
        node_mapping_ch = clustered.node_mapping
    }

    split = SPLIT_POSITIVES(ppis_ch, data.sequences, partition_ch, node_mapping_ch)

    neg = SAMPLE_NEGATIVES(
        split.train_ppis, split.val_ppis, split.test_ppis,
        data.species, data.go_annotations,
        datasets_ch.map { meta, ppis, sequences, go_annotations, species, blast_results, candidate_network, partition, node_mapping -> tuple(meta, candidate_network) }
    )

    // --split_only stops here: SOLVE_ILP (via SPLIT_POSITIVES) + CDHIT2D +
    // REMOVE_REDUNDANT + SAMPLE_NEGATIVES_ILP have already produced and
    // published the four split files; TRAIN_BASELINE/QC add nothing this
    // mode asks for.
    if (!params.split_only) {
        baseline = TRAIN_BASELINE(
            split.train_fasta, split.val_fasta, split.test_fasta,
            neg.train, neg.val, neg.test_balanced, neg.test_realistic
        )

        QC(
            neg.train, neg.val, neg.test_balanced, neg.test_realistic,
            clustered.blast_out, baseline.embeddings, data.go_annotations, data.species,
            split.train_fasta, split.val_fasta, split.test_fasta,
            split.sorted_mqc, split.nr_mqc, neg.mqc, baseline.mqc
        )
    }
}
