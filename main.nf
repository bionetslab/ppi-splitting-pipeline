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
        "embedding_model", "cdhit_identity", "cdhit_wordsize", "split_method", "edge_weight",
        "kahip_k", "ilp_kahip_k", "train_split", "val_split", "test_split", "ilp_epsilon", "ilp_max_sec",
        "negative_sampling_method",
        "neg_ilp_alpha_confidence", "neg_ilp_time_limit", "neg_ilp_alpha_bias", "neg_ilp_lambda_degree",
        "neg_ilp_lambda_taxon_pair", "neg_ilp_lambda_self_loop", "neg_ilp_lambda_jaccard",
    ]
    def rows = samplesheetToList(params.samplesheet, "${projectDir}/assets/schema_input.json")

    return channel.fromList(rows).map { rowList ->
        def row = [fields, rowList].transpose().collectEntries { k, v -> [(k): v] }

        def meta = [
            id                       : row.id,
            embedding_model          : isGiven(row.embedding_model)          ? row.embedding_model          : params.embedding_model,
            cdhit_identity           : isGiven(row.cdhit_identity)           ? row.cdhit_identity           : params.cdhit_identity,
            cdhit_wordsize           : isGiven(row.cdhit_wordsize)           ? row.cdhit_wordsize           : params.cdhit_wordsize,
            split_method             : isGiven(row.split_method)             ? row.split_method             : params.split_method,
            edge_weight              : isGiven(row.edge_weight)              ? row.edge_weight              : params.edge_weight,
            kahip_k                  : isGiven(row.kahip_k)                  ? row.kahip_k                  : params.kahip_k,
            ilp_kahip_k              : isGiven(row.ilp_kahip_k)              ? row.ilp_kahip_k              : params.ilp_kahip_k,
            train_split              : isGiven(row.train_split)              ? row.train_split              : params.train_split,
            val_split                : isGiven(row.val_split)                ? row.val_split                : params.val_split,
            test_split               : isGiven(row.test_split)               ? row.test_split               : params.test_split,
            ilp_epsilon              : isGiven(row.ilp_epsilon)              ? row.ilp_epsilon              : params.ilp_epsilon,
            ilp_max_sec              : isGiven(row.ilp_max_sec)              ? row.ilp_max_sec              : params.ilp_max_sec,
            negative_sampling_method : isGiven(row.negative_sampling_method) ? row.negative_sampling_method : params.negative_sampling_method,
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
        )
    }
}

workflow {
    datasets_ch = buildDatasetsChannel()

    ppis_ch = datasets_ch.map { meta, ppis, sequences, go_annotations, species, blast_results, candidate_network -> tuple(meta, ppis) }

    data = DATA_PREP(datasets_ch)

    clustered = CLUSTERING(
        data.sequences, data.lengths,
        datasets_ch.map { meta, ppis, sequences, go_annotations, species, blast_results, candidate_network -> tuple(meta, blast_results) }
    )

    split = SPLIT_POSITIVES(ppis_ch, data.sequences, clustered.partition, clustered.node_mapping)

    neg = SAMPLE_NEGATIVES(
        split.train_ppis, split.val_ppis, split.test_ppis,
        data.species, data.go_annotations,
        datasets_ch.map { meta, ppis, sequences, go_annotations, species, blast_results, candidate_network -> tuple(meta, candidate_network) }
    )

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
