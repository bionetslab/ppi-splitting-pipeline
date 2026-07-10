#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { samplesheetToList } from 'plugin/nf-schema'

include { DATA_PREP }        from './workflows/data_prep'
include { CLUSTERING }       from './workflows/clustering'
include { SPLIT_POSITIVES }  from './workflows/split_positives'
include { SAMPLE_NEGATIVES } from './workflows/sample_negatives'
include { TRAIN_BASELINE }   from './workflows/train_baseline'
include { QC }               from './workflows/qc'

// One row per PPI dataset. Anything left blank in the samplesheet falls
// back to the corresponding default in nextflow.config, so a single run
// can process several datasets in parallel, each with its own overrides.
def buildDatasetsChannel() {
    def rows = samplesheetToList(params.samplesheet, "${projectDir}/assets/schema_input.json")

    return channel.fromList(rows).map { row ->
        def meta = [
            id                       : row.id,
            embedding_model          : row.embedding_model          ?: params.embedding_model,
            cdhit_identity           : row.cdhit_identity           != null ? row.cdhit_identity : params.cdhit_identity,
            cdhit_wordsize           : row.cdhit_wordsize           != null ? row.cdhit_wordsize : params.cdhit_wordsize,
            split_method             : row.split_method             ?: params.split_method,
            edge_weight              : row.edge_weight              ?: params.edge_weight,
            kahip_k                  : row.kahip_k                  != null ? row.kahip_k : params.kahip_k,
            ilp_kahip_k              : row.ilp_kahip_k              != null ? row.ilp_kahip_k : params.ilp_kahip_k,
            train_split              : row.train_split              != null ? row.train_split : params.train_split,
            val_split                : row.val_split                != null ? row.val_split : params.val_split,
            test_split               : row.test_split               != null ? row.test_split : params.test_split,
            ilp_epsilon              : row.ilp_epsilon              != null ? row.ilp_epsilon : params.ilp_epsilon,
            negative_sampling_method : row.negative_sampling_method ?: params.negative_sampling_method,
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
