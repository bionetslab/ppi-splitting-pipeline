#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.ppis                   = "ppis.csv"
params.outdir                 = "results"
params.edge_weight            = "normalized_bitscore"  // "bitscore" or "normalized_bitscore"
params.kahip_seed             = 1234
params.kahip_k                = 3
params.ilp_kahip_k            = 100
params.kahip_preconfiguration = "strong"
params.cdhit_identity         = 0.4
params.cdhit_wordsize         = 2
params.embedding_model        = "esm2"  // "none" (one-hot), "esm2", "prot_t5", or path to pre-computed .npz
params.seed                   = 42

include { DATA_PREP }        from './workflows/data_prep'
include { CLUSTERING }       from './workflows/clustering'
include { SPLIT_POSITIVES }  from './workflows/split_positives'
include { SAMPLE_NEGATIVES } from './workflows/sample_negatives'
include { TRAIN_BASELINE }   from './workflows/train_baseline'
include { QC }               from './workflows/qc'

workflow {
    ppis_ch = channel.value(file(params.ppis, checkIfExists: true))

    data = DATA_PREP(ppis_ch)

    clustered = CLUSTERING(data.sequences, data.lengths)

    split = SPLIT_POSITIVES(ppis_ch, data.sequences, clustered.partition, clustered.node_mapping)

    neg = SAMPLE_NEGATIVES(
        split.train_ppis, split.val_ppis, split.test_ppis,
        data.species, data.go_annotations
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