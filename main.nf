#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.ppis                   = "ppis.csv"
params.outdir                 = "results"
params.edge_weight            = "normalized_bitscore"  // "bitscore" or "normalized_bitscore"
params.kahip_seed             = 1234
params.kahip_k                = 3
params.kahip_preconfiguration = "strong"
params.cdhit_identity         = 0.4
params.cdhit_wordsize         = 2

include {
    FETCH_SEQUENCES
    GET_LENGTHS
    RUN_BLAST
    MAKE_METIS
    RUN_KAHIP
    SORT_PPIS
    CDHIT as CDHIT_TRAIN_VAL
    CDHIT as CDHIT_TRAIN_TEST
    REMOVE_REDUNDANT
    SAMPLE_NEGATIVES
    MULTIQC
} from './modules/processes'

workflow {
    ppis_ch = Channel.value(file(params.ppis, checkIfExists: true))

    sequences = FETCH_SEQUENCES(ppis_ch)
    lengths   = GET_LENGTHS(sequences)
    blast_out = RUN_BLAST(sequences)
    metis_out = MAKE_METIS(blast_out, lengths)
    partition = RUN_KAHIP(metis_out.graph)

    sorted = SORT_PPIS(
        ppis_ch,
        partition,
        sequences,
        metis_out.node_mapping
    )

    sim_tv = CDHIT_TRAIN_VAL(sorted.train_fasta, sorted.val_fasta)
    sim_tt = CDHIT_TRAIN_TEST(sorted.train_fasta, sorted.test_fasta)

    nr = REMOVE_REDUNDANT(
        sorted.train_ppis,
        sorted.val_ppis,
        sorted.test_ppis,
        sorted.train_fasta,
        sorted.val_fasta,
        sorted.test_fasta,
        sim_tv,
        sim_tt
    )

    SAMPLE_NEGATIVES(nr.train_ppis, nr.val_ppis, nr.test_ppis)

    mqc_files = sorted.mqc
        .mix(nr.mqc)
        .mix(SAMPLE_NEGATIVES.out.mqc)
        .collect()

    MULTIQC(mqc_files)
}
