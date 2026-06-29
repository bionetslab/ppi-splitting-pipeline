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
params.embedding_model        = "esm2"  // "none" (one-hot), "esm2", "prot_t5", or path to pre-computed .npz
params.seed                   = 42

include {
    FETCH_DATA
    GET_LENGTHS
    RUN_BLAST
    MAKE_METIS
    RUN_KAHIP
    SORT_PPIS
    CDHIT as CDHIT_TRAIN_VAL
    CDHIT as CDHIT_TRAIN_TEST
    REMOVE_REDUNDANT
    SAMPLE_NEGATIVES
    EMBED_SEQUENCES
    TRAIN_CLASSIFIER
    BIAS_ANALYSIS
    COLLECT_BIAS
    SIMILARITY_HEATMAP
    MULTIQC
} from './modules/processes'

workflow {
    ppis_ch = Channel.value(file(params.ppis, checkIfExists: true))

    fetched        = FETCH_DATA(ppis_ch)
    lengths        = GET_LENGTHS(fetched.sequences)
    blast_out = RUN_BLAST(fetched.sequences)
    metis_out = MAKE_METIS(blast_out, lengths)
    partition = RUN_KAHIP(metis_out.graph)

    sorted = SORT_PPIS(
        ppis_ch,
        partition,
        fetched.sequences,
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

    neg = SAMPLE_NEGATIVES(nr.train_ppis, nr.val_ppis, nr.test_ppis)

    if (params.embedding_model in ["none", "esm2", "prot_t5"]) {
        embeddings = EMBED_SEQUENCES(nr.train_fasta, nr.val_fasta, nr.test_fasta)
    } else {
        embeddings = Channel.value(file(params.embedding_model, checkIfExists: true))
    }

    clf      = TRAIN_CLASSIFIER(neg.train, neg.val, neg.test_balanced, neg.test_realistic, embeddings)

    bias_attributes = [
        "sequence_similarity",
        "embedding_similarity",
        "functional_relatedness_BP",
        "functional_relatedness_MF",
        "functional_relatedness_CC",
        "self_interactions",
    ]
    bias = BIAS_ANALYSIS(
        bias_attributes,
        neg.train.first(),
        neg.val.first(),
        neg.test_balanced.first(),
        neg.test_realistic.first(),
        blast_out.first(),
        embeddings.first(),
        fetched.go_annotations.first()
    )
    scatter  = COLLECT_BIAS(bias.mqc.collect())
    heatmap  = SIMILARITY_HEATMAP(nr.train_fasta, nr.val_fasta, nr.test_fasta, blast_out)

    mqc_files = sorted.mqc
        .mix(nr.mqc)
        .mix(neg.mqc)
        .mix(clf.mqc)
        .mix(bias.mqc)
        .mix(scatter.mqc)
        .mix(heatmap)
        .collect()

    MULTIQC(mqc_files)
}
