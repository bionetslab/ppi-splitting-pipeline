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

include {
    FETCH_DATA
    GET_LENGTHS
    RUN_BLAST
    MAKE_METIS
    RUN_KAHIP
    SORT_PPIS
    CDHIT2D
    SOLVE_ILP
    REMOVE_REDUNDANT
    SAMPLE_NEGATIVES
    SAMPLE_NEGATIVES_ILP
    EMBED_SEQUENCES
    TRAIN_CLASSIFIER
    BIAS_ANALYSIS
    COLLECT_BIAS
    SIMILARITY_HEATMAP
    MULTIQC
} from './modules/processes'

workflow {
    ppis_ch = channel.value(file(params.ppis, checkIfExists: true))

    if (params.sequences && params.go_annotations && params.species) {
        sequences_ch      = channel.value(file(params.sequences,      checkIfExists: true))
        go_annotations_ch = channel.value(file(params.go_annotations, checkIfExists: true))
        species_ch        = channel.value(file(params.species,        checkIfExists: true))
    } else {
        fetched           = FETCH_DATA(ppis_ch)
        sequences_ch      = fetched.sequences
        go_annotations_ch = fetched.go_annotations
        species_ch        = fetched.species
    }

    if (params.blast_results) {
        blast_out = channel.value(file(params.blast_results, checkIfExists: true))
    } else {
        blast_out = RUN_BLAST(sequences_ch)
    }

    lengths   = GET_LENGTHS(sequences_ch)
    metis_out = MAKE_METIS(blast_out, lengths)

    if (params.split_method == "ilp") {
        gurobi_license_ch = params.gurobi_license
            ? channel.value(file(params.gurobi_license, checkIfExists: true))
            : channel.value([])
        partition = RUN_KAHIP(metis_out.graph, params.ilp_kahip_k)
        sorted    = SOLVE_ILP(ppis_ch, sequences_ch, partition, metis_out.node_mapping, gurobi_license_ch)
    } else {
        partition = RUN_KAHIP(metis_out.graph, params.kahip_k)
        sorted    = SORT_PPIS(ppis_ch, partition, sequences_ch, metis_out.node_mapping)
    }

    // One channel, one process: each (label, fasta1, fasta2) item becomes its
    // own task, so Nextflow runs both CD-HIT-2D comparisons in parallel.
    cdhit_inputs_ch = sorted.train_fasta.combine(sorted.val_fasta).map { t, v -> tuple("train_val", t, v) }
        .mix(sorted.train_fasta.combine(sorted.test_fasta).map { t, te -> tuple("train_test", t, te) })

    cdhit_out = CDHIT2D(cdhit_inputs_ch)

    cdhit_branched = cdhit_out.sim.branch {
        label, f ->
            train_val:  label == "train_val"
            train_test: label == "train_test"
    }
    sim_tv = cdhit_branched.train_val.map { label, f -> f }
    sim_tt = cdhit_branched.train_test.map { label, f -> f }

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

    if (params.negative_sampling_method == "ilp") {
        candidate_network_ch = params.candidate_network
            ? channel.value(file(params.candidate_network, checkIfExists: true))
            : channel.value([])
        neg_gurobi_license_ch = params.gurobi_license
            ? channel.value(file(params.gurobi_license, checkIfExists: true))
            : channel.value([])

        // One channel, one process: each item (label, positives, neg_ratio)
        // becomes its own task, so Nextflow runs all four splits in parallel.
        neg_splits_ch = nr.train_ppis.map { f -> tuple("train", f, 1.0) }
            .mix(nr.val_ppis.map  { f -> tuple("val", f, 1.0) })
            .mix(nr.test_ppis.map { f -> tuple("test_balanced", f, 1.0) })
            .mix(nr.test_ppis.map { f -> tuple("test_realistic", f, 10.0) })

        neg_out = SAMPLE_NEGATIVES_ILP(
            neg_splits_ch, species_ch, go_annotations_ch, candidate_network_ch, neg_gurobi_license_ch
        )

        neg_branched = neg_out.labelled.branch {
            label, f ->
                train:          label == "train"
                val:            label == "val"
                test_balanced:  label == "test_balanced"
                test_realistic: label == "test_realistic"
        }

        neg = [
            train:          neg_branched.train.map { label, f -> f },
            val:            neg_branched.val.map { label, f -> f },
            test_balanced:  neg_branched.test_balanced.map { label, f -> f },
            test_realistic: neg_branched.test_realistic.map { label, f -> f },
            mqc: neg_out.mqc,
        ]
    } else {
        // One channel, one process: each item (label, positives, ratio, uniform)
        // becomes its own task, so Nextflow runs all four splits in parallel.
        neg_splits_ch = nr.train_ppis.map { f -> tuple("train", f, 1.0, false) }
            .mix(nr.val_ppis.map  { f -> tuple("val", f, 1.0, false) })
            .mix(nr.test_ppis.map { f -> tuple("test_balanced", f, 1.0, false) })
            .mix(nr.test_ppis.map { f -> tuple("test_realistic", f, 10.0, true) })

        neg_out = SAMPLE_NEGATIVES(neg_splits_ch)

        neg_branched = neg_out.labelled.branch {
            label, f ->
                train:          label == "train"
                val:            label == "val"
                test_balanced:  label == "test_balanced"
                test_realistic: label == "test_realistic"
        }

        neg = [
            train:          neg_branched.train.map { label, f -> f },
            val:            neg_branched.val.map { label, f -> f },
            test_balanced:  neg_branched.test_balanced.map { label, f -> f },
            test_realistic: neg_branched.test_realistic.map { label, f -> f },
            mqc: neg_out.mqc,
        ]
    }

    if (params.embedding_model in ["none", "esm2", "prot_t5"]) {
        embeddings = EMBED_SEQUENCES(nr.train_fasta, nr.val_fasta, nr.test_fasta)
    } else {
        embeddings = channel.value(file(params.embedding_model, checkIfExists: true))
    }

    clf      = TRAIN_CLASSIFIER(neg.train, neg.val, neg.test_balanced, neg.test_realistic, embeddings)

    same_species_ch = species_ch
        .splitCsv(header: true, sep: '\t')
        .map    { row -> row.taxon_id }
        .collect()
        .map    { ids -> ids.unique() }
        .filter { ids -> ids.size() > 1 }
        .map    { "same_species" }

    bias = BIAS_ANALYSIS(
        channel.of("sequence_similarity", "embedding_similarity",
                   "functional_relatedness_BP", "functional_relatedness_MF",
                   "functional_relatedness_CC", "self_interactions")
               .mix(same_species_ch).collect(),
        neg.train,
        neg.val,
        neg.test_balanced,
        neg.test_realistic,
        blast_out,
        embeddings,
        go_annotations_ch,
        species_ch
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
