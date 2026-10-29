#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameter types only (defaults live in conf/params.config): Nextflow 26 passes CLI values
// as Strings, so `--split_only false` would be truthy. Only the entry script's block counts.
params {
    split_only: Boolean
    bias_only: Boolean
    seed: Integer
    heatmap_max_per_split: Integer

    // clustering
    cdhit_identity: Float
    cdhit_wordsize: Integer
    kahip_seed: Integer
    kahip_k: Integer
    ilp_kahip_k: Integer

    // positive split
    train_split: Float
    val_split: Float
    test_split: Float
    ilp_epsilon: Float
    ilp_max_sec: Integer

    // negative sampling
    neg_ilp_lambda_degree: Float
    neg_ilp_lambda_taxon_pair: Float
    neg_ilp_lambda_self_loop: Float
    neg_ilp_lambda_jaccard: Float
    neg_ilp_time_limit: Integer
    neg_ilp_mip_gap: Float

    // DDI mode
    ddi_mode: Boolean
    ddi_examples_target: Integer
    ddi_examples_pool_factor: Integer
}

include { samplesheetToList } from 'plugin/nf-schema'

include { DATA_PREP }        from './subworkflows/data_prep'
include { DATA_PREP_DDI }    from './subworkflows/data_prep_ddi'
include { CLUSTERING }       from './subworkflows/clustering'
include { SPLIT_POSITIVES }  from './subworkflows/split_positives'
include { SAMPLE_NEGATIVES } from './subworkflows/sample_negatives'
include { TRAIN_BASELINE }   from './subworkflows/train_baseline'
include { QC }               from './subworkflows/qc'
include { BIAS_DIAGNOSTICS } from './subworkflows/bias_diagnostics'

include { mqcLabel }        from './helpers/mqc_labels'

def isGiven(v) {
    !(v == null || v == [])
}

// val or test may be 0, but not both; train never.
def checkSplitFractions(meta) {
    def (train, val, test) = [meta.train_split, meta.val_split, meta.test_split].collect { it.toString().toBigDecimal() }
    if (train == 0 || (val == 0 && test == 0) || (train + val + test - 1).abs() > 1e-6) {
        error("Dataset '${meta.id}': train_split=${train}, val_split=${val}, test_split=${test} -- train_split must be > 0, val_split and test_split must not both be 0, and the three must sum to 1.")
    }
}

// Parse and validate one row's comma-separated negative-set list. Called from PPI_SPLITTING
// so embedding pipelines get the same cross-column checks.
def parseNegsets(meta, candidate_network) {
    def known = ["default", "ilp", "ilp_candidates", "uniform"]
    def names = meta.negative_sampling_method.toString().tokenize(',')*.trim().findAll { n -> n }
    if (!names) {
        error("Dataset '${meta.id}': negative_sampling_method is empty -- name at least one of ${known.join(', ')}.")
    }
    def unknown = names.findAll { n -> !(n in known) }
    if (unknown) {
        error("Dataset '${meta.id}': unknown negative_sampling_method ${unknown.join(', ')} (known: ${known.join(', ')}).")
    }
    if (names.size() != new HashSet(names).size()) {
        error("Dataset '${meta.id}': negative_sampling_method '${meta.negative_sampling_method}' names a method more than once.")
    }
    if ("ilp_candidates" in names && !candidate_network) {
        error("Dataset '${meta.id}': negative_sampling_method lists 'ilp_candidates' but the row supplies no candidate_network.")
    }
    return names
}

// Warn when meta.mqc_labels does not cover exactly the row's negative sets.
def checkMqcLabels(meta, negsets) {
    if (!(meta.mqc_labels instanceof Map)) {
        return
    }
    def given    = meta.mqc_labels.keySet().collect { k -> k.toString() }
    def expected = negsets.collect { n -> n.toString() }
    if (given.toSorted() != expected.toSorted()) {
        log.warn "Dataset '${meta.id}': mqc_labels covers ${given.toSorted().join(', ')} but the row asks for negative sets ${expected.toSorted().join(', ')} -- missing sets fall back to '<id>_<negset>'."
    }
}

// Resolve the run's report order from every display label it will produce.
// params.mqc_order may be a list (config) or a comma-separated string (CLI).
def resolveMqcOrder(labels) {
    def present   = labels.collect { l -> l.toString() }.unique()
    def requested = (params.mqc_order instanceof CharSequence)
        ? params.mqc_order.toString().tokenize(',')*.trim().findAll { r -> r }
        : (params.mqc_order ?: []).collect { r -> r.toString().trim() }.findAll { r -> r }

    def ordered = []
    requested.each { r ->
        if (!(r in present)) {
            log.warn "mqc_order names '${r}', which is not a display label of this run -- ignoring it. Labels present: ${present.toSorted().join(', ')}"
        }
        else if (r in ordered) {
            log.warn "mqc_order names '${r}' more than once -- keeping the first occurrence."
        }
        else {
            ordered << r
        }
    }
    // Labels the list did not mention, appended alphabetically.
    def leftover = (present - ordered).toSorted()
    if (requested && !ordered) {
        log.warn "mqc_order was given but matched none of this run's display labels, so the report order is alphabetical. Labels present: ${present.toSorted().join(', ')}"
    }
    return (ordered + leftover).join(',')
}

// One row per PPI dataset. Anything left blank in the samplesheet falls back to the corresponding default in
// conf/params.config, so a single run can process several datasets in parallel, each with its own overrides.
def buildDatasetsChannel() {
    // samplesheetToList() returns each row as a positional list, not a map.
    // Order here must match assets/schema_input.json's `properties`.
    def fields = [
        "id", "ppis", "sequences", "go_annotations", "species", "domain_instances", "blast_results", "candidate_network",
        "partition", "node_mapping",
        "train_ppis", "val_ppis", "test_balanced_ppis", "test_realistic_ppis",
        "embedding_model", "cdhit_identity", "cdhit_wordsize", "split_method", "edge_weight",
        "kahip_k", "ilp_kahip_k", "train_split", "val_split", "test_split", "ilp_epsilon", "ilp_max_sec",
        "negative_sampling_method",
        "neg_ilp_time_limit", "neg_ilp_lambda_degree",
        "neg_ilp_lambda_taxon_pair", "neg_ilp_lambda_self_loop", "neg_ilp_lambda_jaccard",
    ]
    // moduleDir, not projectDir: when included, projectDir is the embedding pipeline's.
    def rows = samplesheetToList(params.samplesheet, "${moduleDir}/assets/schema_input.json")

    return channel.fromList(rows).map { rowList ->
        def row = [fields, rowList].transpose().collectEntries { k, v -> [(k): v] }

        // Every mode except --bias_only reads the raw PPI CSV (DATA_PREP/
        // SPLIT_POSITIVES); --bias_only never touches it, so it's the one
        // mode where ppis may be left blank in the samplesheet.
        if (!params.bias_only && !isGiven(row.ppis)) {
            error("ppis is required for every samplesheet row unless --bias_only is set (row '${row.id}' is missing it).")
        }

        // bias_only skips everything upstream of the bias diagnostics, so the splits, BLAST hits,
        // GO/species tables, and embeddings must all be supplied precomputed. Mirrored in PPI_SPLITTING.
        if (params.bias_only) {
            def required = ["train_ppis", "val_ppis", "test_balanced_ppis", "test_realistic_ppis", "go_annotations", "species", "blast_results"]
            def missing = required.findAll { !row[it] }
            if (missing) {
                error("--bias_only requires every samplesheet row to supply ${required.join(', ')} (row '${row.id}' is missing ${missing.join(', ')}).")
            }
            if (!isGiven(row.embedding_model) || row.embedding_model in ["none", "esm2", "prot_t5"]) {
                error("--bias_only requires embedding_model to be a path to a precomputed .npz (row '${row.id}' has '${row.embedding_model}') -- there's no embedding step to compute it.")
            }
        }

        // split_only skips FETCH_DATA/CLUSTERING/TRAIN_BASELINE/QC entirely,
        // split/negative-sampling method choice is no longer per-dataset; it's always the ILP path.
        // Mirrored in PPI_SPLITTING; DDI mode needs domain_instances instead of GO annotations.
        if (params.split_only) {
            def required = params.ddi_mode
                ? ["sequences", "species", "domain_instances", "partition", "node_mapping"]
                : ["sequences", "go_annotations", "species", "partition", "node_mapping"]
            def missing = required.findAll { !row[it] }
            if (missing) {
                error("--split_only requires every samplesheet row to supply ${required.join(', ')} (row '${row.id}' is missing ${missing.join(', ')}).")
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
            // --split_only keeps ilp/ilp_candidates and forces anything else to ilp.
            negative_sampling_method : params.split_only
                ? ((isGiven(row.negative_sampling_method) && row.negative_sampling_method.toString().tokenize(',')*.trim().every { n -> n in ["ilp", "ilp_candidates"] }) ? row.negative_sampling_method : "ilp")
                : (isGiven(row.negative_sampling_method) ? row.negative_sampling_method : params.negative_sampling_method),
            neg_ilp_time_limit       : isGiven(row.neg_ilp_time_limit)       ? row.neg_ilp_time_limit       : params.neg_ilp_time_limit,
            neg_ilp_lambda_degree    : isGiven(row.neg_ilp_lambda_degree)     ? row.neg_ilp_lambda_degree     : params.neg_ilp_lambda_degree,
            neg_ilp_lambda_taxon_pair: isGiven(row.neg_ilp_lambda_taxon_pair) ? row.neg_ilp_lambda_taxon_pair : params.neg_ilp_lambda_taxon_pair,
            neg_ilp_lambda_self_loop : isGiven(row.neg_ilp_lambda_self_loop)  ? row.neg_ilp_lambda_self_loop  : params.neg_ilp_lambda_self_loop,
            neg_ilp_lambda_jaccard   : isGiven(row.neg_ilp_lambda_jaccard)    ? row.neg_ilp_lambda_jaccard    : params.neg_ilp_lambda_jaccard,
        ]
        // An absent optional file is [], never null: a `path` input accepts [] but not null.
        def files = [
            ppis             : row.ppis              ? file(row.ppis,              checkIfExists: true) : [],
            sequences        : row.sequences         ? file(row.sequences,         checkIfExists: true) : [],
            go_annotations   : row.go_annotations    ? file(row.go_annotations,    checkIfExists: true) : [],
            species          : row.species           ? file(row.species,           checkIfExists: true) : [],
            domain_instances : row.domain_instances  ? file(row.domain_instances,  checkIfExists: true) : [],
            blast_results    : row.blast_results     ? file(row.blast_results,     checkIfExists: true) : [],
            candidate_network: row.candidate_network ? file(row.candidate_network, checkIfExists: true) : [],
            partition        : row.partition         ? file(row.partition,         checkIfExists: true) : [],
            node_mapping     : row.node_mapping      ? file(row.node_mapping,      checkIfExists: true) : [],
            // --bias_only's four precomputed splits
            train_ppis         : row.train_ppis         ? file(row.train_ppis,         checkIfExists: true) : [],
            val_ppis           : row.val_ppis           ? file(row.val_ppis,           checkIfExists: true) : [],
            test_balanced_ppis : row.test_balanced_ppis ? file(row.test_balanced_ppis, checkIfExists: true) : [],
            test_realistic_ppis: row.test_realistic_ppis ? file(row.test_realistic_ppis, checkIfExists: true) : [],
        ]
        tuple(meta, files)
    }
}

// The whole pipeline, as a named workflow so another project can `include` it.
workflow PPI_SPLITTING {
    take:
    datasets_ch   // tuple(meta, filesMap); filesMap keys as in buildDatasetsChannel(), absent = [].
                  // meta may carry an optional mqc_labels map (negset -> MultiQC label).

    main:
    // Validate every caller's datasets here (embedding pipelines skip buildDatasetsChannel()).
    // Validate only, never rewrite meta: it is the join key.
    checked_ch = datasets_ch.map { meta, f ->
        if (!params.bias_only && !f.ppis) {
            error("Dataset '${meta.id}' supplies no ppis file, which every mode except --bias_only reads.")
        }
        if (!params.bias_only) {
            checkSplitFractions(meta)
        }
        // Same checks as buildDatasetsChannel(); change the two together.
        if (params.bias_only) {
            def required = ["train_ppis", "val_ppis", "test_balanced_ppis", "test_realistic_ppis", "go_annotations", "species", "blast_results"]
            def missing = required.findAll { k -> !f[k] }
            if (missing) {
                error("--bias_only requires every dataset to supply ${required.join(', ')} (dataset '${meta.id}' is missing ${missing.join(', ')}).")
            }
            if (!meta.embedding_model || meta.embedding_model.toString() in ["none", "esm2", "prot_t5"]) {
                error("--bias_only requires embedding_model to be a path to a precomputed .npz (dataset '${meta.id}' has '${meta.embedding_model}') -- there's no embedding step to compute it.")
            }
        }
        if (params.split_only) {
            def required = params.ddi_mode
                ? ["sequences", "species", "domain_instances", "partition", "node_mapping"]
                : ["sequences", "go_annotations", "species", "partition", "node_mapping"]
            def missing = required.findAll { k -> !f[k] }
            if (missing) {
                error("--split_only requires every dataset to supply ${required.join(', ')} (dataset '${meta.id}' is missing ${missing.join(', ')}).")
            }
            def neg_names = meta.negative_sampling_method.toString().tokenize(',')*.trim().findAll { n -> n }
            if (meta.split_method != "ilp" || !neg_names || !neg_names.every { n -> n in ["ilp", "ilp_candidates"] }) {
                error("--split_only runs the ILP path only, but dataset '${meta.id}' asks for split_method='${meta.split_method}', negative_sampling_method='${meta.negative_sampling_method}' (allowed: ilp, ilp_candidates or both).")
            }
        }
        tuple(meta, f)
    }

    ppis_ch = checked_ch.map { meta, f -> tuple(meta, f.ppis) }

    // Gate the candidate network on ilp_candidates once, here.
    prepared_ch = checked_ch.map { meta, f ->
        def negsets = parseNegsets(meta, f.candidate_network)
        if (params.bias_only && negsets.size() != 1) {
            error("Dataset '${meta.id}': --bias_only takes one negative set per row, but negative_sampling_method names ${negsets.size()} (${negsets.join(', ')}).")
        }
        checkMqcLabels(meta, negsets)
        if (f.candidate_network && !("ilp_candidates" in negsets)) {
            log.warn "${meta.id}: candidate_network supplied but negative_sampling_method does not list 'ilp_candidates' -- ignoring the network"
        }
        tuple(meta, negsets, ("ilp_candidates" in negsets) ? f.candidate_network : [])
    }
    // tuple(meta, [negset, ...]) -- the negative sets this row wants, in list order.
    negsets_ch = prepared_ch.map { meta, negsets, cand -> tuple(meta, negsets) }
    // tuple(meta, candidate_network or [])
    candidate_network_ch = prepared_ch.map { meta, negsets, cand -> tuple(meta, cand) }

    // One comma-joined report order string; only MULTIQC takes it, so it stays out of other task hashes.
    mqc_order_ch = prepared_ch
        .flatMap { meta, negsets, cand -> negsets.collect { ns -> mqcLabel(meta, ns) } }
        .collect()
        .map { labels -> resolveMqcOrder(labels) }

    if (params.bias_only) {
        // --bias_only: the four labelled splits come straight from the files map, one negative set per row.
        bias_in = checked_ch.join(negsets_ch).map { meta, f, negsets -> tuple(meta, negsets[0], f) }
        bias = BIAS_DIAGNOSTICS(
            bias_in.map { meta, ns, f -> tuple(meta, ns, f.train_ppis) },
            bias_in.map { meta, ns, f -> tuple(meta, ns, f.val_ppis) },
            bias_in.map { meta, ns, f -> tuple(meta, ns, f.test_balanced_ppis) },
            bias_in.map { meta, ns, f -> tuple(meta, ns, f.test_realistic_ppis) },
            checked_ch.map { meta, f -> tuple(meta, f.blast_results) },
            checked_ch.map { meta, f -> tuple(meta, file(meta.embedding_model, checkIfExists: true)) },
            checked_ch.map { meta, f -> tuple(meta, f.go_annotations) },
            checked_ch.map { meta, f -> tuple(meta, f.species) }
        )

        // Every emit is defined in both branches so its shape does not depend on the mode.
        instances_ch        = channel.empty()
        sequences_ch        = channel.empty()
        labelled_ch         = channel.empty()
        multiqc_report_ch   = channel.empty()
        dropped_families_ch = channel.empty()
    } else {
        // Both emit the same channel names, so nothing downstream branches on DDI mode.
        if (params.ddi_mode) {
            data = DATA_PREP_DDI(checked_ch)
        } else {
            data = DATA_PREP(checked_ch)
        }

        if (params.split_only) {
            // --split_only: partition/node_mapping are precomputed and required
            // (validated above), so CLUSTERING (FETCH_DATA/
            // RUN_BLAST/MAKE_METIS/RUN_KAHIP) never needs to run at all.
            partition_ch    = checked_ch.map { meta, f -> tuple(meta, f.partition) }
            node_mapping_ch = checked_ch.map { meta, f -> tuple(meta, f.node_mapping) }
        } else {
            clustered = CLUSTERING(
                data.sequences, data.lengths,
                checked_ch.map { meta, f -> tuple(meta, f.blast_results) }
            )
            partition_ch    = clustered.partition
            node_mapping_ch = clustered.node_mapping
        }

        split = SPLIT_POSITIVES(ppis_ch, data.sequences, partition_ch, node_mapping_ch)

        neg = SAMPLE_NEGATIVES(
            split.train_ppis, split.val_ppis, split.test_ppis,
            data.species, data.go_annotations,
            candidate_network_ch,
            negsets_ch
        )

        // --split_only stops here: SOLVE_ILP (via SPLIT_POSITIVES) + CDHIT2D +
        // REMOVE_REDUNDANT + SAMPLE_NEGATIVES_ILP have already produced and
        // published the four split files; TRAIN_BASELINE/QC add nothing this
        // mode asks for.
        multiqc_report_ch = channel.empty()
        if (!params.split_only) {
            baseline = TRAIN_BASELINE(
                split.train_fasta, split.val_fasta, split.test_fasta,
                neg.train, neg.val, neg.test_balanced, neg.test_realistic
            )

            bias = BIAS_DIAGNOSTICS(
                neg.train, neg.val, neg.test_balanced, neg.test_realistic,
                clustered.blast_out, baseline.embeddings, data.go_annotations, data.species
            )

            qc = QC(
                bias.scatter, clustered.blast_out,
                split.train_fasta, split.val_fasta, split.test_fasta,
                split.sorted_mqc, split.nr_mqc, neg.mqc, baseline.mqc,
                mqc_order_ch
            )
            multiqc_report_ch = qc.multiqc_report
        }
        instances_ch        = data.instances
        sequences_ch        = data.sequences
        labelled_ch         = neg.labelled
        dropped_families_ch = data.dropped_families
    }

    emit:
    // negset is a tuple field, not a meta key: join these on `by: [0, 1]`.
    // Under --bias_only every channel is empty.
    instances        = instances_ch       // tuple(meta, instances.tsv); (meta, []) in PPI mode
    sequences        = sequences_ch       // tuple(meta, sequences.fasta)
    labelled         = labelled_ch        // tuple(meta, negset, label, csv) at family level
    multiqc_report   = multiqc_report_ch  // path multiqc_report.html; empty under --split_only and --bias_only
    dropped_families = dropped_families_ch  // tuple([id: "_shared"], tsv); at most one item, never join() it
}

workflow {
    PPI_SPLITTING( buildDatasetsChannel() )
}
