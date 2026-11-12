include { ROBUSTNESS_SUMMARY; MULTIQC } from '../processes/qc'

// --robustness_test replaces BIAS_DIAGNOSTICS/QC: every dataset's classifier
// metrics from all (split, sampling, training) seed combinations are
// aggregated into one variance summary per dataset, and the MultiQC report
// shows only those summaries.
workflow ROBUSTNESS {
    take:
    clf_mqc  // tuple(meta, path or [paths]) -- one item per seed combination

    main:
    // TRAIN_CLASSIFIER emits one file per test set; flatten to one
    // (id, seeds, file) item per file before grouping per dataset.
    runs = clf_mqc.flatMap { meta, f ->
        def files = (f instanceof List) ? f : [f]
        files.collect { ff -> tuple(meta.id, [meta.split_seed, meta.sampling_seed, meta.training_seed], ff) }
    }.groupTuple()

    summary = ROBUSTNESS_SUMMARY(runs)

    MULTIQC(summary.mqc.map { id, f -> f }.collect())
}
