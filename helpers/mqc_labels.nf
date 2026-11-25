// Shared MultiQC display-label helpers (functions only, safe for every layer).
// A separate module because processes/ cannot include subworkflows/sample_negatives.nf
// without an include cycle.

// The MultiQC display label for one (dataset row, negative set): meta.mqc_labels[negset]
// if set, else "<id><negsuffix>". Report only; published paths key off meta.id.
def mqcLabel(meta, negset) {
    def mapped = (meta.mqc_labels instanceof Map) ? meta.mqc_labels[negset] : null
    if (mapped) {
        return mapped.toString()
    }
    def multi = meta.negative_sampling_method.toString().tokenize(',').size() > 1
    // A GString never equals an equal String; these are compared and grouped on.
    return (multi ? "${meta.id}_${negset}" : "${meta.id}").toString()
}

// The label for an artefact of the whole row rather than one negative set: the
// first set's label, so it still sorts by params.mqc_order.
def mqcRowLabel(meta) {
    def names = meta.negative_sampling_method.toString().tokenize(',')*.trim().findAll { n -> n }
    def first = names[0]
    def mapped = (meta.mqc_labels instanceof Map) ? meta.mqc_labels[first] : null
    if (mapped) {
        return mapped.toString()
    }
    return (names.size() > 1 ? "${meta.id}_${first}" : "${meta.id}").toString()
}
