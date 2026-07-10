include { FETCH_DATA; GET_LENGTHS } from '../processes/data_prep'

// Fetches sequences/GO annotations/species from UniProt (unless already
// supplied per-dataset via the samplesheet's sequences/go_annotations/species
// columns) and computes per-protein sequence lengths for the downstream
// similarity graph. Each dataset independently decides whether to fetch.
workflow DATA_PREP {
    take:
    datasets_ch  // tuple(meta, ppis, sequences, go_annotations, species, blast_results, candidate_network)

    main:
    branched = datasets_ch.branch { meta, ppis, sequences, go_annotations, species, blast_results, candidate_network ->
        precomputed: sequences && go_annotations && species
            return tuple(meta, sequences, go_annotations, species)
        needs_fetch: true
            return tuple(meta, ppis)
    }

    fetched = FETCH_DATA(branched.needs_fetch)

    sequences_ch      = branched.precomputed.map { meta, sequences, go_annotations, species -> tuple(meta, sequences) }
        .mix(fetched.sequences)
    go_annotations_ch = branched.precomputed.map { meta, sequences, go_annotations, species -> tuple(meta, go_annotations) }
        .mix(fetched.go_annotations)
    species_ch        = branched.precomputed.map { meta, sequences, go_annotations, species -> tuple(meta, species) }
        .mix(fetched.species)

    lengths = GET_LENGTHS(sequences_ch)

    emit:
    sequences      = sequences_ch
    go_annotations = go_annotations_ch
    species        = species_ch
    lengths        = lengths
}
