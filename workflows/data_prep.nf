include { FETCH_DATA; GET_LENGTHS } from '../modules/data_prep'

// Fetches sequences/GO annotations/species from UniProt (unless already
// supplied via params.sequences/go_annotations/species) and computes
// per-protein sequence lengths for the downstream similarity graph.
workflow DATA_PREP {
    take:
    ppis_ch

    main:
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

    lengths = GET_LENGTHS(sequences_ch)

    emit:
    sequences      = sequences_ch
    go_annotations = go_annotations_ch
    species        = species_ch
    lengths        = lengths
}