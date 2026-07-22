include { FETCH_DATA; GET_LENGTHS; GET_LENGTHS as GET_LENGTHS_SHARED; SUBSET_FETCHED_DATA } from '../processes/data_prep'

// Fetches sequences/GO annotations/species from UniProt (unless already
// supplied per-dataset via the samplesheet) and computes per-protein
// sequence lengths for the downstream similarity graph.
//
// Datasets needing a fetch are deduplicated into one union of protein IDs,
// fetched once, then split back out per dataset (SUBSET_FETCHED_DATA).
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

    // Precomputed datasets: use the supplied files as-is; each still needs
    // its own GET_LENGTHS since a precomputed sequences.fasta differs
    // dataset to dataset.
    precomputed_sequences = branched.precomputed.map { meta, sequences, go_annotations, species -> tuple(meta, sequences) }
    precomputed_go        = branched.precomputed.map { meta, sequences, go_annotations, species -> tuple(meta, go_annotations) }
    precomputed_species   = branched.precomputed.map { meta, sequences, go_annotations, species -> tuple(meta, species) }
    precomputed_lengths   = GET_LENGTHS(precomputed_sequences)

    // Extract protein1/protein2 from every needs-fetch dataset's PPI CSV,
    // dedupe, and fetch the union once under a synthetic "_shared" meta
    // (publishes to results/_shared/data/).
    proteins_list = branched.needs_fetch
        .flatMap { meta, ppis -> ppis.splitCsv(header: true).collectMany { row -> [row.protein1.trim(), row.protein2.trim()] } }
        .unique()
        .collectFile(name: 'proteins.txt', newLine: true, sort: true)

    shared_fetch   = FETCH_DATA(proteins_list.map { proteins -> tuple([id: "_shared"], proteins) })
    shared_lengths = GET_LENGTHS_SHARED(shared_fetch.sequences)

    subset_out = SUBSET_FETCHED_DATA(
        branched.needs_fetch,
        shared_fetch.sequences.map      { meta, f -> f }.first(),
        shared_fetch.go_annotations.map { meta, f -> f }.first(),
        shared_fetch.species.map        { meta, f -> f }.first(),
        shared_lengths.map               { meta, f -> f }.first(),
    )

    sequences_ch      = precomputed_sequences.mix(subset_out.sequences)
    go_annotations_ch = precomputed_go.mix(subset_out.go_annotations)
    species_ch        = precomputed_species.mix(subset_out.species)
    lengths_ch        = precomputed_lengths.mix(subset_out.lengths)

    emit:
    sequences      = sequences_ch
    go_annotations = go_annotations_ch
    species        = species_ch
    lengths        = lengths_ch
}
