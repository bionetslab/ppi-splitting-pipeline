include { EMPTY_GO_ANNOTATIONS; FETCH_DOMAIN_META; GET_LENGTHS as GET_LENGTHS_DDI; GET_LENGTHS as GET_LENGTHS_SHARED_DDI; SUBSET_DOMAIN_DATA } from '../processes/data_prep'

// DATA_PREP's DDI-mode counterpart: turns the Pfam family accessions in the
// interaction file into sampled domain instances (sequences, species, lengths)
// plus instances.tsv, mapping every instance to its family, clan and parent
// protein. Families are fetched once under a shared meta, then split back out.
workflow DATA_PREP_DDI {
    take:
    datasets_ch  // tuple(meta, files): ppis, sequences, species, domain_instances

    main:
    branched = datasets_ch.branch { meta, f ->
        precomputed: f.sequences && f.species && f.domain_instances
            return tuple(meta, f.sequences, f.species, f.domain_instances)
        needs_fetch: true
            return tuple(meta, f.ppis)
    }

    precomputed_sequences = branched.precomputed.map { meta, sequences, species, domain_instances -> tuple(meta, sequences) }
    precomputed_species   = branched.precomputed.map { meta, sequences, species, domain_instances -> tuple(meta, species) }
    precomputed_instances = branched.precomputed.map { meta, sequences, species, domain_instances -> tuple(meta, domain_instances) }
    precomputed_lengths   = GET_LENGTHS_DDI(precomputed_sequences)
    precomputed_go        = EMPTY_GO_ANNOTATIONS(branched.precomputed.map { meta, sequences, species, domain_instances -> meta })

    // protein1/protein2 hold Pfam family accessions in DDI mode.
    families_list = branched.needs_fetch
        .flatMap { meta, ddis -> ddis.splitCsv(header: true).collectMany { row -> [row.protein1.trim(), row.protein2.trim()] } }
        .unique()
        .collectFile(name: 'families.txt', newLine: true, sort: true)

    // Optional local copies of the release-wide Pfam reference files.
    clans_ch = params.pfam_clans
        ? channel.value(file(params.pfam_clans, checkIfExists: true))
        : channel.value([])

    regions_ch = params.pfam_regions
        ? channel.value(file(params.pfam_regions, checkIfExists: true))
        : channel.value([])

    // UniProt flat files supplying sequence, taxon and review flag, Swiss-Prot
    // first; a List or a comma-separated string.
    def dat_raw = !params.uniprot_dat
        ? []
        : (params.uniprot_dat instanceof List ? params.uniprot_dat : params.uniprot_dat.toString().split(','))
    def dat_files = dat_raw.collect { it.toString().trim() }.findAll { it }.collect { file(it, checkIfExists: true) }
    dat_ch = channel.value(dat_files)

    // The cache is written to, so it cannot be staged; absolutised so it does not
    // resolve inside the task's work dir. It must be on a shared filesystem.
    cache_dir = params.interpro_cache ? file(params.interpro_cache).toAbsolutePath().toString() : ''
    if (cache_dir && !file(cache_dir).exists()) {
        log.warn "--interpro_cache ${cache_dir} does not exist yet; FETCH_DOMAIN_META will create it."
    }

    // Validated here so a typo fails before any task runs. Must match
    // fetch_domains.py's TIERS.
    def validTiers = ['human_reviewed', 'other_reviewed', 'human_unreviewed', 'other_unreviewed']
    def tiersRaw = params.instance_tiers instanceof List ? params.instance_tiers : params.instance_tiers.toString().split(',')
    def tiers = tiersRaw.collect { it.toString().trim() }.findAll { it }
    if (!tiers) {
        error("--instance_tiers is empty; name at least one of ${validTiers.join(', ')}.")
    }
    def unknownTiers = tiers - validTiers
    if (unknownTiers) {
        error("--instance_tiers has unknown stratum name(s) ${unknownTiers.join(', ')}; valid names are ${validTiers.join(', ')}.")
    }
    log.info "--instance_tiers ${tiers.join(',')}: families without instances in these strata are dropped (see _shared/data/dropped_families.tsv)."
    // A warning only: precomputed datasets never reach FETCH_DOMAIN_META.
    if (tiers.any { it.endsWith('unreviewed') } && !dat_files) {
        log.warn "--instance_tiers names an unreviewed stratum but no --uniprot_dat was given; FETCH_DOMAIN_META will fail without a TrEMBL flat file."
    }

    shared_fetch   = FETCH_DOMAIN_META(
        families_list.map { families -> tuple([id: "_shared"], families) },
        clans_ch,
        regions_ch,
        dat_ch,
        channel.value(cache_dir),
    )
    shared_lengths = GET_LENGTHS_SHARED_DDI(shared_fetch.sequences)

    subset_out = SUBSET_DOMAIN_DATA(
        branched.needs_fetch,
        shared_fetch.sequences.map       { meta, f -> f }.first(),
        shared_fetch.go_annotations.map  { meta, f -> f }.first(),
        shared_fetch.species.map         { meta, f -> f }.first(),
        shared_lengths.map                { meta, f -> f }.first(),
        shared_fetch.instances.map       { meta, f -> f }.first(),
    )

    emit:
    sequences      = precomputed_sequences.mix(subset_out.sequences)
    go_annotations = precomputed_go.mix(subset_out.go_annotations)
    species        = precomputed_species.mix(subset_out.species)
    lengths        = precomputed_lengths.mix(subset_out.lengths)
    instances      = precomputed_instances.mix(subset_out.instances)
    // At most one item, under the "_shared" meta: never join() it per dataset.
    dropped_families = shared_fetch.dropped_families
}
