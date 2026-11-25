process FETCH_DATA {
    publishDir(path: { "${params.outdir}/${meta.id}/data" }, mode: 'copy')
    tag "${meta.id}"

    input:
    tuple val(meta), path(proteins)  // plain text file, one protein ID per line

    output:
    tuple val(meta), path("sequences.fasta"),    emit: sequences
    tuple val(meta), path("go_annotations.tsv"), emit: go_annotations
    tuple val(meta), path("species.tsv"),        emit: species

    script:
    """
    fetch_data.py ${proteins} sequences.fasta go_annotations.tsv species.tsv \\
        --exclude-evidence "${params.go_exclude_evidence}"
    """
}

// DDI mode's FETCH_DATA: turns Pfam family accessions into sampled domain
// instances plus the same sequences/species/GO files, from bulk Pfam/UniProt data.
process FETCH_DOMAIN_META {
    publishDir(path: { "${params.outdir}/${meta.id}/data" }, mode: 'copy')
    tag "${meta.id}"
    // A mid-stream truncation of the Pfam download kills the task; retry it.
    label 'error_retry'

    input:
    tuple val(meta), path(families)  // plain text, one Pfam family accession per line
    path clans                       // Pfam-A.clans.tsv(.gz), or [] to download it
    path pfam_regions                // local Pfam-A.regions.tsv.gz, or [] to stream it from EBI
    path uniprot_dats                // UniProt flat file(s), Swiss-Prot first, or [] to download Swiss-Prot
    val cache_dir                    // absolute cache directory, or '' for no cache

    output:
    tuple val(meta), path("sequences.fasta"),        emit: sequences
    tuple val(meta), path("go_annotations.tsv"),     emit: go_annotations
    tuple val(meta), path("species.tsv"),            emit: species
    tuple val(meta), path("instances.tsv"),          emit: instances
    // Always written, header-only when nothing was dropped.
    tuple val(meta), path("dropped_families.tsv"),   emit: dropped_families

    script:
    def pool_size   = (params.ddi_examples_target as int) * (params.ddi_examples_pool_factor as int)
    def clans_arg   = clans      ? "--clans ${clans}"           : ''
    def regions_arg = pfam_regions ? "--pfam-regions ${pfam_regions}" : ''
    // Order matters (first file wins on a duplicate accession). A single staged
    // path arrives as a Path, not a list.
    def dats        = uniprot_dats instanceof List ? uniprot_dats : (uniprot_dats ? [uniprot_dats] : [])
    def dat_arg     = dats.collect { "--uniprot-dat ${it}" }.join(' ')
    // params.instance_tiers may be a List (config) or a comma-separated String (CLI).
    def tiers       = params.instance_tiers instanceof List ? params.instance_tiers.join(',') : params.instance_tiers
    def cache_arg   = cache_dir  ? "--interpro-cache ${cache_dir}" : ''
    def release_arg = params.pfam_release ? "--pfam-release ${params.pfam_release}" : ''
    """
    fetch_domains.py \\
        --families  ${families} \\
        --pool-size ${pool_size} \\
        --seed      ${params.seed} \\
        --tiers     ${tiers} \\
        ${clans_arg} \\
        ${cache_arg} \\
        ${regions_arg} \\
        ${dat_arg} \\
        ${release_arg}
    """
}

process GET_LENGTHS {
    tag "${meta.id}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("lengths.tsv")

    script:
    """
    { printf 'protein_id\\tlength\\n'; \
      awk '/^>/{if(acc) print acc"\\t"len; acc=substr(\$1,2); len=0; next} {len+=length(\$0)} END{if(acc) print acc"\\t"len}' ${fasta} \
          | sort; \
    } > lengths.tsv
    """
}

// Splits the shared FETCH_DATA/GET_LENGTHS batch back out per dataset, so
// downstream steps (e.g. BLAST) see only this dataset's own proteins.
process SUBSET_FETCHED_DATA {
    publishDir(path: { "${params.outdir}/${meta.id}/data" }, mode: 'copy', saveAs: { f -> f == 'lengths.tsv' ? null : f })
    tag "${meta.id}"

    input:
    tuple val(meta), path(ppis)
    path shared_sequences,      stageAs: 'shared_sequences.fasta'
    path shared_go_annotations, stageAs: 'shared_go_annotations.tsv'
    path shared_species,        stageAs: 'shared_species.tsv'
    path shared_lengths,        stageAs: 'shared_lengths.tsv'

    output:
    tuple val(meta), path("sequences.fasta"),    emit: sequences
    tuple val(meta), path("go_annotations.tsv"), emit: go_annotations
    tuple val(meta), path("species.tsv"),        emit: species
    tuple val(meta), path("lengths.tsv"),        emit: lengths

    script:
    """
    subset_fetched_data.py \\
        --ppis ${ppis} \\
        --sequences ${shared_sequences} \\
        --go_annotations ${shared_go_annotations} \\
        --species ${shared_species} \\
        --lengths ${shared_lengths} \\
        --out_sequences sequences.fasta \\
        --out_go_annotations go_annotations.tsv \\
        --out_species species.tsv \\
        --out_lengths lengths.tsv
    """
}

// DDI mode's SUBSET_FETCHED_DATA: the CSV holds families but sequences are
// instance-keyed, so the keep-set is resolved through instances.tsv.
process SUBSET_DOMAIN_DATA {
    publishDir(path: { "${params.outdir}/${meta.id}/data" }, mode: 'copy', saveAs: { f -> f == 'lengths.tsv' ? null : f })
    tag "${meta.id}"

    input:
    tuple val(meta), path(ddis)
    path shared_sequences,      stageAs: 'shared_sequences.fasta'
    path shared_go_annotations, stageAs: 'shared_go_annotations.tsv'
    path shared_species,        stageAs: 'shared_species.tsv'
    path shared_lengths,        stageAs: 'shared_lengths.tsv'
    path shared_instances,      stageAs: 'shared_instances.tsv'

    output:
    tuple val(meta), path("sequences.fasta"),    emit: sequences
    tuple val(meta), path("go_annotations.tsv"), emit: go_annotations
    tuple val(meta), path("species.tsv"),        emit: species
    tuple val(meta), path("lengths.tsv"),        emit: lengths
    tuple val(meta), path("instances.tsv"),      emit: instances

    script:
    """
    subset_fetched_data.py \\
        --ppis ${ddis} \\
        --sequences ${shared_sequences} \\
        --go_annotations ${shared_go_annotations} \\
        --species ${shared_species} \\
        --lengths ${shared_lengths} \\
        --instances ${shared_instances} \\
        --out_sequences sequences.fasta \\
        --out_go_annotations go_annotations.tsv \\
        --out_species species.tsv \\
        --out_lengths lengths.tsv \\
        --out_instances instances.tsv
    """
}

// Header-only GO table for DDI datasets with precomputed domain data, since
// --go_annotations stays a required flag downstream.
process EMPTY_GO_ANNOTATIONS {
    publishDir(path: { "${params.outdir}/${meta.id}/data" }, mode: 'copy')
    tag "${meta.id}"

    input:
    val meta

    output:
    tuple val(meta), path("go_annotations.tsv")

    script:
    """
    printf 'protein_id\\tgo_bp\\tgo_mf\\tgo_cc\\n' > go_annotations.tsv
    """
}
