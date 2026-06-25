process FETCH_SEQUENCES {
    publishDir "${params.outdir}", mode: 'copy'
    tag "fetch"

    input:
    path ppis

    output:
    path "sequences.fasta"

    script:
    """
    fetch_sequences.py ${ppis} sequences.fasta
    """
}

process GET_LENGTHS {
    tag "lengths"

    input:
    path fasta

    output:
    path "lengths.tsv"

    script:
    """
    { printf 'protein_id\\tlength\\n'; \
      awk '/^>/{if(acc) print acc"\\t"len; acc=substr(\$1,2); len=0; next} {len+=length(\$0)} END{if(acc) print acc"\\t"len}' ${fasta} \
          | sort; \
    } > lengths.tsv
    """
}

process RUN_BLAST {
    publishDir "${params.outdir}/similarities", mode: 'copy'
    tag "blast"

    input:
    path fasta

    output:
    path "all_vs_all.tsv"

    script:
    """
    makeblastdb -dbtype prot -in ${fasta} -out blastdb -parse_seqids
    blastp \\
        -query ${fasta} \\
        -db blastdb \\
        -outfmt "6 qseqid sseqid evalue bitscore" \\
        -max_hsps 1 \\
        -num_threads ${task.cpus} \\
        -out all_vs_all.tsv
    """
}

process MAKE_METIS {
    publishDir "${params.outdir}/similarities", mode: 'copy'
    tag "metis"

    input:
    path blast_results
    path lengths

    output:
    path "similarity.graph", emit: graph
    path "node_mapping.tsv", emit: node_mapping

    script:
    """
    make_metis.py \\
        ${blast_results} \\
        ${lengths} \\
        similarity.graph \\
        node_mapping.tsv \\
        --edge_weight ${params.edge_weight}
    """
}

process RUN_KAHIP {
    publishDir "${params.outdir}/similarities", mode: 'copy'
    tag "kahip"

    input:
    path graph

    output:
    path "partitioned_proteome.txt"

    script:
    """
    kaffpa \\
        ${graph} \\
        --seed=${params.kahip_seed} \\
        --output_filename=partitioned_proteome.txt \\
        --k=${params.kahip_k} \\
        --preconfiguration=${params.kahip_preconfiguration}
    """
}

process SORT_PPIS {
    tag "sort"

    input:
    path ppis
    path partition
    path fasta
    path node_mapping

    output:
    path "train.csv",          emit: train_ppis
    path "val.csv",            emit: val_ppis
    path "test.csv",           emit: test_ppis
    path "train.fasta",        emit: train_fasta
    path "val.fasta",          emit: val_fasta
    path "test.fasta",         emit: test_fasta
    path "sort_ppis_mqc.json", emit: mqc

    script:
    """
    sort_ppis.py \\
        --ppis         ${ppis} \\
        --partition    ${partition} \\
        --fasta        ${fasta} \\
        --node_mapping ${node_mapping}
    """
}

process CDHIT {
    input:
    path db1_fasta
    path db2_fasta

    output:
    path "cdhit.out"

    script:
    """
    cd-hit-2d \\
        -i  ${db1_fasta} \\
        -i2 ${db2_fasta} \\
        -o  cdhit.out \\
        -c  ${params.cdhit_identity} \\
        -n  ${params.cdhit_wordsize} \\
        -T  ${task.cpus} \\
        -M  4000
    """
}

process REMOVE_REDUNDANT {
    publishDir "${params.outdir}", mode: 'copy'
    tag "remove_redundant"

    input:
    path train_ppis
    path val_ppis
    path test_ppis
    path train_fasta
    path val_fasta
    path test_fasta
    path sim_train_val,  stageAs: 'sim_train_val.out'
    path sim_train_test, stageAs: 'sim_train_test.out'

    output:
    path "train_nr.csv",              emit: train_ppis
    path "val_nr.csv",                emit: val_ppis
    path "test_nr.csv",               emit: test_ppis
    path "train_nr.fasta",            emit: train_fasta
    path "val_nr.fasta",              emit: val_fasta
    path "test_nr.fasta",             emit: test_fasta
    path "remove_redundant_mqc.json", emit: mqc

    script:
    """
    remove_redundant.py \\
        --train_ppis     ${train_ppis} \\
        --val_ppis       ${val_ppis} \\
        --test_ppis      ${test_ppis} \\
        --train_fasta    ${train_fasta} \\
        --val_fasta      ${val_fasta} \\
        --test_fasta     ${test_fasta} \\
        --sim_train_val  ${sim_train_val} \\
        --sim_train_test ${sim_train_test}
    """
}

process SAMPLE_NEGATIVES {
    publishDir "${params.outdir}", mode: 'copy'
    tag "negatives"

    input:
    path train_ppis
    path val_ppis
    path test_ppis

    output:
    path "train_negatives.csv"
    path "val_negatives.csv"
    path "test_negatives.csv"
    path "test_negatives_random.csv"
    path "sample_negatives_mqc.json", emit: mqc

    script:
    """
    sample_negatives.py \\
        --train ${train_ppis} \\
        --val   ${val_ppis} \\
        --test  ${test_ppis}
    """
}

process PREPARE_MQC {
    input:
    path mqc_files

    output:
    path "section_*_mqc.json"

    script:
    """
    python3 -c "
import json, pathlib
sections = []
for f in sorted(pathlib.Path('.').glob('*_mqc.json')):
    sections.extend(json.loads(f.read_text()))
for i, sec in enumerate(sections):
    pathlib.Path(f'section_{i:03d}_mqc.json').write_text(json.dumps(sec, indent=2))
"
    """
}

process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    tag "multiqc"

    input:
    path mqc_files

    output:
    path "multiqc_report.html"
    path "multiqc_report_data"

    script:
    """
    multiqc . --title "PPI Splitting Pipeline" --filename multiqc_report.html
    """
}
