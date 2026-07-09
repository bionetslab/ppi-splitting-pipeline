process FETCH_DATA {
    publishDir "${params.outdir}/data", mode: 'copy'
    tag "fetch"

    input:
    path ppis

    output:
    path "sequences.fasta",    emit: sequences
    path "go_annotations.tsv", emit: go_annotations
    path "species.tsv",        emit: species

    script:
    """
    fetch_data.py ${ppis} sequences.fasta go_annotations.tsv species.tsv
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
