#!/usr/bin/env bash
# Compare two published results trees. Text files are compared order-insensitively
# (sorted), everything else byte-for-byte. MultiQC HTML and Nextflow's own reports
# are skipped -- they carry timestamps.
set -u
A=${1:?usage: cmp_results.sh <tree-a> <tree-b>}
B=${2:?usage: cmp_results.sh <tree-a> <tree-b>}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT

for side in A B; do
    root=$([ $side = A ] && echo "$A" || echo "$B")
    ( cd "$root" && find . -type f ! -path '*/multiqc/*' ! -path './pipeline_info/*' ) | sort > "$T/$side"
done

comm -23 "$T/A" "$T/B" | sed 's|^\./|ONLY IN A     |'
comm -13 "$T/A" "$T/B" | sed 's|^\./|ONLY IN B     |'

# One line per FASTA record (">id<TAB>sequence"), so sorting keeps each header
# with its own sequence -- sorting raw lines would compare a re-paired file equal.
records() {
    awk '/^>/ { if (h != "") print h "\t" s; h = $0; s = ""; next } { s = s $0 } END { if (h != "") print h "\t" s }' "$1"
}

same=0; differ=0
while IFS= read -r f; do
    case "$f" in
        *.fasta)
            if [ "$(records "$A/$f" | sort | md5sum)" = "$(records "$B/$f" | sort | md5sum)" ]
                then same=$((same+1)); else printf 'DIFF          %s\n' "${f#./}"; differ=$((differ+1)); fi ;;
        *.csv|*.tsv|*.txt|*.graph|*.json|*.yml)
            if [ "$(sort -- "$A/$f" | md5sum)" = "$(sort -- "$B/$f" | md5sum)" ]
                then same=$((same+1)); else printf 'DIFF          %s\n' "${f#./}"; differ=$((differ+1)); fi ;;
        *)
            if cmp -s -- "$A/$f" "$B/$f"
                then same=$((same+1)); else printf 'DIFF (bytes)  %s\n' "${f#./}"; differ=$((differ+1)); fi ;;
    esac
done < <(comm -12 "$T/A" "$T/B")

printf '\n%d file(s) identical, %d differing, %d only in A, %d only in B\n' \
    "$same" "$differ" "$(comm -23 "$T/A" "$T/B" | wc -l)" "$(comm -13 "$T/A" "$T/B" | wc -l)"
[ "$differ" -eq 0 ] && [ "$(comm -3 "$T/A" "$T/B" | wc -l)" -eq 0 ]
