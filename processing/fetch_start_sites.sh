#!/usr/bin/env bash
# usage: fetch_start_sites.sh /path/to/samtools input.bam output.bed
# dependence: awk

set -euo pipefail

samtools_bin="$1"
input_bam="$2"
output_bed="$3"

"$samtools_bin" view "$input_bam" | \
awk 'BEGIN{OFS="\t"}{
    read_len = length($10)
    if (and($2,16)) {
        strand = "-"
        end_pos = $4 + read_len - 1
        start = end_pos - 1
        end   = end_pos
    } else {
        strand = "+"
        start = $4 - 1
        end   = $4
    }
    print $3, start, end, strand
}' | \
sort -k1,1 -k2,2n -k4,4 | \
uniq -c | \
awk 'BEGIN{OFS="\t"}{
    print $2, $3, $4, $1, $5
}' > "$output_bed"
