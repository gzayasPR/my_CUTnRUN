#!/bin/bash
IFS=$'\n\t'

usage() {
    echo "Usage: $0 <RESULT_DIR> <ECOLI_OUTPUT_TSV> <MOUSE_OUTPUT_TSV> <TOTAL_OUTPUT_TSV>"
    exit 1
}

# --- PARSE ARGS ---
[[ $# -eq 4 ]] || usage
RESULT_DIR=$1
ECOLI_OUTPUT_TSV=$2
MOUSE_OUTPUT_TSV=$3
TOTAL_OUTPUT_TSV=$4


# --- PREP OUTPUT DIRS ---
mkdir -p "$(dirname "$ECOLI_OUTPUT_TSV")"
mkdir -p "$(dirname "$MOUSE_OUTPUT_TSV")"
mkdir -p "$(dirname "$TOTAL_OUTPUT_TSV")"

# --- FUNCTION TO AGGREGATE ---
# Args: suffix (e.g. "ecoli_read_count.txt"), out_tsv, header_label
aggregate_counts() {
    local suffix=$1
    local out_tsv=$2
    local header=$3

    printf "Sample\t%s\n" "$header" > "$out_tsv"

    find "$RESULT_DIR" -type f -name "*_${suffix}" -print0 \
    | while IFS= read -r -d '' file; do
        # strip off _suffix from basename
        sample=${file##*/}
        sample=${sample%"_$suffix"}
        # read the count
        count=$(<"$file")
        printf "%s\t%s\n" "$sample" "$count"
    done >> "$out_tsv"
}

# --- RUN IT ---
aggregate_counts "ecoli_read_count.txt" "$ECOLI_OUTPUT_TSV" "Ecoli_Reads"
aggregate_counts "mouse_read_count.txt" "$MOUSE_OUTPUT_TSV" "Mouse_Reads"
aggregate_counts "total_read_count.txt" "$TOTAL_OUTPUT_TSV" "Total_Reads"