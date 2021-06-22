#!/bin/bash
set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
INPUT_PATH="${SCRIPT_DIR}/data/seq.fasta"
OUTPUT_PATH="${SCRIPT_DIR}/data/output"
ORIG_OUTPUT_PATH="${SCRIPT_DIR}/data/original_output"
rm -rf "${OUTPUT_PATH}/*"

$SCRIPT_DIR/../detect_covid_mutations.py -i $INPUT_PATH -o $OUTPUT_PATH

cmp --silent $OUTPUT_PATH/report.json $ORIG_OUTPUT_PATH/report.json && \
    echo 'report.json files are similar' || \
    echo 'report.json files are NOT similar'
