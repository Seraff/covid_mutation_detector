#!/bin/bash
set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
INPUT_PATH="${SCRIPT_DIR}/data/seq.fasta"
OUTPUT_PATH="${SCRIPT_DIR}/data/output"
rm -rf "${OUTPUT_PATH}/*"

$SCRIPT_DIR/../detect_covid_mutations.py -i $INPUT_PATH -o $OUTPUT_PATH
