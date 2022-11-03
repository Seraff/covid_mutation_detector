#!/bin/bash

set -e

PIPELINE_OUTPUT_PATH=$1
OUTPUT_PATH=${PIPELINE_OUTPUT_PATH}/final_analysis

mkdir -p ${OUTPUT_PATH}


## Get new dataset

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd  )
ROOT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )/../.." &> /dev/null && pwd  )
DATASET_DIR="${ROOT_DIR}/data/nextclade_dataset"

nextclade dataset get --name 'sars-cov-2' --output-dir ${DATASET_DIR}

## Init Nextclade & Pangolin paths

NEXTCLADE_TSV=${PIPELINE_OUTPUT_PATH}/nextclade/nextclade.tsv
PANGOLIN_CSV=${PIPELINE_OUTPUT_PATH}/pangolin/lineage_report.csv

##-- compare nextclade-pango and UShER (columns are isolate, Nextclade-pango, UShER)
diff -y <(cut -f1,3 --output-delimiter="," ${NEXTCLADE_TSV} | sort -r) \
        <(cut -f1,2 -d"," ${PANGOLIN_CSV} | sort -r) | \
	sed -r 's/[[:space:]]+/,/g' | sed 's/,|,/,/g' | awk -F"," '{if (($1 == $3) || ($1 == "seqName")) {print $1 "\t" $2 "\t" $4}}' | \
	sed 's/seqName\tNextclade_pango\tlineage/Isolate\tNextclade_pango\tUShER/g' > ${OUTPUT_PATH}/lineageAssignmentChanges.tsv


##-- compare UShER and Nextclade
diff -y <(cut -f1,2 -d"," ${PANGOLIN_CSV} | sort -r) \
        <(cut -f1,2 ${NEXTCLADE_TSV} | sort -r) | \
	sed -r 's/[[:space:]]+/,/g' | sed 's/,|,/,/g' | \
	awk -F"," '{if (($1 == $3) || ($1 == "taxon")) {print $1 "\t" $2 "\t" $4 "," $5}}' | \
	sed 's/taxon\tlineage\tclade,/Isolate\tUShER\tNextclade/g' | sed 's/"//g' > ${OUTPUT_PATH}/cladeAssignment.tsv


##-- collect sequences for recombinant reanalysis
## unassigned by pango
grep -i unassigned ${PANGOLIN_CSV} | cut -d"," -f1 | sort -u >  ${OUTPUT_PATH}/putativeRecombinants.tmp
## recombinant in nextclade
grep -i recombinant ${NEXTCLADE_TSV} | cut -f1 | sort -u >> ${OUTPUT_PATH}/putativeRecombinants.tmp

sort -u ${OUTPUT_PATH}/putativeRecombinants.tmp > ${OUTPUT_PATH}/putativeRecombinants.txt
rm ${OUTPUT_PATH}/putativeRecombinants.tmp

NEXTCLADE_ALN_FASTA=${PIPELINE_OUTPUT_PATH}/nextclade/nextclade.aligned.fasta
cat ${OUTPUT_PATH}/putativeRecombinants.txt | xargs -I % sh -c "grep -A1 % ${NEXTCLADE_ALN_FASTA}" > ${OUTPUT_PATH}/putativeRecombinants.fasta


# Drawing plots
Rscript ${SCRIPT_DIR}/lineageAssignment.R ${PIPELINE_OUTPUT_PATH}
Rscript ${SCRIPT_DIR}/cladeAssignment.R ${PIPELINE_OUTPUT_PATH}
