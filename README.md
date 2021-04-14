# COVID-19 Mutation Detector

## cluster_sequences.py

Example:

`./cluster_sequences.py -f tmp/all_output/gene.S.fasta -o tmp/gene_s_output`

Notes:

- folder provided with `-o` option will be created if it doesn't exist
- fasta file provided with `-f` option should be aligned protein fasta file from the output of `nextalign` program (it can contain `-` gaps)
