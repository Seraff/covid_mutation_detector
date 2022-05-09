# COVID-19 Mutation Detector

## cluster_sequences.py

Example:

```
./cluster_sequences.py -f tmp/all_output/gene.S.fasta -o tmp/gene_s_output
```

Notes:

- folder provided with `-o` option will be created if it doesn't exist
- fasta file provided with `-f` option should be aligned protein fasta file from the output of `nextalign` program (it can contain `-` gaps)

## print_nextclade_result.py

Get pretty print of `S` gene for `datafreeze-2022-04-22_12_weeks` data for the sequence with id `EPI_ISL_10132849`:

```
./print_nextclade_result.py -r /auto/vestec1-elixir/projects/cogcz/serafim/df-20220422/datafreeze-2022-04-22_12_weeks/ -g S -s EPI_ISL_10132849
```

Example of output:

```
# EPI_ISL_10132849
ref_na: ATG TTT GTT TTT CTT GTT TTA TTG CCA CTA GTC TCT AGT CAG TGT GTT AAT CTT ACA ACC AGA ACT CAA TTA CCC CCT GCA TAC ACT AAT TCT TTC ACA CGT GGT ...
ref_aa:  M   F   V   F   L   V   L   L   P   L   V   S   S   Q   C   V   N   L   T   T   R   T   Q   L   P   P   A   Y   T   N   S   F   T   R   G  ...
mut_id:                                                                                                                                             ...
seq_na: ATG TTT GTT TTT CTT GTT TTA TTG CCA CTA GTT TCT AGT CAG TGT GTT AAT CTT ACA ACC AGA ACT CAA TTA CCC CCT GCA TAC ACT AAT TCT TTC ACA CGT GGT ...
seq_aa:  M   F   V   F   L   V   L   L   P   L   V   S   S   Q   C   V   N   L   T   T   R   T   Q   L   P   P   A   Y   T   N   S   F   T   R   G  ...
```
