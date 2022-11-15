# COVID-19 Mutation Detector

## Main Pipeline

### Generating Per-mutation Statistics

The main pipeline takes **FASTA file** with multiple sequences and generates **JSON file** with mutation statistics in the following format:

```
{
    "E:N66K": {
        "cnt": 2,
        "seq_ids": [
            "EPI_ISL_4882657",
            "EPI_ISL_4882659"
        ]
    },
    ...
}
```
#### Preparation

1. Install [Conda](https://docs.conda.io/en/latest/)

2. Install [Mamba](https://github.com/mamba-org/mamba)

```bash
conda install -n base -c conda-forge mamba
```

#### How To Run

##### 1. Navigate to the repository folder

```bash
cd covid_mutation_detector
```

##### 2. Run the pipeline in the following format

```bash
INPUT_PATH="<path_to_fasta_file>" OUTPUT_PATH="<output_folder_path>" snakemake -c1 --use-conda
```

Example:

```bash
INPUT_PATH="/home/alice/covid/samples.fasta" OUTPUT_PATH="/home/alice/covid/analysis" snakemake -c1 --use-conda
```

As input you can also provide an archived `xz` fasta file (i.e. `/home/alice/covid/samples.fasta.xz`).

##### 3. Fix Suspisious Mutations (if they appear)

Sometimes mutations are predicted in a wrong way, not the way, it's presented on https://outbreak.info/.

The problem occures in replacement mutation after the chain of deletions.

Example:

```
Positions:                          156  157  158  159  160
Reference AA:    S    W    M    E    S    E    F    R    V    Y    S    S    A
Reference NA:   AGT  TGG  ATG  GAA  AGT  GAG  TTC  AGA  GTT  TAT  TCT  AGT  GCG
Mutations:                                *    *    *
Sequence AA      X    W    M    E    S    -    -    G    V    Y    S    S    A
Sequence NA     NGT  TGG  ATG  GAA  AGT  G--  ---  -GA  GTT  TAT  TCT  AGT  GCG
```
Our approach interpretes this situation as `S:E156-`, `S:F157-`, `S:R158G`.

Nextclade puts the replacement instead of first deletion: `S:E156G`, `S:F157-`, `S:R158-`.

When such situation happen, the pipeline stops with an error. The output folder will contain file `suspicious_new.json` will have a detailed list of all the suspicious situations.

The format of a suspicious situation is presented in the follwing format:

```
[
    {
        "predict": "S:E156-,S:F157-,S:R158G",
        "verdict": "",
        "ref_aa": " S   W   M   E   S   E   F   R   V   Y   S   S   A ",
        "ref_na": "AGT TGG ATG GAA AGT GAG TTC AGA GTT TAT TCT AGT GCG",
        "mut_id": "                     *   *   *                     ",
        "seq_aa": " X   W   M   E   S   -   -   G   V   Y   S   S   A ",
        "seq_na": "NGT TGG ATG GAA AGT G-- --- -GA GTT TAT TCT AGT GCG",
        "notes": "",
        "cnt": 2,
        "seq_ids": ["EPI_ISL_4882657", "EPI_ISL_4882659"]
    }
]
```
* `predict` - predicted mutations
* `verdict` - how to interpret it (in the `predict` format)
* `ref_aa` - reference amino acid sequence
* `ref_na` - reference nucleotide sequence
* `mut_id` - highlight the places of mutation
* `seq_aa` - our amino acid sequence
* `seq_na` - our nucleotide sequence
* `notes` - any notes we can add about this situation
* `cnt` - number of sequences with suspicious situation in the dataset
* `seq_ids` - headers of sequences containing suspicious situation

Our software already contains the list of the most common suspicious situation in the file `mut_replacement_rules.json`.

However, if new suspicious situation is found, user should resolve it.

* Open `suspicious_new.json`
* Analyze all the situations and fill `verdict` fields
  * if no changes are needed just copy `predict` value to `verdict` field
* The fixed file should look like this:
```
[
    {
        "predict": "S:E156-,S:F157-,S:R158G",
        "verdict": "S:E156G,S:F157-,S:R158-",
        "ref_aa": " S   W   M   E   S   E   F   R   V   Y   S   S   A ",
        "ref_na": "AGT TGG ATG GAA AGT GAG TTC AGA GTT TAT TCT AGT GCG",
        "mut_id": "                     *   *   *                     ",
        "seq_aa": " X   W   M   E   S   -   -   G   V   Y   S   S   A ",
        "seq_na": "NGT TGG ATG GAA AGT G-- --- -GA GTT TAT TCT AGT GCG",
        "notes": "S:E156G has 3'856'616 entries in outbreak.info, while S:R158G has only 5'211",
        "cnt": 2,
        "seq_ids": ["EPI_ISL_4882657", "EPI_ISL_4882659"]
    }
]
```
* Copy entries from edited `suspicious_new.json` and append them to `mut_replacement_rules.json` file in the repository of the project
  * `mut_replacement_rules.json` contains a list of objects ( fields surrounded by `{` and `}`), the new objects should be appended to the end of the list
  * The new `mut_replacement_rules.json` should look like:
```
[
    ...
    { ... already known rules ... },
    {
        "predict": "S:E156-,S:F157-,S:R158G",
        "verdict": "S:E156G,S:F157-,S:R158-",
        "ref_aa": " S   W   M   E   S   E   F   R   V   Y   S   S   A ",
        "ref_na": "AGT TGG ATG GAA AGT GAG TTC AGA GTT TAT TCT AGT GCG",
        "mut_id": "                     *   *   *                     ",
        "seq_aa": " X   W   M   E   S   -   -   G   V   Y   S   S   A ",
        "seq_na": "NGT TGG ATG GAA AGT G-- --- -GA GTT TAT TCT AGT GCG",
        "notes": "S:E156G has 3'856'616 entries in outbreak.info, while S:R158G has only 5'211",
        "cnt": 2,
        "seq_ids": ["EPI_ISL_4882657", "EPI_ISL_4882659"]
    }
]
```
* Re-run the pipeline. If all the problems are resolved, there should be no error displayed.

##### 4. Obtaining the results

The output folder will contain the following files
```
.
├── data                                - unarchived sequences (same as input)
│   └── sequences.fasta
├── nextclade                           - Nextclade files
│   ├── sequences.aligned.fasta
│   ├── sequences.errors.csv
│   ├── sequences.gene.E.fasta
│   ├── sequences.gene.M.fasta
│   ├── sequences.gene.N.fasta
│   ├── sequences.gene.ORF1a.fasta
│   ├── sequences.gene.ORF1b.fasta
│   ├── sequences.gene.ORF3a.fasta
│   ├── sequences.gene.ORF6.fasta
│   ├── sequences.gene.ORF7a.fasta
│   ├── sequences.gene.ORF7b.fasta
│   ├── sequences.gene.ORF8.fasta
│   ├── sequences.gene.ORF9b.fasta
│   ├── sequences.gene.S.fasta
│   └── sequences.insertions.csv
├── nextclade.json                      - Nextclade report
├── report.json                         - Report without suspicious situations fixed
├── report_fixed.json                   - Report with suspicious situations fixed
├── suspicious_all.json                 - All suspicious situations found in the dataset
└── suspicious_new.json                 - New suspicious situations that should be fixed (contains empty list if everything is correct)
```

## Additional Utilities

### cluster_sequences.py

Example:

```
scripts/cluster_sequences.py -f tmp/all_output/gene.S.fasta -o tmp/gene_s_output
```

Notes:

- folder provided with `-o` option will be created if it doesn't exist
- fasta file provided with `-f` option should be aligned protein fasta file from the output of `nextalign` program (it can contain `-` gaps)

## print_nextclade_result.py

Get pretty print of `S` gene for `datafreeze-2022-04-22_12_weeks` data for the sequence with id `EPI_ISL_10132849`:

```
scripts/print_nextclade_result.py -r /auto/vestec1-elixir/projects/cogcz/serafim/df-20220422/datafreeze-2022-04-22_12_weeks/ -g S -s EPI_ISL_10132849
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
