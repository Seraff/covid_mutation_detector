## Example:
## INPUT_PATH="some/path.xz" OUTPUT_PATH="path/to/some/folder" snakemake -c1
## Input can be either xz archive or fasta file

import os
import json
import platform
from pathlib import Path


envvars:
    "INPUT_PATH",
    "OUTPUT_PATH"


INPUT_PATH = os.environ["INPUT_PATH"]
OUTPUT_PATH = os.environ["OUTPUT_PATH"]

MACOS = platform.system() == 'Darwin'
ROOT_PATH = os.getcwd()
COV_DATA_PATH = os.path.join(ROOT_PATH, 'data', 'nextclade_dataset')

INPUT_FASTA_FILENAME = os.path.basename(INPUT_PATH).rsplit('.')[0]
INPUT_FASTA_PATH = os.path.join(OUTPUT_PATH, 'data', f'{INPUT_FASTA_FILENAME}.fasta')


rule all:
    input:
        f'{OUTPUT_PATH}/report.json'
    run:
        readlink_cmd = 'greadlink' if MACOS else 'readlink'
        shell(f'{readlink_cmd} -f {" ".join(input)}')


## Unpacking
rule unxz:
    input:
        INPUT_PATH
    output:
        INPUT_FASTA_PATH
    run:
        if INPUT_PATH.split('.')[-1] == 'xz':
            shell(f"unxz -c {input} > {output}")
        else:
            shell(f"cp {input} {output}")


## Nextclading
rule nextclade:
    input:
        INPUT_FASTA_PATH
    output:
        f'{OUTPUT_PATH}/nextclade.json'
    shell:
        "nextclade run --in-order "
        "--input-fasta {input} --input-dataset {COV_DATA_PATH} "
        "--output-dir {OUTPUT_PATH}/nextclade "
        "--output-json {output} "

## Converting reports
rule report:
    input:
        f'{OUTPUT_PATH}/nextclade.json'
    output:
        f'{OUTPUT_PATH}/report_raw.json'
    run:
        converter = f'{ROOT_PATH}/parse_nextclade_report.py'
        shell(f'{converter} -i {input} -o {output}')


## Finding all suspicious mutations
rule find_suspicious_mutations:
    input:
        f'{OUTPUT_PATH}/report_raw.json'
    output:
        f'{OUTPUT_PATH}/suspicious_all.json'
    run:
        script = f'{ROOT_PATH}/find_suspicious_mutations.py'
        shell(f"{script} -i {input} -d {OUTPUT_PATH}/nextclade -o {output}")


## Finding new suspicious mutations
rule reporting_new_suspicious_mutations:
    input:
        f'{OUTPUT_PATH}/suspicious_all.json'

    output:
        f'{OUTPUT_PATH}/suspicious_new.json'
    run:
        script = f'{ROOT_PATH}/report_new_suspicious_mutations.py'
        shell(f"{script} -s {input} -o {output}")


rule add_suspicious_mutations:
    input:
        f'{OUTPUT_PATH}/suspicious_new.json'
    output:
        f'{OUTPUT_PATH}/report_with_suspicious.json'
    run:
        script = f'{ROOT_PATH}/apply_mut_replacements.py'
        shell(f"{script} -i {OUTPUT_PATH}/report_raw.json -o {output}")

rule add_unknown_sequences:
    input:
        f'{OUTPUT_PATH}/report_with_suspicious.json'
    output:
        f'{OUTPUT_PATH}/report_with_suspicious_with_unknown.json'
    run:
        script = f'{ROOT_PATH}/add_unknown_to_report.py'
        shell(f"{script} -r {input} -n {OUTPUT_PATH}/nextclade -o {output}")

rule finalize_report:
    input:
        f'{OUTPUT_PATH}/report_with_suspicious_with_unknown.json'
    output:
        f'{OUTPUT_PATH}/report.json'
    shell:
        "cp {input} {output}"
