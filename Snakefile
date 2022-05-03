## Example:
## DF_SHORT="df-20220114" DF_LONG="datafreeze-2022-01-14" snakemake -c1

import os
import json

DIR_PREFIX = "/storage/vestec1-elixir/projects/cogcz"
COV_DATA_PATH = f"{DIR_PREFIX}/serafim/sars-cov-2"
DETECTOR_PATH = f"{DIR_PREFIX}/serafim/covid_mutation_detector"


envvars:
    "DF_SHORT",
    "DF_LONG"


def nextclade_12_weeks_all_output_path():
    return f'{DIR_PREFIX}/serafim/{os.environ["DF_SHORT"]}/{os.environ["DF_LONG"]}_12_weeks'


def nextclade_12_weeks_good_output_path():
    return f'{nextclade_12_weeks_all_output_path()}_good'


def nextclade_only_current_output_path():
    return f'{DIR_PREFIX}/serafim/{os.environ["DF_SHORT"]}/{os.environ["DF_LONG"]}'


rule all:
    input:
        f'{DIR_PREFIX}/serafim/{os.environ["DF_SHORT"]}/{os.environ["DF_LONG"]}_12_weeks_good/report_fixed.json',
        f'{DIR_PREFIX}/serafim/{os.environ["DF_SHORT"]}/{os.environ["DF_LONG"]}_12_weeks/report_fixed.json',
        f'{DIR_PREFIX}/serafim/{os.environ["DF_SHORT"]}/{os.environ["DF_LONG"]}/report_fixed.json'
    run:
        shell(f'readlink -f {" ".join(input)}')

## Unpacking
rule unxz:
    input:
        f'{DIR_PREFIX}/seq/{os.environ["DF_SHORT"]}/{os.environ["DF_LONG"]}_12_weeks_good.fasta.xz',
        f'{DIR_PREFIX}/seq/{os.environ["DF_SHORT"]}/{os.environ["DF_LONG"]}_12_weeks.fasta.xz',
        f'{DIR_PREFIX}/seq/{os.environ["DF_SHORT"]}/{os.environ["DF_LONG"]}.fasta.xz'
    output:
        f'{DIR_PREFIX}/serafim/{os.environ["DF_SHORT"]}/data/{os.environ["DF_LONG"]}_12_weeks_good.fasta',
        f'{DIR_PREFIX}/serafim/{os.environ["DF_SHORT"]}/data/{os.environ["DF_LONG"]}_12_weeks.fasta',
        f'{DIR_PREFIX}/serafim/{os.environ["DF_SHORT"]}/data/{os.environ["DF_LONG"]}.fasta'
    run:
        for i in range(3):
            cmd = f"unxz -c {input[i]} > {output[i]}"
            shell(cmd)

## Nextclading
rule nextclade:
    input:
        f'{DIR_PREFIX}/serafim/{os.environ["DF_SHORT"]}/data/{os.environ["DF_LONG"]}_12_weeks_good.fasta',
        f'{DIR_PREFIX}/serafim/{os.environ["DF_SHORT"]}/data/{os.environ["DF_LONG"]}_12_weeks.fasta',
        f'{DIR_PREFIX}/serafim/{os.environ["DF_SHORT"]}/data/{os.environ["DF_LONG"]}.fasta'
    output:
        f'{ nextclade_12_weeks_good_output_path() }/nextclade.json',
        f'{ nextclade_12_weeks_all_output_path() }/nextclade.json',
        f'{ nextclade_only_current_output_path() }/nextclade.json'
    run:
        for i, out_path in enumerate(output):
            shell(f"echo working with file {i+1}...")
            cmd = f"nextclade --in-order \
                    --input-fasta {input[i]} --input-dataset {COV_DATA_PATH} \
                    --output-dir {os.path.dirname(out_path)} \
                    --output-json {output[i]}"
            shell(cmd)

## Converting reports
rule report:
    input:
        f'{ nextclade_12_weeks_good_output_path() }/nextclade.json',
        f'{ nextclade_12_weeks_all_output_path() }/nextclade.json',
        f'{ nextclade_only_current_output_path() }/nextclade.json'
    output:
        f'{ nextclade_12_weeks_good_output_path() }/report.json',
        f'{ nextclade_12_weeks_all_output_path() }/report.json',
        f'{ nextclade_only_current_output_path() }/report.json'
    run:
        converter = f'{DETECTOR_PATH}/parse_nextclade_report.py'

        for i in range(3):
            shell(f'{converter} -i {input[i]} -o {output[i]}')

## Finding suspicious mutations
rule find_suspicious_mutations:
    input:
        f'{ nextclade_12_weeks_good_output_path() }/report.json',
        f'{ nextclade_12_weeks_all_output_path() }/report.json',
        f'{ nextclade_only_current_output_path() }/report.json'

    output:
        f'{ nextclade_12_weeks_good_output_path() }/all_suspicious.json',
        f'{ nextclade_12_weeks_all_output_path() }/all_suspicious.json',
        f'{ nextclade_only_current_output_path() }/all_suspicious.json'
    run:
        script = f'{DETECTOR_PATH}/find_suspicious_mutations.py'

        for i, input_path in enumerate(input):
            shell(f"{script} -i {input[i]} -d {os.path.dirname(input_path)} -o {output[i]}")

## Finding suspicious mutations
rule reporting_new_suspicious_mutations:
    input:
        f'{ nextclade_12_weeks_good_output_path() }/all_suspicious.json',
        f'{ nextclade_12_weeks_all_output_path() }/all_suspicious.json',
        f'{ nextclade_only_current_output_path() }/all_suspicious.json'

    output:
        f'{ nextclade_12_weeks_good_output_path() }/suspicious_new.json',
        f'{ nextclade_12_weeks_all_output_path() }/suspicious_new.json',
        f'{ nextclade_only_current_output_path() }/suspicious_new.json'
    run:
        script = f'{DETECTOR_PATH}/report_new_suspicious_mutations.py'

        for i in range(3):
            shell(f"{script} -s {input[i]} -o {output[i]}")

rule fix_report:
    input:
        f'{ nextclade_12_weeks_good_output_path() }/suspicious_new.json',
        f'{ nextclade_12_weeks_all_output_path() }/suspicious_new.json',
        f'{ nextclade_only_current_output_path() }/suspicious_new.json'
    output:
        f'{ nextclade_12_weeks_good_output_path() }/report_fixed.json',
        f'{ nextclade_12_weeks_all_output_path() }/report_fixed.json',
        f'{ nextclade_only_current_output_path() }/report_fixed.json'
    run:
        script = f'{DETECTOR_PATH}/apply_mut_replacements.py'
        reports = [os.path.join(os.path.dirname(e), 'report.json') for e in input]

        for i in range(3):
            shell(f"{script} -i {reports[i]} -o {output[i]}")
