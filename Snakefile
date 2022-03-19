## Example:
## DF_SHORT="df-20220114" DF_LONG="datafreeze-2022-01-14_12_weeks" snakemake -c1

import json

DIR_PREFIX = "/storage/vestec1-elixir/projects/cogcz"
COV_DATA_PATH = f"{DIR_PREFIX}/serafim/sars-cov-2"
DETECTOR_PATH = f"{DIR_PREFIX}/serafim/covid_mutation_detector"

envvars:
    "DF_SHORT",
    "DF_LONG"

def nextclade_all_output_path():
    return f'{DIR_PREFIX}/serafim/{os.environ["DF_SHORT"]}/{os.environ["DF_LONG"]}'

def nextclade_good_output_path():
    return f'{nextclade_all_output_path()}_good'

rule all:
    input:
        f'{DIR_PREFIX}/serafim/{os.environ["DF_SHORT"]}/{os.environ["DF_LONG"]}_good/report_fixed.json',
        f'{DIR_PREFIX}/serafim/{os.environ["DF_SHORT"]}/{os.environ["DF_LONG"]}/report_fixed.json'
    run:
        shell(f'readlink -f {" ".join(input)}')

## Unpacking
rule unxz:
    input:
        f'{DIR_PREFIX}/seq/{os.environ["DF_SHORT"]}/{os.environ["DF_LONG"]}_good.fasta.xz',
        f'{DIR_PREFIX}/seq/{os.environ["DF_SHORT"]}/{os.environ["DF_LONG"]}.fasta.xz'
    output:
        f'{DIR_PREFIX}/serafim/{os.environ["DF_SHORT"]}/data/{os.environ["DF_LONG"]}_12_weeks_good.fasta',
        f'{DIR_PREFIX}/serafim/{os.environ["DF_SHORT"]}/data/{os.environ["DF_LONG"]}_12_weeks.fasta'
    run:
        for i in range(2):
            cmd = f"unxz -c {input[i]} > {output[i]}"
            shell(cmd)

## Nextclading
rule nextclade:
    input:
        f'{DIR_PREFIX}/serafim/{os.environ["DF_SHORT"]}/data/{os.environ["DF_LONG"]}_12_weeks_good.fasta',
        f'{DIR_PREFIX}/serafim/{os.environ["DF_SHORT"]}/data/{os.environ["DF_LONG"]}_12_weeks.fasta'
    output:
        f'{ nextclade_good_output_path() }/nextclade.json',
        f'{ nextclade_all_output_path() }/nextclade.json'
    run:
        base_folder_name = f'{DIR_PREFIX}/serafim/{os.environ["DF_SHORT"]}/{os.environ["DF_LONG"]}'
        folders = [base_folder_name + '_good', base_folder_name]

        for i, folder in enumerate(folders):
            shell(f"echo working with file {i+1}...")
            cmd = f"nextclade --in-order \
                    --input-fasta {input[i]} --input-dataset {COV_DATA_PATH} \
                    --output-dir {folder} \
                    --output-json {output[i]}"
            shell(cmd)

## Converting reports
rule report:
    input:
        f'{ nextclade_good_output_path() }/nextclade.json',
        f'{ nextclade_all_output_path() }/nextclade.json'
    output:
        f'{ nextclade_good_output_path() }/report.json',
        f'{ nextclade_all_output_path() }/report.json'
    run:
        converter = f'{DETECTOR_PATH}/parse_nextclade_report.py'

        for i in range(2):
            shell(f'{converter} -i {input[i]} -o {output[i]}')

## Finding suspicious mutations
rule find_suspicious_mutations:
    input:
        f'{ nextclade_good_output_path() }/report.json',
        f'{ nextclade_all_output_path() }/report.json'

    output:
        f'{ nextclade_good_output_path() }/all_suspicious.json',
        f'{ nextclade_all_output_path() }/all_suspicious.json'
    run:
        script = f'{DETECTOR_PATH}/find_suspicious_mutations.py'
        dirs = [nextclade_good_output_path(), nextclade_all_output_path()]

        for i in range(2):
            shell(f"{script} -i {input[i]} -d {dirs[i]} -o {output[i]}")

## Finding suspicious mutations
rule reporting_new_suspicious_mutations:
    input:
        f'{ nextclade_good_output_path() }/all_suspicious.json',
        f'{ nextclade_all_output_path() }/all_suspicious.json'

    output:
        f'{ nextclade_good_output_path() }/suspicious_new.json',
        f'{ nextclade_all_output_path() }/suspicious_new.json'
    run:
        script = f'{DETECTOR_PATH}/report_new_suspicious_mutations.py'

        for i in range(2):
            shell(f"{script} -s {input[i]} -o {output[i]}")

rule fix_report:
    input:
        f'{ nextclade_good_output_path() }/suspicious_new.json',
        f'{ nextclade_all_output_path() }/suspicious_new.json'
    output:
        f'{ nextclade_good_output_path() }/report_fixed.json',
        f'{ nextclade_all_output_path() }/report_fixed.json'
    run:
        script = f'{DETECTOR_PATH}/apply_mut_replacements.py'
        dirs = [nextclade_good_output_path(), nextclade_all_output_path()]
        reports = [os.path.join(e, 'report.json') for e in dirs]

        for i in range(2):
            shell(f"{script} -i {reports[i]} -o {output[i]}")
