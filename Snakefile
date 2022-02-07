## Example:
## DF_SHORT="df-20220114" DF_LONG="datafreeze-2022-01-14_12_weeks" snakemake -c1

DIR_PREFIX = "/storage/vestec1-elixir/projects/cogcz"
COV_DATA_PATH = f"{DIR_PREFIX}/serafim/sars-cov-2"
DETECTOR_PATH = f"{DIR_PREFIX}/serafim/covid_mutation_detector"

envvars:
    "DF_SHORT",
    "DF_LONG"

rule all:
    input:
        f'{DIR_PREFIX}/serafim/{os.environ["DF_SHORT"]}/{os.environ["DF_LONG"]}_good/report.json',
        f'{DIR_PREFIX}/serafim/{os.environ["DF_SHORT"]}/{os.environ["DF_LONG"]}/report.json'
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
        f'{DIR_PREFIX}/serafim/{os.environ["DF_SHORT"]}/{os.environ["DF_LONG"]}_good/nextclade.json',
        f'{DIR_PREFIX}/serafim/{os.environ["DF_SHORT"]}/{os.environ["DF_LONG"]}/nextclade.json'
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
        f'{DIR_PREFIX}/serafim/{os.environ["DF_SHORT"]}/{os.environ["DF_LONG"]}_good/nextclade.json',
        f'{DIR_PREFIX}/serafim/{os.environ["DF_SHORT"]}/{os.environ["DF_LONG"]}/nextclade.json'
    output:
        f'{DIR_PREFIX}/serafim/{os.environ["DF_SHORT"]}/{os.environ["DF_LONG"]}_good/report.json',
        f'{DIR_PREFIX}/serafim/{os.environ["DF_SHORT"]}/{os.environ["DF_LONG"]}/report.json'
    run:
        converter = f'{DETECTOR_PATH}/parse_nextclade_report.py'

        for i in range(2):
            shell(f'{converter} -i {input[i]} -o {output[i]}')
