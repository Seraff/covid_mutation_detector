#!/usr/bin/env python3
# @Author: Serafim Nenarokov (serafim.nenarokov@gmail.com)

from ast import arguments
import os
import argparse
import subprocess
import json

def parse_arguments():
    usage = "proccess_datafreeze.py"
    description = """
    Proccess specified datafreeze. Exits with error is any of the pipelines fails.
    """
    formatter = argparse.RawTextHelpFormatter

    parser = argparse.ArgumentParser(usage,
                                     description=description,
                                     formatter_class=formatter)

    parser.add_argument("--df-short",
                        required=True,
                        help="Short datafreeze name (i.e. df-20220826)")
    parser.add_argument("--df-long",
                        required=True,
                        help="Long datafreeze name (i.e. datafreeze-2022-08-26)")
    return parser.parse_args()


def try_to_print_suspicious(output_path):
    suspicious_new_path = os.path.join(output_path, f'suspicious_new.json')

    if os.path.exists(suspicious_new_path):
        with open(suspicious_new_path) as f:
            raw_content = f.read()
            content = json.loads(raw_content)

        if len(content) > 0:
            print("Suspicious mutations found. Please resolve the following situations.\n")
            print(raw_content)


def main():
    arguments = parse_arguments()

    cog_root = "/storage/vestec1-elixir/projects/cogcz"
    df_path_prefix = f"{cog_root}/seq/{arguments.df_short}/{arguments.df_long}"
    postfixes = ('', '_12_weeks', '_12_weeks_good')

    input_paths = [f"{df_path_prefix}{pfx}.fasta.xz" for pfx in postfixes]
    metadata_paths = [f"{df_path_prefix}{pfx}.tsv" for pfx in postfixes]

    out_root = "/auto/vestec1-elixir/projects/cogcz/serafim"
    out_prefix = f"{out_root}/{arguments.df_short}/{arguments.df_long}"
    output_paths = [f"{out_prefix}{pfx}" for pfx in postfixes]

    for i, input_path in enumerate(input_paths):
        output_path = output_paths[i]
        metadata_path = metadata_paths[i]

        cmd = f"INPUT_PATH={input_path} OUTPUT_PATH={output_path} METADATA_PATH={metadata_path}"
        cmd += " snakemake -c1 --use-conda"
        rtn = subprocess.call(cmd, shell=True)

        if rtn != 0:
            print(f"One of the pipelines failed (`{cmd}`)")
            try_to_print_suspicious(output_path)
            exit(-1)


    print('Everything finished successfully!')
    print('Slack message:\n')
    msg = f""":bar_chart: Statistics for the datafreeze `{arguments.df_short}` is ready.

*Datafreeze*
`{ f"{output_paths[0]}/report.json" }`
`{ f"{output_paths[0]}/ratio_table_weeks.reduced.csv" }`
`{ f"{output_paths[0]}/ratio_table_regions.reduced.csv" }`

*12 Weeks*
`{ f"{output_paths[1]}/report.json" }`
`{ f"{output_paths[1]}/ratio_table_weeks.reduced.csv" }`
`{ f"{output_paths[1]}/ratio_table_regions.reduced.csv" }`

*12 Weeks Good*
`{ f"{output_paths[2]}/report.json" }`
`{ f"{output_paths[2]}/ratio_table_weeks.reduced.csv" }`
`{ f"{output_paths[2]}/ratio_table_regions.reduced.csv" }`
    """

    print(msg)


if __name__ == '__main__':
    main()
