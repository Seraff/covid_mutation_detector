#!/usr/bin/env python3
# @Author  : Serafim Nenarokov (serafim.nenarokov@gmail.com)

import os
import argparse
import re
import json

from tqdm import tqdm

AA_SUBSTITUTIONS_KEY = 'aaSubstitutions'
AA_DELETIONS_KEY = 'aaDeletions'
AA_INSERTIONS_KEY = 'aaInsertions'

def parse_arguments():
    usage = "parse_nextclade_report.py"
    description = """
    Converts nextclade report to the report in COG format.
    """
    description = re.sub(r"(?<=\w)\n(?=\w)", "\n", description)
    formatter = argparse.RawTextHelpFormatter

    parser = argparse.ArgumentParser(usage,
                                     description=description,
                                     formatter_class=formatter)

    out_path_help_str = """Path to the output file
    (report.json file in nextclade output path will be created if not specified)"""

    parser.add_argument("-i", "--nextclade_json_path",
                        required=True,
                        help="Nextclade output JSON file")
    parser.add_argument('-o', '--output_file',
                        help=out_path_help_str)
    return parser.parse_args()

def append_mutation_to_report(report, mut_code, seq_id):
    if mut_code not in report:
        report[mut_code] = {'cnt': 0, 'seq_ids': []}

    report[mut_code]['seq_ids'].append(seq_id)
    report[mut_code]['seq_ids'] = list(set(report[mut_code]['seq_ids']))
    report[mut_code]['cnt'] += 1

def main():
    arguments = parse_arguments()

    with open(arguments.nextclade_json_path) as f:
        nextclade_report = json.load(f)

    cog_report = {}

    for seq_data in tqdm(nextclade_report['results']):
        seq_id = seq_data['seqName']

        if AA_SUBSTITUTIONS_KEY in seq_data:
            for sub in seq_data[AA_SUBSTITUTIONS_KEY]:
                codon = int(sub['codon'])+1
                mut_code = f"{sub['gene']}:{sub['refAA']}{codon}{sub['queryAA']}"
                append_mutation_to_report(cog_report, mut_code, seq_id)

        if AA_DELETIONS_KEY in seq_data:
            for deletion in seq_data[AA_DELETIONS_KEY]:
                codon = int(deletion['codon'])+1
                mut_code = f"{deletion['gene']}:{deletion['refAA']}{codon}-"
                append_mutation_to_report(cog_report, mut_code, seq_id)

    out_path = arguments.output_file

    if not out_path:
        dir_path = os.path.dirname(arguments.nextclade_json_path)
        out_path = os.path.join(dir_path, 'report.json')

    with open(out_path, 'w') as f:
        json.dump(cog_report, f, indent=4, sort_keys=True)

    print("Done.")

if __name__ == "__main__":
    main()
