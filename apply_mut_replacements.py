#!/usr/bin/env python3
# @Author  : Serafim Nenarokov (serafim.nenarokov@gmail.com)

import os
import argparse
import re
import json

from lib.utils import MUT_REPLACEMENT_RULES_PATH
from lib.utils import get_suspicious_mutations

def parse_arguments():
    usage = "apply_mut_replacements.py"

    description = """
    Apply replacement rules to report file.
    """

    description = re.sub(r"(?<=\w)\n(?=\w)", "\n", description)
    formatter = argparse.RawTextHelpFormatter

    parser = argparse.ArgumentParser(usage,
                                     description=description,
                                     formatter_class=formatter)

    parser.add_argument("-i", "--report_path",
                        required=True,
                        help="Report in COG json format")
    parser.add_argument("-r", "--rules_path",
                        help="Rules json file (optional)")
    parser.add_argument('-o', '--output_path',
                        required=True,
                        help="Modified report path")
    return parser.parse_args()

def main():
    arguments = parse_arguments()

    with open(arguments.report_path) as f:
        report = json.load(f)

    rules_path = arguments.rules_path

    if not rules_path:
        rules_path = MUT_REPLACEMENT_RULES_PATH

    with open(rules_path) as f:
        rules = json.load(f)

    rules_by_mut_codes = {}
    for rule in rules:
        rules_by_mut_codes[rule['predict']] = rule

    ## Find suspicious mutations
    susp_mutations = get_suspicious_mutations(report)

    undefined_rules = []

    for mut_codes, data in susp_mutations.items():
        if mut_codes not in rules_by_mut_codes:
            undefined_rules.append(mut_codes)
        else:
            rule = rules_by_mut_codes[mut_codes]
            seq_ids = data['seq_ids']

            old_mut_list = rule['predict'].split(',')
            new_mut_list = rule['verdict'].split(',')

            # Removing from report
            for mut in old_mut_list:
                for seq_id in seq_ids:
                    report[mut]['seq_ids'].remove(seq_id)
                    report[mut]['cnt'] -= 1

                    if len(report[mut]['seq_ids']) == 0:
                        del report[mut]

            # Adding to report
            for mut in new_mut_list:
                if mut not in report:
                    report[mut] = {'seq_ids': [], 'cnt': 0}

                for seq_id in seq_ids:
                    if seq_id in report[mut]['seq_ids']:
                        raise(f'Something wrong, {seq_id} is already in report')
                    else:
                        report[mut]['seq_ids'].append(seq_id)
                        report[mut]['cnt'] += 1

    if len(undefined_rules) == 0:
        print('All the suspicions satisfied')
        with open(arguments.output_path, 'w') as out_f:
            json.dump(report, out_f, indent=4, sort_keys=True)
    else:
        print(f'There are {len(undefined_rules)} unsatisfied suspicious patterns. Please, check.')
        for mut_codes in undefined_rules:
            print(mut_codes)


if __name__ == '__main__':
    main()
