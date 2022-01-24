#!/usr/bin/env python3

import json
import re
import os
import argparse

def parse_arguments():
    usage = "./generate_heatmap.py"
    description = """
    Generates mutations heatmap from input table file.
    """
    description = re.sub(r"(?<=\w)\n(?=\w)", "\n", description)
    formatter = argparse.RawTextHelpFormatter

    parser = argparse.ArgumentParser(usage,
                                     description=description,
                                     formatter_class=formatter)

    parser.add_argument("-f", "--first",
                        required=True,
                        help="First report .json")
    parser.add_argument('-s', '--second',
                        required=True,
                        help="Second report .json")
    parser.add_argument('-o', '--output',
                        required=True,
                        help="Output report file")
    return parser.parse_args()

def main():
    arguments = parse_arguments()

    with open(arguments.first) as f:
        first_report = json.load(f)

    with open(arguments.second) as f:
        second_report = json.load(f)

    result = {}

    for mut_code, mut_data in first_report.items():
        if mut_code not in result:
            result[mut_code] = {'first': 0, 'second': 0}
        result[mut_code]['first'] = mut_data['cnt']

    for mut_code, mut_data in second_report.items():
        if mut_code not in result:
            result[mut_code] = {'first': 0, 'second': 0}
        result[mut_code]['second'] = mut_data['cnt']

    with open(arguments.output, 'w') as out_f:
        for mut_code, mut_data in result.items():
            out_f.write(f"{mut_code},{mut_data['first']},{mut_data['second']}\n")

if __name__ == '__main__':
    main()
