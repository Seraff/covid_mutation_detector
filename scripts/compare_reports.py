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

    result = []
    for mut_code, mut_data in first_report.items():
        row = [mut_code, mut_data['cnt']]
        if mut_code in second_report:
            row.append(second_report[mut_code]['cnt'])
        else:
            row.append(0)
        result.append(row)

    with open(arguments.output, 'w') as out_f:
        for row in result:
            out_f.write(','.join([str(e) for e in row]) + '\n')

if __name__ == '__main__':
    main()
