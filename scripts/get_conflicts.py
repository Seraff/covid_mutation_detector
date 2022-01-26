#!/usr/bin/env python3

import json
import re
import os
import argparse
from slow_test import load_genes

MUTATION_REGEX = re.compile(
    '(?P<gene_name>\w+):(?P<source_aa>[a-zA-Z\-\*])(?P<index>\d+)(?P<target_aa>[a-zA-Z\-\*])')

def parse_arguments():
    usage = "./generate_heatmap.py"
    description = """
    Generates report about differences in two mutation predictions
    """
    description = re.sub(r"(?<=\w)\n(?=\w)", "\n", description)
    formatter = argparse.RawTextHelpFormatter

    parser = argparse.ArgumentParser(usage,
                                     description=description,
                                     formatter_class=formatter)
    parser.add_argument("-r", "--ref",
                        required=True,
                        help="Reference sequences folder")
    parser.add_argument("--our_rep",
                        required=True,
                        help="Our report .json")
    parser.add_argument("--our_seqs",
                        required=True,
                        help="Our analysis sequences folder")
    parser.add_argument('--nextclade_rep',
                        required=True,
                        help="Nextclade report .json")
    parser.add_argument("--nextclade_seqs",
                        required=True,
                        help="Nextclade analysis sequences folder")
    parser.add_argument("--nextclade_seq_prefix",
                        required=False,
                        default='seq.gene.',
                        help="Nextclade analysis sequences folder")
    parser.add_argument('-o', '--output',
                        required=True,
                        help="Output report file in Markdown format")
    return parser.parse_args()

def main():
    arguments = parse_arguments()

    with open(arguments.our_rep) as f:
        our_report = json.load(f)

    with open(arguments.nextclade_rep) as f:
        nextclade_report = json.load(f)

    merged_report = {}

    all_muts = set(our_report.keys()).union(set(nextclade_report.keys()))

    for mut_code in all_muts:
        data = {'our_seq_ids': [], 'nextclade_seq_ids': []}

        if mut_code in our_report:
            data['our_seq_ids'] = our_report[mut_code]['seq_ids']

        if mut_code in nextclade_report:
            data['nextclade_seq_ids'] = nextclade_report[mut_code]['seq_ids']

        merged_report[mut_code] = data

    # Getting sequences
    # { 'SEQ0001': { 'ORF3a': 'MADSNGTI...', ...}, ...}
    ref_seqs = load_genes(arguments.ref, '')
    ref_seqs = ref_seqs[list(ref_seqs.keys())[0]]

    our_seqs = load_genes(arguments.our_seqs, 'gene.')
    nextclade_seqs = load_genes(
        arguments.nextclade_seqs, arguments.nextclade_seq_prefix)


    conflicts = []

    for mut_code, data in merged_report.items():
        us_not_them = set(data['our_seq_ids']) - set(data['nextclade_seq_ids'])
        them_not_us = set(data['nextclade_seq_ids']) - set(data['our_seq_ids'])

        if us_not_them != them_not_us:
            mut_matched = re.match(MUTATION_REGEX, mut_code)
            mut_index = mut_matched.group('index')
            gene_name = mut_matched.group('gene_name')

            if gene_name in {'ORF10', 'ORF14'}:
                continue

            for group in [us_not_them, them_not_us]:
                tp = 'ONLY_US' if group == us_not_them else 'ONLY_THEM'

                if len(group) > 0:
                    example_seq_id = list(group)[0]

                    conflict = {'mut_code': mut_code}
                    conflict['type'] = tp
                    conflict['seq_ids'] = list(group)
                    conflict['seq_cnt'] = len(group)
                    conflict['gene_name'] = gene_name
                    conflict['site_id'] = mut_index
                    conflict['ref_seq'] = ref_seqs[gene_name]
                    conflict['our_seq'] = our_seqs[example_seq_id][gene_name]
                    conflict['nxt_seq'] = nextclade_seqs[example_seq_id][gene_name]

                    conflicts.append(conflict)

    conflicts = sorted(conflicts, key=lambda x: (x['gene_name'], x['seq_cnt']), reverse=True)

    with open(arguments.output, 'w') as out_f:
        for conflict in conflicts:
            out_f.write(f"## {conflict['mut_code']}\n")
            out_f.write('\n')

            if conflict['type'] == 'ONLY_US':
                tp = '**Only we** have'
            else:
                tp = '**Only Nextclade** has'

            out_f.write(f"{tp} this mutation in {conflict['seq_cnt']} sequence(s).\n")
            out_f.write('```\n')
            out_f.write(f"REF: {conflict['ref_seq']}\n")
            out_f.write(f"OUR: {conflict['our_seq']}\n")
            out_f.write(f"NXT: {conflict['nxt_seq']}\n")
            out_f.write('```\n\n')


if __name__ == '__main__':
    main()
