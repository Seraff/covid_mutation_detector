#!/usr/bin/env python3

import json
import re
import os
import argparse
from slow_test import load_genes
from pathlib import Path
from pprint import pprint

from Bio import SeqIO

MUTATION_REGEX = re.compile(
    '(?P<gene_name>\w+):(?P<source_aa>[a-zA-Z\-\*])(?P<index>\d+)(?P<target_aa>[a-zA-Z\-\*])')
ROOT_PATH = str(Path(os.path.dirname(os.path.realpath(__file__))).parent)
REF_AA_PATH = os.path.join(ROOT_PATH, 'data', 'reference')
REF_NA_PATH = os.path.join(ROOT_PATH, 'data', 'sars_cov_reference.fasta')
GENE_COORDS = {
    "E": (26245, 26472),
    "M": (26523, 27191),
    "N": (28274, 29533),
    "ORF1a": (266, 13468),
    "ORF1b": (13468, 21555),
    "ORF3a": (25393, 26220),
    "ORF6": (27202, 27387),
    "ORF7a": (27394, 27759),
    "ORF7b": (27756, 27887),
    "ORF8": (27894, 28259),
    "ORF9b": (28284, 28577),
    "S": (21563, 25384)
}

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
    parser.add_argument("--nextclade_alignments_fasta",
                        required=True,
                        help="Nextclade analysis alignment fasta file")
    parser.add_argument('-o', '--output',
                        required=True,
                        help="Output report file in Markdown format")
    return parser.parse_args()

def expand_aa_seq(seq):
    '''
    MVANAH ->   M  V  A  N  A  H
    '''
    return ' ' + '  '.join(list(seq)) + ' '

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
    ref_aa_seqs = load_genes(REF_AA_PATH, '')
    ref_aa_seqs = ref_aa_seqs[list(ref_aa_seqs.keys())[0]]

    our_seqs = load_genes(arguments.our_seqs, 'gene.')
    nextclade_seqs = load_genes(
        arguments.nextclade_seqs, arguments.nextclade_seq_prefix)

    ref_na_seqs = {}
    # Get reference NA sequences
    ref_rec = SeqIO.read(REF_NA_PATH, "fasta")
    for gene_name, coords in GENE_COORDS.items():
        ref_na_seqs[gene_name] = str(ref_rec.seq[coords[0]-1:coords[1]])

    nextclade_na_seqs = {}
    for rec in SeqIO.parse(arguments.nextclade_alignments_fasta, 'fasta'):
        if rec.id not in nextclade_na_seqs:
            nextclade_na_seqs[rec.id] = {}

        for gene_name, coords in GENE_COORDS.items():
            sequence = str(rec.seq[coords[0]-1:coords[1]])
            nextclade_na_seqs[rec.id][gene_name] = sequence

    # Checking conflicts

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
                tp = 'ONLY_US' if group == us_not_them else 'ONLY_NEXTCLADE'

                if len(group) > 0:
                    example_seq_id = list(group)[0]

                    conflict = {'mut_code': mut_code}
                    conflict['type'] = tp
                    conflict['seq_ids'] = list(group)
                    conflict['seq_cnt'] = len(group)
                    conflict['gene_name'] = gene_name
                    conflict['site_id'] = int(mut_index)
                    conflict['ref_seq'] = ref_aa_seqs[gene_name]
                    conflict['our_seq'] = our_seqs[example_seq_id][gene_name]
                    conflict['nxt_seq'] = nextclade_seqs[example_seq_id][gene_name]

                    conflict['ref_na_seq'] = ref_na_seqs[gene_name]
                    conflict['nxt_na_seq'] = nextclade_na_seqs[example_seq_id][gene_name]

                    conflicts.append(conflict)

    conflicts = sorted(conflicts, key=lambda x: (
        x['gene_name'], x['seq_cnt']), reverse=True)

    with open(arguments.output, 'w') as out_f:
        top_conflicts = sorted(conflicts, key=lambda x: x['seq_cnt'], reverse=True)[:15]
        out_f.write('## Top conflicts\n')
        out_f.write('\n')
        for c in top_conflicts:
            out_f.write(f"* {c['mut_code']} ({c['type'].lower().replace('_', ' ')}) - {c['seq_cnt']} sequences \n")

        out_f.write('\n')

        for conflict in conflicts:
            left_cut_len = len(conflict['our_seq']) - len(conflict['our_seq'].lstrip('-'))

            if left_cut_len > 0 and conflict['site_id'] <= left_cut_len and conflict['mut_code'][-1] == '-':
                continue

            out_f.write(f"## Mutation {conflict['mut_code']}\n")
            out_f.write('\n')

            if conflict['type'] == 'ONLY_US':
                tp = '**Only we** have'
            else:
                tp = '**Only Nextclade** has'

            out_f.write(f"{tp} this mutation in {conflict['seq_cnt']} sequence(s).\n")
            out_f.write('```\n')
            out_f.write(f"REF: {conflict['ref_na_seq']}\n")
            out_f.write(f"REF: {expand_aa_seq(conflict['ref_seq'])}\n")
            before = conflict['site_id'] - 1
            after = len(conflict['ref_seq']) - before - 1
            out_f.write(f"{' '*5}{' '*before*3} ^ {' '*after*3}\n")
            out_f.write(f"OUR: {expand_aa_seq(conflict['our_seq'])}\n")
            out_f.write("\n")
            out_f.write(f"NXT: {conflict['nxt_na_seq']}\n")
            out_f.write(f"NXT: {expand_aa_seq(conflict['nxt_seq'])}\n")
            out_f.write('```\n\n')


if __name__ == '__main__':
    main()
