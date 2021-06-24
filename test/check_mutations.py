#!/usr/bin/env python3

import os
import json
from pathlib import Path

from termcolor import colored
from Bio import SeqIO

GENE_LIST = ("E", "M", "N", "ORF10",
             "ORF14", "ORF1a", "ORF1b",
             "ORF3a", "ORF6", "ORF7a",
             "ORF7b", "ORF8", "ORF9b", "S")

ROOT_PATH = str(Path(os.path.dirname(os.path.realpath(__file__))))
OUTPUT_PATH = os.path.join(ROOT_PATH, 'data', 'output')
REPORT_PATH = os.path.join(OUTPUT_PATH, 'report.json')

def load_gene_seqs():
    result = {}

    for gene in GENE_LIST:
        seqs = SeqIO.parse(os.path.join(OUTPUT_PATH, f"gene.{gene}.fasta"), 'fasta')
        result[gene] = [str(seq.seq) for seq in seqs]

    return result

def get_position_stats(position, seq_list):
    '''
    Position in mutation format, starts with 1
    '''
    if position <= 0:
        return None

    site = [seq[position-1:position] for seq in seq_list]

    result = {}
    for aa in site:
        if aa not in result:
            result[aa] = 0

        result[aa] += 1

    return result

def parse_mutation(mut_id):
    result = {}

    gene, data = mut_id.split(':')
    result['gene'] = gene
    result['from'] = data[0]
    result['to'] = data[-1]
    result['pos'] = data[1:-1]

    return result

def main():
    gene_seqs = load_gene_seqs()

    report = json.loads(open(REPORT_PATH).read())
    for mut_id, data in report.items():
        mut_data = parse_mutation(mut_id)
        report_cnt = data['cnt']

        stats = get_position_stats(int(mut_data['pos']), gene_seqs[mut_data['gene']])
        success = report_cnt == stats[mut_data['to']]
        if not success:
            print(colored(f"{mut_id}: BAD", 'red'))
        else:
            print(colored(f"{mut_id}: OK", 'green'))
    # print(get_position_stats(203, gene_seqs['N']))

if __name__ == '__main__':
    main()
