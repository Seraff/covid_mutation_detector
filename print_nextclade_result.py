#!/usr/bin/env python3
# @Author  : Serafim Nenarokov (serafim.nenarokov@gmail.com)

import argparse
from lib.utils import *


def parse_arguments():
    usage = "print_nextclade_result.py"
    description = """
    Pretty print reference gene & actual nextclade gene by it's COG df-number and sequence ID.
    Very useful for displaying the concrete situation, if it looks suspicious.

    Example:
    ./print_nextclade_result.py -r /disk/cogcz/serafim/df-20220422/datafreeze-2022-04-22_12_weeks/ -g S -s EPI_ISL_10132849
    """

    parser = argparse.ArgumentParser(usage,
                                     description=description,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-r", "--report_folder",
                        required=True,
                        help="Report folder (everythin before `report_fixed.json` in final report path.")
    parser.add_argument("-g", "--gene_id",
                        required=True,
                        help="Id of the gene (ex. `S`)")
    parser.add_argument("-s", "--seq_id",
                        required=True,
                        help="Id of the sequence (ex. `EPI_ISL_10132849`)")
    return parser.parse_args()


def main():
    arguments = parse_arguments()

    if arguments.gene_id not in COVID_GENES:
        print("Gene should be a valid COVID-19 gene.")
        exit(-1)

    # Loading sequences for the sample

    na_alignments = load_nextclade_na_alignments(arguments.report_folder)

    if arguments.seq_id not in na_alignments:
        print(f'Sequence {arguments.seq_id} not found \
                in Nextclade nucleotide alignments')
        exit(-1)

    sample_na_seq = na_alignments[arguments.seq_id][arguments.gene_id]

    aa_alignments = load_nextclade_aa_genes(arguments.report_folder)

    if arguments.seq_id not in aa_alignments:
        print(
            f'Sequence {arguments.seq_id} not found \
              in Nextclade nucleotide alignments')
        exit(-1)

    sample_aa_seq = aa_alignments[arguments.seq_id][arguments.gene_id]

    # Loading sequences for reference

    ref_aa_seq = load_reference_aa_genes()[arguments.gene_id]
    ref_na_seq = load_reference_na_genes()[arguments.gene_id]

    # preparing mutations line (with '*' symbol)

    mut_line = ''
    for i, ref_aa in enumerate(ref_aa_seq):
        if ref_aa == sample_aa_seq[i]:
            symbol = ' '
        else:
            symbol = '*'

        mut_line += symbol

    # Printhing the output

    print('# ' + arguments.seq_id)
    print('ref_na: ' + pprint_na_seq(ref_na_seq))
    print('ref_aa: ' + pprint_aa_seq(ref_aa_seq))
    print('mut_id: ' + pprint_aa_seq(mut_line))
    print('seq_na: ' + pprint_na_seq(sample_na_seq))
    print('seq_aa: ' + pprint_aa_seq(sample_aa_seq))


if __name__ == '__main__':
    main()
