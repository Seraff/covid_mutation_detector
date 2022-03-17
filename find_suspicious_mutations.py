#!/usr/bin/env python3
# @Author  : Serafim Nenarokov (serafim.nenarokov@gmail.com)

import os
import argparse
import re
import json

from tqdm import tqdm

from lib.utils import MUTATION_REGEX
from lib.utils import load_nextclade_aa_genes
from lib.utils import load_nextclade_na_alignments
from lib.utils import load_reference_aa_genes
from lib.utils import load_reference_na_genes
from lib.utils import get_suspicious_mutations

def parse_arguments():
    usage = "parse_nextclade_report.py"
    description = """
    Parse nextclade report in COG json format and provide a list of suspicious mutations.
    """
    description = re.sub(r"(?<=\w)\n(?=\w)", "\n", description)
    formatter = argparse.RawTextHelpFormatter

    parser = argparse.ArgumentParser(usage,
                                     description=description,
                                     formatter_class=formatter)

    parser.add_argument("-i", "--report_path",
                        required=True,
                        help="Report in COG json format.")
    parser.add_argument("-d", "--nextclade_out_path",
                        required=True,
                        help="Path to nextclade output folder.")
    parser.add_argument('-o', '--output_path',
                        required=True,
                        help="Json file with a list of suspicious mutations")
    return parser.parse_args()


def expand_aa_seq(seq):
    '''
    MVANAH ->   M  V  A  N  A  H
    '''
    return ' ' + '   '.join(list(seq)) + ' '


def expand_na_seq(seq):
    '''
    GTACCT -> GTA CCT
    '''
    l = list(seq)
    grouped = [l[i:i+3] for i in range(0, len(l), 3)]
    return ' '.join([''.join(e) for e in grouped])


def main():
    arguments = parse_arguments()

    with open(arguments.report_path) as f:
        report = json.load(f)

    ## Find suspicious mutations
    susp_mutations = get_suspicious_mutations(report)

    ref_aa_genes = load_reference_aa_genes()
    ref_na_genes = load_reference_na_genes()
    aa_genes = load_nextclade_aa_genes(arguments.nextclade_out_path)
    na_genes = load_nextclade_na_alignments(arguments.nextclade_out_path)

    result = []
    ## Adding sequences
    for mut_codes, data in susp_mutations.items():
        mut_codes_str = mut_codes
        mut_codes = mut_codes.split(',')

        mut_left_parsed = re.match(MUTATION_REGEX, mut_codes[0])
        mut_right_parsed = re.match(MUTATION_REGEX, mut_codes[-1])
        gene_name = mut_left_parsed.group('gene_name')
        left_mut_id = int(mut_left_parsed.group('index'))
        right_mut_id = int(mut_right_parsed.group('index'))

        seq_id = data['seq_ids'][0]

        item = {}

        radius = 5
        left_id = left_mut_id - radius - 1

        if left_id < 0:
            left_id = 0

        right_id = right_mut_id + radius

        if right_id > len(aa_genes[seq_id][gene_name]):
            right_id = len(aa_genes[seq_id][gene_name]) - 1

        ref_aa_seq = ref_aa_genes[gene_name][left_id:right_id]
        ref_na_seq = ref_na_genes[gene_name][left_id*3:right_id*3]
        aa_seq = aa_genes[seq_id][gene_name][left_id:right_id]
        na_seq = na_genes[seq_id][gene_name][left_id*3:right_id*3]

        item['predict'] = mut_codes_str
        item['verdict'] = ''
        item['ref_aa'] = expand_aa_seq(ref_aa_seq)
        item['ref_na'] = expand_na_seq(ref_na_seq)
        stars = ['*' for _i in range(len(mut_codes))]

        left_aa_num = left_mut_id - left_id - 1
        right_aa_num = right_id - right_mut_id
        item['mut_id'] = f"{ ' '*left_aa_num*4 } { (' '*3).join(stars) } { ' '*right_aa_num*4 }"

        item['seq_aa'] = expand_aa_seq(aa_seq)
        item['seq_na'] = expand_na_seq(na_seq)

        item['notes'] = ''
        item['cnt'] = data['cnt']
        item['seq_ids'] = data['seq_ids']
        result.append(item)

    ## Output
    with open(arguments.output_path, 'w') as out_f:
        json.dump(result, out_f, indent=4, sort_keys=False)

    print(f'Finished. {len(susp_mutations)} situations found.')


if __name__ == '__main__':
    main()
