#!/usr/bin/env python3
# @Author  : Serafim Nenarokov (serafim.nenarokov@gmail.com)

import os
import argparse
import re
import json

from tqdm import tqdm

from lib.utils import load_nextclade_aa_genes
from lib.utils import load_nextclade_na_alignments
from lib.utils import load_reference_aa_genes
from lib.utils import load_reference_na_genes

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

    out_path_help_str = """Path to the output file
    (report.json file in nextclade output path will be created if not specified)"""

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

MUTATION_REGEX = re.compile(
    '(?P<gene_name>\w+):(?P<source_aa>[a-zA-Z\-\*])(?P<index>\d+)(?P<target_aa>[a-zA-Z\-\*])')

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

def group_numbers_by_neighbours(numbers):
    groups = []

    current_group = []
    for i, current_num in enumerate(numbers):
        # first item
        if len(current_group) == 0:
            current_group.append(current_num)
            continue

        if current_num == current_group[-1]+1:
            current_group.append(current_num)
        else:  # we met a new number
            if len(current_group) > 1:
                groups.append(current_group)

            current_group = [current_num]

        if i == len(numbers) - 1:
            if len(current_group) > 1:
                groups.append(current_group)

    return groups

def main():
    arguments = parse_arguments()

    with open(arguments.report_path) as f:
        report = json.load(f)

    by_seq_report = {}
    # {'seq_id': {
    #   'gene_id': { site_id: {
    #                   code: "MUT_CODE",
    #                   type: "DEL/SUB"}, ...},
    #    ...},
    # ...}
    for mut_code, data in report.items():
        mut_parsed = re.match(MUTATION_REGEX, mut_code)

        if not mut_parsed:
            raise(f'Mutation `{mut_code}` cannot be parsed!')

        site_id = int(mut_parsed.group('index'))
        gene_name = mut_parsed.group('gene_name')

        if mut_parsed.group('target_aa') == '-':
            mut_type = 'DEL'
        else:
            mut_type = 'SUB'

        for seq_id in data['seq_ids']:
            if seq_id not in by_seq_report:
                by_seq_report[seq_id] = {}

            if gene_name not in by_seq_report[seq_id]:
                by_seq_report[seq_id][gene_name] = {}

            by_seq_report[seq_id][gene_name][site_id] = {
                'code': mut_code, 'type': mut_type}

    ## Find suspicious mutations
    susp_mutations = {}

    # {'mut_code': {'cnt': N, 'seq_ids': [...]}, ...}
    for seq_id, data in tqdm(by_seq_report.items()):
        for gene_name, mutations in data.items():
            site_ids = sorted(list(mutations.keys()))

            groups = group_numbers_by_neighbours(site_ids)

            susp_groups = []
            for group in groups:
                del_found = False
                sub_found = False

                for idx in group:
                    mut_type = mutations[idx]['type']
                    mut_code = mutations[idx]['code']

                    if mut_type == 'DEL':
                        del_found = True

                    if mut_type == 'SUB' and del_found:
                        sub_found = True
                        break

                if del_found and sub_found:
                    susp_groups.append(group)

            for susp_group in susp_groups:
                key = ','.join([mutations[idx]['code'] for idx in susp_group])


                if key not in susp_mutations:
                    susp_mutations[key] = {'seq_ids': [], 'cnt': 0}

                susp_mutations[key]['seq_ids'].append(seq_id)
                susp_mutations[key]['cnt'] += 1


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
        json.dump(susp_mutations, out_f, indent=4, sort_keys=False)

    print(f'Finished. {len(susp_mutations)} mutations found.')


if __name__ == '__main__':
    main()
