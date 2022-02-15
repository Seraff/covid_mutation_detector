#!/usr/bin/env python3
# @Author  : Serafim Nenarokov (serafim.nenarokov@gmail.com)

import os
import argparse
import re
import json

from tqdm import tqdm

from lib.utils import load_nextclade_aa_genes
from lib.utils import load_nextclade_na_alignments
from lib.utils import load_reference_genes

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

    susp_mutations = {}
    # {'mut_code': {'cnt': N, 'seq_ids': [...]}, ...}
    for seq_id, data in tqdm(by_seq_report.items()):
        for gene_name, mutations in data.items():
            site_ids = sorted(list(mutations.keys()))

            del_found = False
            last_id = None
            for site_id in site_ids:
                mut_type = mutations[site_id]['type']

                # If first or new cluster of mutations
                if last_id == None or site_id != last_id+1:
                    del_found = mut_type == 'DEL'

                elif mut_type == 'DEL':
                    del_found = True

                elif mut_type == 'SUB' and del_found:
                    code = mutations[site_id]['code']

                    if code not in susp_mutations:
                        susp_mutations[code] = {'cnt': 0, 'seq_ids': []}

                    susp_mutations[code]['seq_ids'].append(seq_id)
                    susp_mutations[code]['cnt'] += 1

                last_id = site_id

    ref_aa_genes = load_reference_genes()
    aa_genes = load_nextclade_aa_genes(arguments.nextclade_out_path)
    na_genes = load_nextclade_na_alignments(arguments.nextclade_out_path)

    ## Adding sequences
    for mut_code, data in susp_mutations.items():
        mut_parsed = re.match(MUTATION_REGEX, mut_code)
        gene_name = mut_parsed.group('gene_name')
        site_id = int(mut_parsed.group('index'))
        seq_id = data['seq_ids'][0]

        result = {}

        radius = 10
        left_id = site_id - radius - 1
        if left_id < 0:
            left_id = 0

        ref_aa_seq = ref_aa_genes[gene_name][left_id:site_id+radius]
        aa_seq = aa_genes[seq_id][gene_name][left_id:site_id+radius]
        na_seq = na_genes[seq_id][gene_name][left_id*3:(site_id+radius)*3]

        result['ref'] = expand_aa_seq(ref_aa_seq)
        result['gen'] = expand_aa_seq(aa_seq)
        result['mut'] = f"{' '*radius*4} ^ {' '*radius*4}"
        result['nuc'] = expand_na_seq(na_seq)

        result['cnt'] = data['cnt']
        result['seq_ids'] = data['seq_ids']
        susp_mutations[mut_code] = result

    ## Output
    with open(arguments.output_path, 'w') as out_f:
        json.dump(susp_mutations, out_f, indent=4, sort_keys=False)

    print(f'Finished. {len(susp_mutations)} mutations found.')


if __name__ == '__main__':
    main()
