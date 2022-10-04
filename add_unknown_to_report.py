#!/usr/bin/env python3
# @Author  : Serafim Nenarokov (serafim.nenarokov@gmail.com)

import os
import argparse
import re
import json

from tqdm import tqdm

from lib.utils import load_nextclade_aa_genes
from lib.utils import mut_id_to_dict

UNKNOWN_SYMBOL = 'X'

def parse_arguments():
    usage = "add_unknown_to_report.py"
    description = """
    Takes nextclade output, adds to report the list of `seq_ids` with
    X symbols on the mutation position. We are not sure, if there is a mutation
    or not.
    """
    description = re.sub(r"(?<=\w)\n(?=\w)", "\n", description)
    formatter = argparse.RawTextHelpFormatter

    parser = argparse.ArgumentParser(usage,
                                     description=description,
                                     formatter_class=formatter)

    output_report_path_help = """
    Path to the output report JSON file, report will be overwritten if not specified
    """

    parser.add_argument("-n", "--nextclade_dir_path",
                        required=True,
                        help="Nextclade output directory (with per-gene files)")
    parser.add_argument("-r", "--report_path",
                        required=True,
                        help="Path to the report JSON file")
    parser.add_argument("-o", "--output_report_path",
                        help=output_report_path_help.strip())

    return parser.parse_args()


def main():
    arguments = parse_arguments()

    with open(arguments.report_path) as f:
        report = json.load(f)

    # { record_id: { gene_name: seq, ... }, ... }
    sequences = load_nextclade_aa_genes(arguments.nextclade_dir_path)

    # Invert the sequences dictionary
    # to get { gene_name: { record_id: seq, ... }, ... }
    sequences_by_gene = {}
    for seq_id, data in sequences.items():
        for gene_name, seq in data.items():
            if gene_name not in sequences_by_gene:
                sequences_by_gene[gene_name] = {}
            sequences_by_gene[gene_name][seq_id] = seq

    # Inspecting mutation by mutation
    for mut_id, data in tqdm(report.items()):
        mut_dict = mut_id_to_dict(mut_id)

        if not mut_dict:
            raise ValueError(f"Cannot parse mutation {mut_id}")

        gene = mut_dict['gene_name']
        position = int(mut_dict['index'])

        unknown_seq_ids = []
        for seq_id, seq in sequences_by_gene[gene].items():
            if seq[position-1] == UNKNOWN_SYMBOL:
                unknown_seq_ids.append(seq_id)

        data['unknown_seq_ids'] = unknown_seq_ids

    out_path = arguments.output_report_path
    if not out_path:
        out_path = arguments.report_path

    with open(out_path, 'w') as out_f:
        json.dump(report, out_f, indent=4, sort_keys=True)


if __name__ == "__main__":
    main()
