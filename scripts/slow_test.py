#!/usr/bin/env python3

import os
import glob
import argparse
import re
import json

from Bio import SeqIO
from tqdm import tqdm

MUTATION_REGEX = re.compile(
    '(?P<gene_name>\w+):(?P<source_aa>[a-zA-Z\-\*])(?P<index>\d+)(?P<target_aa>[a-zA-Z\-\*])')
MASK_SYMBOL = 'X'

def parse_arguments():
    usage = "./slow_test.py"
    description = """
    The script checks report json file and nextalign sequences
    for correctness of mutation prediction.
    """
    description = re.sub(r"(?<=\w)\n(?=\w)", "\n", description)
    formatter = argparse.RawTextHelpFormatter

    parser = argparse.ArgumentParser(usage,
                                     description=description,
                                     formatter_class=formatter)

    parser.add_argument("-r", "--report",
                        required=True,
                        help="Report file (.json)")
    parser.add_argument('-s', '--seq_folder',
                        required=True,
                        help="Folder with reported sequences (per-gen files, NOT nextalign one)")
    return parser.parse_args()


def load_genes(dir_path, filename_prefix='gene.'):
    file_paths = glob.glob(os.path.join(
        dir_path, '', f"{filename_prefix}*.fasta"))

    if len(file_paths) == 0:
        print("No nextalign files found")
        exit(-1)

    GENE_FILENAME_REGEX = re.compile(re.escape(filename_prefix) + '(.*)\.fasta')

    result = {}

    for path in file_paths:
        filename = path.split('/')[-1]
        match = re.match(GENE_FILENAME_REGEX, filename)

        if not match:
            print(f"Cannot gene a gene name from path: {path}")
            exit(-1)

        gene_name = match.group(1)

        for record in SeqIO.parse(path, "fasta"):
            if record.id not in result:
                result[record.id] = {}

            result[record.id][gene_name] = str(record.seq)

    return result


def main():
    arguments = parse_arguments()

    with open(arguments.report) as f:
        report = json.loads(f.read())

    # { 'SEQ0001': { 'ORF3a': ['ORF3a:L123P', ...] }, ...}
    per_seq_mutations = {}

    for mut_str, data in report.items():
        for seq_id in data['seq_ids']:
            if seq_id not in per_seq_mutations:
                per_seq_mutations[seq_id] = {}

            gene = mut_str.split(':')[0]

            if gene not in per_seq_mutations[seq_id]:
                per_seq_mutations[seq_id][gene] = []

            per_seq_mutations[seq_id][gene].append(mut_str)

    # { 'SEQ0001': { 'ORF3a': 'MADSNGTI...', ...}, ...}
    genes = load_genes(arguments.seq_folder)

    ref_genes = load_genes('data/reference', '')
    ref_genes = ref_genes[list(ref_genes.keys())[0]]

    for seq_id, mutations in tqdm(per_seq_mutations.items()):
        for gene_id, mut_list in mutations.items():
            gene_seq = genes[seq_id][gene_id]
            ref_seq = ref_genes[gene_id]

            mutation_idxs = []

            # check if all mutations in the list present in the sequence
            for mut in mut_list:
                m = re.match(MUTATION_REGEX, mut)

                if not m:
                    print(f"Mutation {mut} cannot be parsed...")
                    exit(-1)

                mut_index = int(m.group('index'))
                source_aa = m.group('source_aa')
                target_aa = m.group('target_aa')

                mutation_idxs.append(mut_index-1)

                if target_aa.upper() != gene_seq[mut_index-1].upper():
                    print(f"Problem in mutation {mut} in sequence {seq_id}")
                    exit(-1)

            # check if the rest of the gene sequence is ok
            for i, aa in enumerate(gene_seq):
                aa = aa.upper()
                if aa == MASK_SYMBOL or i in mutation_idxs:
                    continue

                if aa != ref_seq[i]:
                    print(f"Problem with sequence {seq_id} on position {i+1}, not like in reference")
                    exit(-1)

    print("Done.")


if __name__ == '__main__':
    main()
