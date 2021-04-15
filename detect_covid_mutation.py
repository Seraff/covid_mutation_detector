#!/usr/bin/env python3
# @Date    : 2021-04-12
# @Author  : Serafim Nenarokov (serafim.nenarokov@gmail.com)

import os
import argparse
import re
from pathlib import Path
import pathlib

from Bio import SeqIO
from lib.nextaligner import Nextaligner
from cluster_sequences import Clusterizer

ROOT_PATH = str(Path(os.path.dirname(os.path.realpath(__file__))))
TMP_PATH = os.path.join(ROOT_PATH, "tmp")
REFERENCE_PATH = os.path.join(ROOT_PATH, "data", "sars_cov_reference_prot.fasta")
NEXTALIGN_BIN_PATH = os.getenv('NEXTALIGN_PATH') or 'nextalign'

def parse_arguments():
    usage = "./detect_covid_mutation.py"
    description = """
    The script filters finds COVID mutation in sequence from fasta file.
    """
    description = re.sub(r"(?<=\w)\n(?=\w)", "\n", description)
    formatter = argparse.RawTextHelpFormatter

    parser = argparse.ArgumentParser(usage,
                                     description=description,
                                     formatter_class=formatter)

    parser.add_argument("-f", "--fasta",
                        required=True,
                        help="Input fasta file (.fasta)")
    return parser.parse_args()

def check_mutations(record):
    pass

def main():
    args = parse_arguments()

    fasta_path = pathlib.Path(args.fasta)

    if not fasta_path.exists():
        print("Fasta file not found")
        exit(1)

    if not fasta_path.is_file():
        print("Fasta file is not a file")
        exit(1)

    # align sequences with Nextaligner

    clusterizer = Clusterizer(fasta_path)
    clusterizer.perform()
    clusters = clusterizer.get_cluster_info_list()

    for cluster in clusters:
        pass

    mutation_detector = MutationDetector()

    # run cluster_sequences.py, get JSON [{ id: 'blabla', seqs: ['MTITA...', ...]}]
    # detect and apply mutations to clusters, get a list of mutations per cluster

    # compare mutations with known mutations


if __name__ == '__main__':
    main()
