#!/usr/bin/env python3
# @Date    : 2021-04-14
# @Author  : Serafim Nenarokov (serafim.nenarokov@gmail.com)

import os
import argparse
import re
import pathlib
import json

from Bio import SeqIO
from lib.clusterizer import Clusterizer

ROOT_PATH = str(pathlib.Path(os.path.dirname(os.path.realpath(__file__))))


def parse_arguments():
    usage = "./detect_covid_mutation.py"
    description = """
    Sequence clustering with 100%.
    """
    description = re.sub(r"(?<=\w)\n(?=\w)", "\n", description)
    formatter = argparse.RawTextHelpFormatter

    parser = argparse.ArgumentParser(usage,
                                     description=description,
                                     formatter_class=formatter)
    o_hlp = "Folder for output files (will be created if doesn't exist)"

    parser.add_argument("-f", "--fasta",
                        required=True,
                        help="Input fasta file (.fasta)")
    parser.add_argument("-o", "--output_folder_path",
                        required=True,
                        help=o_hlp)
    return parser.parse_args()



def main():
    args = parse_arguments()
    fasta_path = pathlib.Path(args.fasta)

    if not fasta_path.exists():
        print("Fasta file not found")
        exit(1)

    if not fasta_path.is_file():
        print("Fasta file is not a file")
        exit(1)

    clusterizer = Clusterizer(fasta_path)
    clusterizer.perform()
    clusterizer.save_files(args.output_folder_path)


if __name__ == '__main__':
    main()
