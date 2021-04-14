#!/usr/bin/env python3
# @Date    : 2021-04-14
# @Author  : Serafim Nenarokov (serafim.nenarokov@gmail.com)

import os
import argparse
import re
import pathlib
import json

from Bio import SeqIO


MASK_SYMBOL = 'X'


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

    parser.add_argument("-f", "--fasta",
                        required=True,
                        help="Input fasta file (.fasta)")
    parser.add_argument("-o", "--output_folder_path",
                        required=True,
                        help="Folder for output files (will be created if doesn't exist)")
    return parser.parse_args()


def belongs_to_one_cluster(query, target):
    if query == target:
        return True

    for i, aa in enumerate(query):
        query_aa = query[i].upper()
        target_aa = target[i].upper()

        if query_aa == 'X' or target_aa == 'X':
            continue

        if query_aa != target_aa:
            return False

    return True


def main():
    args = parse_arguments()

    fasta_path = pathlib.Path(args.fasta)

    if not fasta_path.exists():
        print("Fasta file not found")
        exit(1)

    if not fasta_path.is_file():
        print("Fasta file is not a file")
        exit(1)

    clusters = []

    with open(fasta_path) as fasta_f:
        for record in SeqIO.parse(fasta_f, 'fasta'):
            cluster_found = False
            seq_str = str(record.seq)

            for cluster in clusters:
                cluster_seq_str = str(cluster[0].seq)

                if belongs_to_one_cluster(cluster_seq_str, seq_str):
                    cluster.append(record)
                    cluster_found = True
                    break

            if not cluster_found:
                clusters.append([record])

    print(f"Found {len(clusters)} clusters")

    clusters.sort(key=lambda x: len(x), reverse=True)

    pathlib.Path(args.output_folder_path).mkdir(parents=True, exist_ok=True)
    txt_path = os.path.join(args.output_folder_path, 'clusters.txt')
    json_path = os.path.join(args.output_folder_path, 'clusters.json')

    with open(txt_path, 'w') as out_f:
        for i, cluster in enumerate(clusters):
            out_f.write(f"Cluster {i+1}\n")

            for record in cluster:
                out_f.write(f"{str(record.seq)} - {record.id}\n")

    result = {'clusters': []}

    with open(json_path, 'w') as out_f:
        for i, cluster in enumerate(clusters):
            min_X_seq = None

            for seq in cluster:
                seq_str = str(seq)
                sec_X_cnt = str(seq.seq).count('X')

                if str(min_X_seq) == 'None':
                    min_X_seq = seq
                else:
                    if sec_X_cnt < str(min_X_seq.seq).count('X'):
                        min_X_seq = seq

            entry = {'id': min_X_seq.id, 'seq': str(min_X_seq.seq)}
            result['clusters'].append(entry)

        out_f.write(json.dumps(result, indent=4))


if __name__ == '__main__':
    main()
