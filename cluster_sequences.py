#!/usr/bin/env python3
# @Date    : 2021-04-14
# @Author  : Serafim Nenarokov (serafim.nenarokov@gmail.com)

import os
import argparse
import re
import pathlib
import json

from Bio import SeqIO

ROOT_PATH = str(pathlib.Path(os.path.dirname(os.path.realpath(__file__))))
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
    o_hlp = "Folder for output files (will be created if doesn't exist)"

    parser.add_argument("-f", "--fasta",
                        required=True,
                        help="Input fasta file (.fasta)")
    parser.add_argument("-o", "--output_folder_path",
                        required=True,
                        help=o_hlp)
    return parser.parse_args()


class Clusterizer:
    def __init__(self, fasta_path):
        self.fasta_path = fasta_path

        self.clusters = []

    def perform(self):
        with open(self.fasta_path) as fasta_f:
            for record in SeqIO.parse(fasta_f, 'fasta'):
                seq_str = str(record.seq).upper()

                seq_str = self._remove_gaps(seq_str)

                cluster_found = False
                for cluster in self.clusters:
                    cluster_seq_str = self._remove_gaps(str(cluster[0].seq))

                    if self._belongs_to_one_cluster(cluster_seq_str, seq_str):
                        cluster.append(record)
                        cluster_found = True
                        break

                if not cluster_found:
                    self.clusters.append([record])

        print(f"Found {len(self.clusters)} clusters")

        self.clusters.sort(key=lambda x: len(x), reverse=True)

    def save_files(self, output_folder_path):
        if not self._is_performed():
            return False

        pathlib.Path(output_folder_path).mkdir(parents=True, exist_ok=True)
        txt_path = os.path.join(output_folder_path, 'clusters.txt')
        json_path = os.path.join(output_folder_path, 'clusters.json')

        # making cluster output file:
        # Cluster <id>
        # <sequence> - <seq id>
        # ...

        with open(txt_path, 'w') as out_f:
            for i, cluster in enumerate(self.clusters):
                out_f.write(f"Cluster {i+1}\n")

                for record in cluster:
                    out_f.write(f"{str(record.seq)} - {record.id}\n")



        result = {'clusters': self.get_cluster_info_list()}

        with open(json_path, 'w') as out_f:
            out_f.write(json.dumps(result, indent=4))

    def get_cluster_info_list(self):
        '''
        making list with clusters, in each cluster id and seq belong to
        a guy with smallest amount of X
        '''

        result = []

        for i, cluster in enumerate(self.clusters):
            min_X_seq = None

            for seq in cluster:
                seq_str = str(seq)
                sec_X_cnt = str(seq.seq).count(MASK_SYMBOL)

                if str(min_X_seq) == 'None':
                    min_X_seq = seq
                else:
                    if sec_X_cnt < str(min_X_seq.seq).count(MASK_SYMBOL):
                        min_X_seq = seq

            entry = {'id': min_X_seq.id,
                     'seq': self._remove_gaps(str(min_X_seq.seq)),
                     'orig_seq': str(min_X_seq.seq)}
            result.append(entry)

        return result

    def _belongs_to_one_cluster(self, query, target):
        if query == target:
            return True

        for i, aa in enumerate(query):
            query_aa = query[i].upper()
            target_aa = target[i].upper()

            if query_aa == MASK_SYMBOL or target_aa == MASK_SYMBOL:
                continue

            if query_aa != target_aa:
                return False

        return True

    def _is_performed(self):
        return len(self.clusters) != 0

    def _remove_gaps(self, sequence):
        return sequence.replace('-', '')


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
