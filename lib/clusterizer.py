#!/usr/bin/env python3
# @Date    : 2021-04-14
# @Author  : Serafim Nenarokov (serafim.nenarokov@gmail.com)

import os
import pathlib
import json

from Bio import SeqIO

MASK_SYMBOL = 'X'
GAP_SYMBOL = '-'
END_SYMBOL = '*'

class Clusterizer:
    def __init__(self, fasta_path):
        self.fasta_path = fasta_path

        self.clusters = []

    def perform(self):
        with open(self.fasta_path) as fasta_f:
            processed_seq_headers = []

            for record in SeqIO.parse(fasta_f, 'fasta'):
                if record.description in processed_seq_headers:
                    continue

                record.seq = self._fix_incomplete_tail(record.seq)

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

                processed_seq_headers.append(record.description)

        self.clusters.sort(key=lambda x: len(x), reverse=True)

    def save_files(self, output_folder_path):
        if not self._is_performed():
            return False

        pathlib.Path(output_folder_path).mkdir(parents=True, exist_ok=True)
        csv_path = os.path.join(output_folder_path, 'clusters.csv')
        json_path = os.path.join(output_folder_path, 'clusters.json')

        with open(csv_path, 'w') as out_f:
            out_f.write("cluster_id,seq_id,seq")
            for i, cluster in enumerate(self.clusters):
                for record in cluster:
                    out_f.write(f"{i+1},{record.description},{str(record.seq)}\n")

        result = {'clusters': self.get_cluster_info_list()}

        with open(json_path, 'w') as out_f:
            out_f.write(json.dumps(result, indent=4))

    def get_cluster_info_list(self):
        '''
        making list with clusters, in each cluster id and seq belong to
        a sequence with smallest amount of X
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

            all_seqs = [s.description for s in cluster]

            entry = {
                     'name': f"Cluster {i}",
                     'seq_id': min_X_seq.description,
                     'seq': self._remove_gaps(str(min_X_seq.seq)),
                     'orig_seq': str(min_X_seq.seq),
                     'all_seq_ids': all_seqs,
                     'size': len(cluster)}
            result.append(entry)

        return result


    def _belongs_to_one_cluster(self, query, target):
        if query == target:
            return True

        if len(query) != len(target):
            return False

        for i, aa in enumerate(query):
            query_aa = query[i].upper()
            target_aa = target[i].upper()

            if query_aa == MASK_SYMBOL or target_aa == MASK_SYMBOL:
                continue

            if query_aa != target_aa:
                return False

        return True

    def _fix_incomplete_tail(self, seq):
        '''
        Gets string sequence, replaces --- at the end with XXX
        '''

        stripped = seq.rstrip(GAP_SYMBOL)
        cnt = len(seq) - len(stripped)

        return stripped + (MASK_SYMBOL*cnt)

    def _is_performed(self):
        return len(self.clusters) != 0

    def _remove_gaps(self, sequence):
        return sequence.replace(GAP_SYMBOL, '')
