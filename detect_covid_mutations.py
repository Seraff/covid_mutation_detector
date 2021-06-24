#!/usr/bin/env python3
# @Date    : 2021-04-12
# @Author  : Serafim Nenarokov (serafim.nenarokov@gmail.com)

import os
import argparse
import re
from pathlib import Path
import pathlib
import json

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import FastaIO
from tqdm import tqdm

from lib.nextaligner import Nextaligner
from lib.clusterizer import Clusterizer
from lib.mutation_detector import MutationDetector
from lib.utils import *

ROOT_PATH = str(Path(os.path.dirname(os.path.realpath(__file__))))
TMP_PATH = os.path.join(ROOT_PATH, "tmp")
REFERENCE_PATH = os.path.join(ROOT_PATH, "data", "sars_cov_reference_prot.fasta")
NEXTALIGN_BIN_PATH = os.getenv('NEXTALIGN_PATH') or 'nextalign'
REPORT_JSON_NAME = 'report.json'
REPORT_CSV_NAME = 'report.csv'

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

    out_path_help_str = """Path to the output folder
    (will be created if doesn't exist)"""

    parser.add_argument("-i", "--input",
                        required=True,
                        help="Input fasta file (.fasta)")
    parser.add_argument('-o', '--output_folder',
                        required=True,
                        help=out_path_help_str)
    return parser.parse_args()

def get_reference_path(gene_name):
    file_name = f"{gene_name}.fasta"
    return os.path.join(ROOT_PATH, 'data', 'reference', file_name)

def main():
    args = parse_arguments()

    assure_file_exists(args.input)


    # Align sequences with Nextaligner

    aligner = Nextaligner()
    aligner.perform(args.input, args.output_folder)


    # Cluster sequences
    print("Clustering")

    gene_clusters = {}

    nextalign_files = aligner.get_ouput_files()

    for gene_name, path in tqdm(nextalign_files.items()):
        clusterizer = Clusterizer(path)
        clusterizer.perform()
        gene_clusters[gene_name] = clusterizer.get_cluster_info_list()


    # Detect mutations
    print("Detecting mutations")

    for gene_name, clusters in gene_clusters.items():
        ref_path = get_reference_path(gene_name)
        ref_record = SeqIO.read(ref_path, 'fasta')
        ref_seq = ref_record.seq

        cluster_iter = tqdm(clusters)

        for cl in cluster_iter:
            cluster_iter.set_description(gene_name)
            query_sec = cl['seq']
            mut_detector = MutationDetector(ref_seq, query_sec, gene_name)
            cl['mutations'] = mut_detector.get_mutations()
            cl['aligned_seq'] = mut_detector.aligned_sequence


    # Preparing the list of aligned sequences
    print("Preparing the list of aligned sequences")

    aligned_seqs = {}

    for gene_name, clusters in gene_clusters.items():
        aligned_seqs[gene_name] = []

        for cl in clusters:
            for header in cl['all_seq_ids']:
                record = SeqRecord(Seq(cl['aligned_seq']),
                                   id=header,
                                   description='')
                aligned_seqs[gene_name].append(record)

    for gene_name, seq_recs in aligned_seqs.items():
        path = os.path.join(args.output_folder, f"gene.{gene_name}.fasta")
        with open(path, 'w') as out_f:
            fasta_out = FastaIO.FastaWriter(out_f, wrap=None)
            fasta_out.write_file(seq_recs)


    # collect per-mutation statistics
    print("Preparing statistics")

    stats = {}

    for gene_name, clusters in gene_clusters.items():
        cluster_iter = tqdm(clusters)

        for cluster in cluster_iter:
            cluster_iter.set_description(gene_name)

            for mutation in cluster['mutations']:
                mut_code = mutation['code']
                mut_seqs = cluster['all_seq_ids']

                if mut_code not in stats:
                    stats[mut_code] = {'seq_ids': []}

                stats[mut_code]['seq_ids'] += mut_seqs

    for _mut_code, data in stats.items():
        uniq_seq_ids = list(set(data['seq_ids']))
        uniq_seq_ids.sort()
        data['seq_ids'] = uniq_seq_ids
        data['cnt'] = len(uniq_seq_ids)


    # Generate reports
    print("Generating reports")
    report_json_path = os.path.join(args.output_folder, REPORT_JSON_NAME)

    with open(report_json_path, 'w') as f:
        f.write(json.dumps(stats, indent=4, sort_keys=True))

    report_csv_path = os.path.join(args.output_folder, REPORT_CSV_NAME)
    with open(report_csv_path, 'w') as f:
        for mut_id, data in stats.items():
            ids = ';'.join(data["seq_ids"])
            s = f'{mut_id},{ids}\n'
            f.write(s)

    print("Done")


if __name__ == '__main__':
    main()
