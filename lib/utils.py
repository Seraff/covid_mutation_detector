import os
import glob
import pathlib
import re
from pathlib import Path

from Bio import SeqIO
from tqdm import tqdm

AA_ALPHABET = {'Y', 'V', 'P', 'C', 'A', 'M', 'G',
               'H', 'R', 'D', 'N', 'K', 'Q', 'W',
               'S', 'I', 'L', 'F', 'E', 'T', "*"}

COVID_GENES = {
    "E": (26245, 26472),
    "M": (26523, 27191),
    "N": (28274, 29533),
    "ORF1a": (266, 13468),
    "ORF1b": (13468, 21555),
    "ORF3a": (25393, 26220),
    "ORF6": (27202, 27387),
    "ORF7a": (27394, 27759),
    "ORF7b": (27756, 27887),
    "ORF8": (27894, 28259),
    "ORF9b": (28284, 28577),
    "S": (21563, 25384)
}

MUTATION_REGEXP = re.compile((f"^({'|'.join(list(COVID_GENES.keys()))}):"
                              f"[{''.join(AA_ALPHABET)}-]"
                              f"\d+"
                              f"[{''.join(AA_ALPHABET)}-]$"))

ROOT_PATH = str(Path(os.path.dirname(os.path.realpath(__file__))).parent)
MUT_REPLACEMENT_RULES_PATH = os.path.join(
    ROOT_PATH, 'mut_replacement_rules.json')

REF_AA_PATH = os.path.join(ROOT_PATH, 'data', 'reference')
REF_NA_PATH = os.path.join(ROOT_PATH, 'data', 'sars_cov_reference.fasta')

MUTATION_REGEX = re.compile(
    '(?P<gene_name>\w+):(?P<source_aa>[a-zA-Z\-\*])(?P<index>\d+)(?P<target_aa>[a-zA-Z\-\*])')

def create_dir_if_not_exists(path):
    pathlib.Path(path).mkdir(parents=True, exist_ok=True)


def assure_file_exists(path):
    path = pathlib.Path(path)

    if not path.exists():
        print(f"Fasta file not found: {path}")
        exit(1)

    if not path.is_file():
        print(f"Fasta file is not a file: {path}")
        exit(1)


def assure_mutation_correct(name):
    return re.match(MUTATION_REGEX, name) != None


def load_nextclade_aa_genes(out_dir_path):
    file_paths = glob.glob(os.path.join(
        out_dir_path, '', f"*.fasta"))

    genes = list(COVID_GENES.keys())
    GENE_FILENAME_REGEX = re.compile(f".*({ '|'.join(genes) })\.fasta")

    file_paths = [p for p in file_paths if re.match(GENE_FILENAME_REGEX, p)]

    if len(file_paths) == 0:
        print("No nextalign files found")
        exit(-1)

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


def load_nextclade_na_alignments(out_dir_path):
    paths = glob.glob(os.path.join(out_dir_path, '', '*aligned.fasta'))

    if len(paths) == 0:
        print("No nextalign `*.aligned.fasta` file found")
        exit(-1)

    path = paths[0]

    return extract_genes_from_na_fasta(path)


def extract_genes_from_na_fasta(path):
    assure_file_exists(path)

    result = {}
    for rec in SeqIO.parse(path, 'fasta'):
        if rec.id not in result:
            result[rec.id] = {}

        for gene_name, coords in COVID_GENES.items():
            sequence = str(rec.seq[coords[0]-1:coords[1]])
            result[rec.id][gene_name] = sequence

    return result


def load_reference_aa_genes():
    ref_aa_seqs = load_nextclade_aa_genes(REF_AA_PATH)
    ref_aa_seqs = ref_aa_seqs[list(ref_aa_seqs.keys())[0]]
    return ref_aa_seqs


def load_reference_na_genes():
    ref_na_seqs = extract_genes_from_na_fasta(REF_NA_PATH)
    ref_na_seqs = ref_na_seqs[list(ref_na_seqs.keys())[0]]
    return ref_na_seqs


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


def by_mut_to_by_seq(by_mut_report):
    '''
    Receives report in COG format.

    Returns report in format:

    {'seq_id': {
      'gene_id': { site_id: {
                      code: "MUT_CODE",
                      type: "DEL/SUB"}, ...},
       ...},
    ...}
    '''

    by_seq_report = {}

    for mut_code, data in by_mut_report.items():
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

    return by_seq_report


def get_suspicious_mutations(by_mut_report):
    '''
    Receives report in COG format.

    Returns the dict of suspicious mutations in format:
    {'mut_code_1,mut_code_2,...': {'cnt': N, 'seq_ids': [...]}, ...}
    '''

    by_seq_report = by_mut_to_by_seq(by_mut_report)

    susp_mutations = {}

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

    return susp_mutations
