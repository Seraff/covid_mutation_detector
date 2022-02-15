import os
import glob
import pathlib
import re
from pathlib import Path

from Bio import SeqIO

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
REF_AA_PATH = os.path.join(ROOT_PATH, 'data', 'reference')

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
    return re.match(MUTATION_REGEXP, name) != None


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

    result = {}
    for rec in SeqIO.parse(path, 'fasta'):
        if rec.id not in result:
            result[rec.id] = {}

        for gene_name, coords in COVID_GENES.items():
            sequence = str(rec.seq[coords[0]-1:coords[1]])
            result[rec.id][gene_name] = sequence

    return result

def load_reference_genes():
    ref_aa_seqs = load_nextclade_aa_genes(REF_AA_PATH)
    ref_aa_seqs = ref_aa_seqs[list(ref_aa_seqs.keys())[0]]
    return ref_aa_seqs
