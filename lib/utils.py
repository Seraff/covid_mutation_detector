import pathlib
import re

AA_ALPHABET = {'Y', 'V', 'P', 'C', 'A', 'M', 'G',
               'H', 'R', 'D', 'N', 'K', 'Q', 'W',
               'S', 'I', 'L', 'F', 'E', 'T', "*"}

GENE_LIST = ("E", "M", "N", "ORF10",
             "ORF14", "ORF1a", "ORF1b",
             "ORF3a", "ORF6", "ORF7a",
             "ORF7b", "ORF8", "ORF9b", "S")

MUTATION_REGEXP = re.compile((f"^({'|'.join(GENE_LIST)}):"
                              f"[{''.join(AA_ALPHABET)}-]"
                              f"\d+"
                              f"[{''.join(AA_ALPHABET)}-]$"))

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
