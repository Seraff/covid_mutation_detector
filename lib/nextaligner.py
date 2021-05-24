import os
import re
import subprocess
from pathlib import Path

from lib.utils import *

ROOT_PATH = str(Path(os.path.dirname(os.path.realpath(__file__))).parent)
TMP_PATH = os.path.join(ROOT_PATH, "tmp")
GENEMAP_PATH = os.path.join(ROOT_PATH, "data", "sars_cov_genemap.gff")
REFERENCE_PATH = os.path.join(ROOT_PATH, "data", "sars_cov_reference.fasta")

class Nextaligner:
    def __init__(self, binary='nextalign'):
        self.binary = binary

        self.output_path = None
        self.performed = False
        self.genes = None

    def perform(self, input_path, output_path, genes=None):
        self.output_path = output_path
        self.genes = genes or GENE_LIST

        create_dir_if_not_exists(output_path)

        genes = ','.join(self.genes)

        command = f'''{self.binary} --sequences={input_path}\
                                    --reference={REFERENCE_PATH}\
                                    --genemap={GENEMAP_PATH}\
                                    --genes={genes}\
                                    --output-dir={output_path}\
                                    --output-basename=nextalign'''

        command = re.sub('\s+', ' ', command)

        ex_code = subprocess.call(command, shell=True)

        if (ex_code != 0):
            print(f"{command} \n returned status {ex_code}")
            exit(1)

        self.performed = True

    def get_ouput_files(self):
        '''
        Returns the dictionary { 'gene_name': 'path/to/file.fasta'}
        '''

        result = {}

        for gene in self.genes:
            filename = f'nextalign.gene.{gene}.fasta'
            path = os.path.join(self.output_path, filename)

            assure_file_exists(path)
            result[gene] = path

        return result


if __name__ == '__main__':
    aligner = Nextaligner()
    fasta_path = os.path.join(ROOT_PATH, "tmp", "0509", "20210507_short.fasta")
    output_path = os.path.join(ROOT_PATH, "tmp", "nextaligner_test")
    aligner.perform(fasta_path, output_path)
