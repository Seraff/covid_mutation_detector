import os
import shutil
from pathlib import Path

from utils import *

ROOT_PATH = str(Path(os.path.dirname(os.path.realpath(__file__))).parent)
TMP_PATH = os.path.join(ROOT_PATH, "tmp")
REFERENCE_PATH = os.path.join(ROOT_PATH, "data", "sars_cov_genemap.gff")
GENEMAP_PATH = os.path.join(ROOT_PATH, "data", "sars_cov_reference.fasta")

class Nextaligner:
    def __init__(self, binary='nexalign'):
        self.tmp_path = os.path.join(TMP_PATH, record.id)
        self.S_path = os.path.join(self.tmp_path, 'nextalign.gene.S.fasta')
        self.M_path = os.path.join(self.tmp_path, 'nextalign.gene.M.fasta')
        self.records = None

    def perform(self, path):
        create_dir_if_not_exists(self.tmp_path)

        command = f'''nextalign --sequences={path}
                                --reference={REFERENCE_PATH}
                                --genemap={GENEMAP_PATH}
                                --genes=M,S
                                --output-dir={self.tmp_path}
                                --output-basename=nextalign'''

        ex_code = subprocess.call(command, shell=True)

        if (ex_code != 0):
            print(f"{command} \n returned status {ex_code}")
            exit(1)

    def is_performed(self):
        return self.record != None

    def clean_tmp(self):
        if self.is_performed():
            shutil.rmtree(self.tmp_path)
