#!/usr/bin/env python3
# @Date    : 2021-04-14
# @Author  : Serafim Nenarokov (serafim.nenarokov@gmail.com)

import os
import glob
from pathlib import Path

from Bio import SeqIO


ROOT_PATH = str(Path(os.path.dirname(os.path.realpath(__file__))).parent)
REFERENCE_DIR_PATH = os.path.join(ROOT_PATH, 'data', 'reference')


class ReferenceGenes:
    def __init__(self):
        self.genes = {}
        self.loaded = False

    def load(self):
        paths = glob.glob(os.path.join(REFERENCE_DIR_PATH, '*'))

        for path in paths:
            name = os.path.basename(path)
            name = name.split('.')[0]

            rec = SeqIO.read(path, 'fasta')
            self.genes[name] = rec.seq

        self.loaded = True


    def get_aa(self, gene_name):
        if not self.loaded:
            print("ERROR: reference genes are not loaded")
            exit(1)

        return self.genes[gene_name]
