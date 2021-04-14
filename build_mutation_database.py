#!/usr/bin/env python3
# @Date    : 2021-04-12
# @Author  : Serafim Nenarokov (serafim.nenarokov@gmail.com)

import os
from pathlib import Path

ROOT_PATH = str(Path(os.path.dirname(os.path.realpath(__file__))))
TMP_PATH = os.path.join(ROOT_PATH, "tmp")
REFERENCE_PATH = os.path.join(ROOT_PATH, "data", "sars_cov_genemap.gff")
GENEMAP_PATH = os.path.join(ROOT_PATH, "data", "sars_cov_reference.fasta")

if __name__ == '__main__':
    pass
