#!/usr/bin/env python3
# @Date    : 2022-10-12
# @Author  : Serafim Nenarokov (serafim.nenarokov@gmail.com)

import os
import sys
import pathlib
from pprint import pprint

ROOT_PATH = str(pathlib.Path(os.path.dirname(os.path.realpath(__file__))).parent.parent)
sys.path.append(ROOT_PATH)

from lib.outbreak_api import OutbreakApi

OUTPUT_PATH = "./VariantSignaturesNew"

def main():
    lineages = sys.argv[1:]

    if len(lineages) == 0:
        exit(0)

    if not os.path.exists(OUTPUT_PATH):
        os.mkdir(OUTPUT_PATH)

    api = OutbreakApi()

    for lineage in lineages:
        mutations = api.mutations_by_lineage(lineage)

        if mutations == None:
            exit(-1)

        with open(os.path.join(OUTPUT_PATH, f"{lineage}.tsv"), 'w') as out_f:
            for data in mutations[lineage]:
                cnt = f"{data['mutation_count']} / {data['lineage_count']}"
                mutation = data['mutation'].upper()

                replaces = (('/', '-'), ('*', 'STOP'))
                for rule in replaces:
                    mutation = mutation.replace(rule[0], rule[1])

                line = [mutation, data['prevalence'], cnt]
                line = [str(i) for i in line]

                out_f.write("\t".join(line) + '\n')


if __name__ == '__main__':
    main()
