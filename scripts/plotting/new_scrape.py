#!/usr/bin/env python3
# @Date    : 2022-10-12
# @Author  : Serafim Nenarokov (serafim.nenarokov@gmail.com)

import os
import sys
import pathlib
from pprint import pprint
import argparse
from tabnanny import verbose

from click import argument

ROOT_PATH = str(pathlib.Path(os.path.dirname(os.path.realpath(__file__))).parent.parent)
sys.path.append(ROOT_PATH)

from lib.outbreak_api import OutbreakApi


def parse_arguments():
    usage = "./new_scrape.py"
    description = """

    """
    formatter = argparse.RawTextHelpFormatter

    parser = argparse.ArgumentParser(usage,
                                     description=description,
                                     formatter_class=formatter)

    parser.add_argument('-o', '--output',
                        required=True,
                        help="Output folder")

    parser.add_argument('--verbose',
                        action='store_true',
                        help="Detailed output")

    parser.add_argument('lineages', nargs='+', help="List of lineages")

    return parser.parse_args()


def make_alias_dict(raw_list):
    result = {}

    for variant in raw_list:
        who_name = variant['label']

        if who_name:
            result[who_name] = variant['char_muts_parent_query']

    return result


def main():
    arguments = parse_arguments()
    lineages = arguments.lineages
    output_path = arguments.output

    if len(lineages) == 0:
        exit(0)

    if not os.path.exists(output_path):
        os.mkdir(output_path)

    api = OutbreakApi(timeout=64, quiet=not arguments.verbose)

    curated_list = api.get_curated_lineages()
    aliases = make_alias_dict(curated_list)


    for lineage in lineages:

        query_lineage = lineage

        if lineage in aliases:
            query_lineage = aliases[lineage]

        mutations = api.mutations_by_lineage(query_lineage)

        if mutations == None:
            exit(-1)

        if len(mutations) > 0:
            with open(os.path.join(output_path, f"{lineage}.tsv"), 'w') as out_f:
                for data in mutations[query_lineage]:
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
