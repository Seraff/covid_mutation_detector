#!/usr/bin/env python3
# @Date    : 2022-10-12
# @Author  : Serafim Nenarokov (serafim.nenarokov@gmail.com)

import re
import os
import sys
import argparse
import pathlib
import requests

ROOT_PATH = str(pathlib.Path(os.path.dirname(os.path.realpath(__file__))).parent)

sys.path.append(ROOT_PATH)
from lib.outbreak_api import OutbreakApi

API_URL = "https://api.outbreak.info/genomics"
API_METHOD = "prevalence-by-location"

def parse_arguments():
    usage = "./get_mutation_prevalence.py"
    description = """
    Calls outbreak.info API and return the prevalence of given mutation.
    """
    description = re.sub(r"(?<=\w)\n(?=\w)", "\n", description)
    formatter = argparse.RawTextHelpFormatter

    parser = argparse.ArgumentParser(usage,
                                     description=description,
                                     formatter_class=formatter)

    parser.add_argument("-m", "--mutation",
                        required=True,
                        help="Mutation (i.e. `S:R78K`)")
    return parser.parse_args()


def get_prevalence(mut_id):
    params = {
        'cumulative': 'true',
        'mutations': mut_id
    }

    api = OutbreakApi()
    result = api.get_data('prevalence-by-location', params=params)
    return result

def main():
    arguments = parse_arguments()
    data = get_prevalence(arguments.mutation)

    if data:
        print(data['results'][arguments.mutation]['lineage_count'])


if __name__ == "__main__":
    main()
