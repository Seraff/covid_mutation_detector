#!/usr/bin/env python3

import os
import json
from pathlib import Path
import requests
from datetime import datetime

from tqdm import tqdm

from lib.reference_genes import ReferenceGenes
from lib.utils import assure_mutation_correct

VARIATS_LIST_URL = "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/lineage_notes.txt"
API_URL = "https://api.outbreak.info/genomics/"
MUT_INFO_ACTION = "lineage-mutations"

ROOT_PATH = str(Path(os.path.dirname(os.path.realpath(__file__))))
OUTPUT_PATH = os.path.join(ROOT_PATH, 'data', "known_variants.json")

def response_successful(response):
    if response.status_code != 200:
        print(f"ERROR: response failed: {response.url}")
        return False
    else:
        return True

def try_to_get(url, params={}):
    success = False
    response = None

    while not success:
        try:
            response = requests.get(url, params=params, timeout=5)
        except:
            print("Connection problem, reconnecting...")
            success = False
            continue

        if (response_successful(response)):
            success = True
        else:
            print("Trying again...")

    if response == None:
        raise Exception(f'Something wrong happened with url {url}')

    return response

def process_indel_mutation(mutation, reference, type='deletion'):
    mutations = []

    start = int(mut['codon_num'])
    end = int(mut['codon_end']) if mut['codon_end'] != 'None' else start

    for aa_i in range(start, end+1):
        aa = reference.get_aa(mut['gene'])

        if len(aa) <= aa_i:
            print(f'Problem in mutation {mutation}, out of bounds')
            continue

        if type == 'deletion':
            name = f"{mut['gene']}:{aa[aa_i-1]}{aa_i}-"
        else:
            name = f"{mut['gene']}:-{aa_i}{aa[aa_i-1]}"
        mutations.append(name)

    return mutations


if __name__ == "__main__":
    response = try_to_get(VARIATS_LIST_URL)

    variants = []

    for i, line in enumerate(response.text.split('\n')):
        if i == 0:
            continue

        splitted = line.strip().split('\t')

        if len(splitted) > 1:
            desc = splitted[1]
        else:
            desc = ''

        var_id = splitted[0].strip()

        if var_id == '' or var_id[0] == '*':
            continue

        data = { 'id': var_id,
                 'mutations': [],
                 'description': desc }
        variants.append(data)

    reference = ReferenceGenes()
    reference.load()

    for var in tqdm(variants):
        url = f"{API_URL}{MUT_INFO_ACTION}"
        payload = {'pangolin_lineage': var['id'], 'frequency': 0.75}
        response = try_to_get(url, payload)

        response_json = json.loads(response.text)

        for mut in response_json['results'][var['id']]:
            if mut['type'] == 'substitution':
                name = f"{mut['gene']}:"
                name += mut['ref_aa']
                name += str(mut['codon_num'])
                name += mut['alt_aa']
                var['mutations'].append(name)

            elif mut['type'] == 'deletion':
                mutations = process_indel_mutation(mut, reference)
                [var['mutations'].append(n) for n in mutations]

            elif mut['type'] == 'insertion':
                mutations = process_indel_mutation(mut, reference, 'insertion')
                [var['mutations'].append(n) for n in mutations]

            for m in var['mutations']:
                if not assure_mutation_correct(m):
                    print(f"{var['id']}")
                    print(f"Mutation {m} has incorrect format!")
                    exit(-1)

            var['mutations'].sort()

    # keeping variants with mutations only
    variants = [v for v in variants if len(v['mutations']) > 0]

    result = {}
    result['variants'] = variants

    now = datetime.now()
    current_time = now.strftime("%d/%m/%Y %H:%M:%S")

    result['timestamp'] = current_time

    with open(OUTPUT_PATH, 'w') as f:
        f.write(json.dumps(result, indent=4))

    print(f'The results were written to {OUTPUT_PATH}.')





