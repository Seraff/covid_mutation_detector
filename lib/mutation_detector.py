#!/usr/bin/env python3
# @Date    : 2021-04-14
# @Author  : Serafim Nenarokov (serafim.nenarokov@gmail.com)

import json

from Bio.Seq import Seq
from Bio import pairwise2
from tqdm import tqdm

GAP_SYMBOL = '-'
UNKNOWN_SYMBOL = 'X'
STOP_SYMBOL = '*'


class MutationDetector:
    def __init__(self, ref_seq, query_seq, prot_name='S'):
        self.ref_seq = str(ref_seq)
        self.query_seq = str(query_seq)
        self.prot_name = prot_name


    def get_mutations(self):
        alignment = self._align()

        reference = alignment[0].strip()
        query = alignment[1].strip()

        mutations = []

        for i, ref_aa in enumerate(reference):
            query_aa = query[i].upper()

            if ref_aa == query_aa or query_aa == UNKNOWN_SYMBOL:
                continue

            # Insertion
            # ML-NA
            # MLPNA
            elif ref_aa == GAP_SYMBOL and query_aa != GAP_SYMBOL:
                code = f'{self.prot_name}:-{i+1}{query_aa}'
                mutations.append({'type': 'INSERTION',
                                  'from': i+1,
                                  'to': i+1,
                                  'code': code})

            # Deletion
            # MLPNA
            # ML-NA
            elif ref_aa != GAP_SYMBOL and query_aa == GAP_SYMBOL:
                code = f'{self.prot_name}:{ref_aa}{i+1}-'
                mutations.append({'type': 'DELETION',
                                  'from': i+1,
                                  'to': i+1,
                                  'code': code})

            # Substitution
            # MLPNA
            # MLINA
            elif ref_aa != query_aa:
                code = f'{self.prot_name}:{ref_aa}{i+1}{query_aa}'
                mutations.append({'type': 'SUBSTITUTION',
                                  'at': i+1,
                                  'before': ref_aa,
                                  'after': query_aa,
                                  'code': code})

        return mutations


    def _align(self):
        # Identical characters have score of 1, otherwise 0.
        # Same open and extend gap penalties for both sequences.
        # https://biopython.org/docs/1.75/api/Bio.pairwise2.html

        open_gap_score = -0.5
        extend_gap_score = -0.1
        result = pairwise2.align.globalxs(self.ref_seq,
                                          self.query_seq,
                                          open_gap_score,
                                          extend_gap_score)

        return [result[0].seqA, result[0].seqB]


if __name__ == '__main__':
    ref = 'MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT'

    clusters = None
    with open('tmp/final_output/clusters.json') as f:
        clusters = json.loads(f.read())

    for cluster in tqdm(clusters['clusters']):
       detector = MutationDetector(ref, cluster['seq'].replace('*', ''))
       cluster['mutations'] = detector.get_mutations()


    with open('tmp/final_output/clusters_mutations2.json', 'w') as f:
        f.write(json.dumps(clusters, indent=4))
