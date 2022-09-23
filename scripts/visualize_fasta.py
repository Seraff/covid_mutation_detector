#!/usr/bin/env python3
# @Date    : 2021-04-14
# @Author  : Serafim Nenarokov (serafim.nenarokov@gmail.com)

import os
from pathlib import Path
import argparse

from Bio import SeqIO
from tqdm import tqdm
from PIL import Image, ImageDraw, ImageFont

from lib.utils import *

ROOT_PATH = str(Path(os.path.dirname(os.path.realpath(__file__))))
REFERENCE_PATH = os.path.join(ROOT_PATH, "data", "sars_cov_reference.fasta")
REF_SEQ = str(SeqIO.read(REFERENCE_PATH, 'fasta').seq)

def parse_arguments():
    usage = "./visualize_fasta.py"
    description = """
    Draws COVID fasta
    """
    description = re.sub(r"(?<=\w)\n(?=\w)", "\n", description)
    formatter = argparse.RawTextHelpFormatter

    parser = argparse.ArgumentParser(usage,
                                     description=description,
                                     formatter_class=formatter)

    parser.add_argument("-i", "--input",
                        required=True,
                        help="Input fasta file (.fasta)")
    parser.add_argument('-o', '--output',
                        required=True,
                        help="Path to output image")
    return parser.parse_args()

COLORS = ['#A1C9F4',
          '#FFB482',
          '#8DE5A1',
          '#FF9F9B',
          '#D0BBFF',
          '#DEBB9B',
          '#FAB0E4',
          '#faff78',
          '#FFFEA3',
          '#B9F2F0',
          '#A1C9F4',
          '#FFB482',
          '#8DE5A1',
          '#FF9F9B']

GENE_COORDS =  {"E": [26245, 26472],
                "M": [26523, 27191],
                "N": [28274, 29533],
                "ORF10": [29558, 29674],
                "ORF14": [28734, 28955],
                "ORF1a": [266, 13468],
                "ORF1b": [13469, 21555],
                "ORF3a": [25393, 26220],
                "ORF6": [27202, 27387],
                "ORF7a": [27394, 27759],
                "ORF7b": [27756, 27887],
                "ORF8": [27894, 28259],
                "ORF9b": [28284, 28577],
                "S": [21563, 25384]}

GLOBAL_MARGIN_TOP = 100
FONT = ImageFont.truetype("Keyboard.ttf", 26)

def get_ref_color(coord):
    result = {}
    i = 0
    for gene, coords in GENE_COORDS.items():
        if coord >= coords[0] and coord <= coords[1]:
            result['color'] = COLORS[i]
            result['gene'] = gene
            result['first_coord'] = coord == coords[0]
            return result
        i += 1

    result['color'] = (201, 201, 201)
    result['gene'] = None
    return result


def draw_nucleotide(draw, match, h_index, v_index):
    if v_index == 0:
        margin_top = -10
        bottom = 20
    else:
        margin_top = 10
        bottom = 20

    if v_index == 0:
        color_data = get_ref_color(h_index)
        color = color_data['color']
    else:
        if match == 'YES':
            color = (125, 255, 164)
        elif match == 'NO':
            color = (181, 12, 0)
        elif match == 'N':
            color = (201, 201, 201)

    x = 50+h_index
    y0 = 20*v_index+margin_top+GLOBAL_MARGIN_TOP
    y1 = 20*v_index+bottom+GLOBAL_MARGIN_TOP

    if v_index == 0:
        if color_data['gene'] and color_data['first_coord']:
            draw.text((x, y0-30), color_data['gene'], font=FONT, fill=(0, 0, 0))

    draw.line(((x, y0), (x, y1)), fill=color)

def draw(img, ref, seq, v_index):
    draw = ImageDraw.Draw(img)

    length = len(max(REF_SEQ, seq))
    for i in range(length):
        if i >= len(REF_SEQ):
            draw_nucleotide(draw, 'NO', i, v_index)
            continue

        elif i >= len(seq):
            break

        if seq[i] == REF_SEQ[i]:
            draw_nucleotide(draw, 'YES', i, v_index)
        elif seq[i] != REF_SEQ[i] and seq[i].upper() == 'N':
            draw_nucleotide(draw, 'N', i, v_index)
        else:
            draw_nucleotide(draw, 'NO', i, v_index)


def main():
    args = parse_arguments()

    coords = []
    max_width = 0

    seqs = []
    seqs.append(REF_SEQ)

    for rec in SeqIO.parse(args.input, 'fasta'):
        seqs.append(rec.seq)

    width = 30000 # 50px left and right
    height = len(seqs) * 30 + GLOBAL_MARGIN_TOP
    img = Image.new('RGB', (width, height), color=(255, 255, 255))

    i = 0
    for seq in tqdm(seqs):
        draw(img, REF_SEQ, seq, i)
        i += 1

    img.save(args.output)


if __name__ == '__main__':
    main()
