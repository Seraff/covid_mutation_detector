#!/usr/bin/env python3
# @Date    : 2021-05-10
# @Author  : Serafim Nenarokov (serafim.nenarokov@gmail.com)

import os
import re
import argparse

import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt

def parse_arguments():
    usage = "./generate_heatmap.py"
    description = """
    Generates mutations heatmap from input table file.
    """
    description = re.sub(r"(?<=\w)\n(?=\w)", "\n", description)
    formatter = argparse.RawTextHelpFormatter

    parser = argparse.ArgumentParser(usage,
                                     description=description,
                                     formatter_class=formatter)


    parser.add_argument("-i", "--input",
                        required=True,
                        help="Input table file (CSV format)")
    parser.add_argument('-o', '--output',
                        required=True,
                        help="Output image file (PNG format)")
    return parser.parse_args()

def main():
    args = parse_arguments()

    data = pd.read_csv(args.input, index_col=0)
    orf_data = data[data.index.str.startswith('ORF1a:')]

    cmap = sns.color_palette("Blues", as_cmap=True)
    fig, ax = plt.subplots(figsize=(10, 40))
    plt.tight_layout(pad=15)


    heatmap = sns.heatmap(data=data,
                       cmap=cmap,
                       ax=ax,
                       center=0.5,
                       linewidths=1)

    figure = heatmap.figure
    figure.savefig(args.output)

if __name__ == '__main__':
    main()
