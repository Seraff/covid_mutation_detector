#!/usr/bin/env python3
# @Author: Martin Kolisko

import os
import sys
import glob
import json
import pickle
import argparse
import datetime


def parse_options():
    usage = "./make_table_columns_sorted.py"
    description = "Takes JSON file with per genome AA mutations"
    parser = argparse.ArgumentParser(usage, description=description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-j", "--json_path",
                        required=True,
                        help="Json file with mutation info (report)")
    parser.add_argument("-m", "--metadata_path",
                        required=True,
                        help="Metadata file")
    parser.add_argument("-d", "--data_column",
                        required=True,
                        type=int,
                        help="Sampling data column in metadata file")
    parser.add_argument("-o", "--output_path",
                        required=True,
                        help="Output table csv file")
    return parser.parse_args()


if __name__ == '__main__':
    Mutations_To_keep = []
    options = parse_options()
    infile = open(options.metadata_path)
    lines = infile.readlines()
    new_lines = []
    added_seqs = []

    for line in lines:
        if line.split('\t')[0] in added_seqs:
            pass
        else:
            added_seqs.append(line.split('\t')[0])
            new_lines.append(line)

    lines = new_lines
    infile.close()
    SeqId_Date = {}
    Date_Seqs = {}

    for line in lines[1:]:
        Date = line.split('\t')[options.data_column - 1]
        SeqId_Date[line.split('\t')[0]] = Date
        if Date in Date_Seqs:
            Date_Seqs[Date].append(line.split('\t')[0])
        else:
            Date_Seqs[Date] = [line.split('\t')[0]]

    print(len(Date_Seqs))
    print(len(SeqId_Date))

    with open(options.json_path) as f:
        guys = json.loads(f.read())
        new_d1	 = {}
        for mut_id, mut_data in guys.items():
            if len(set(mut_data['seq_ids'])) >= 1:
                new_d1[mut_id] = mut_data['seq_ids']

        new_d = {}

        for mut in new_d1:
            l = []
            for rec in new_d1[mut]:
                newrec = rec.split('|')[0]
                newrec = newrec.replace('Czech_Republic','CzechRepublic')
                l.append(newrec)
            new_d[mut] = l

    missing_seqs_inmeta = []

    for mut in new_d:
        for seq in new_d[mut]:
            if seq not in SeqId_Date:
                missing_seqs_inmeta.append(seq)
                SeqId_Date[seq] = ''
                Date_Seqs[''].append(seq)

    print(len(Date_Seqs))
    print(len(SeqId_Date))
    print(len(set(missing_seqs_inmeta)))

    out = open(options.output_path, 'w')

    for date in sorted(Date_Seqs):
        out.write(',%s_%s' % (date, len(Date_Seqs[date])))

    out.write(',pocet_0')
    out.write('\n')

    for mutation in new_d:
        dates_d = {}
        for seq in new_d[mutation]:
            if seq in SeqId_Date:
                if SeqId_Date[seq] in dates_d:
                    dates_d[SeqId_Date[seq]] = dates_d[SeqId_Date[seq]] + 1
                else:
                    dates_d[SeqId_Date[seq]] = 1
            else:
                print('KURVA %s is still not there' % (seq))

        out.write('%s, ' % (mutation))

        for date in sorted(Date_Seqs):
            if date in dates_d:
                if int(dates_d[date]) >= 3:
                    Mutations_To_keep.append(mutation)
                out.write('%s, ' % (float(dates_d[date])/float(len(set(Date_Seqs[date])))))
            else:
                out.write('0,')
            #	print(date)
        out.write('%s, ' % (len(new_d[mutation])))
        out.write('\n')
    out.close()

    print(len(set(Mutations_To_keep)))

    Mutations_recorded = []

    with open(options.output_path) as infile:
        lines = infile.readlines()

    out = open('%s_reduced.csv' % (options.output_path), 'w')
    out.write(lines[0])

    for line in lines[1:]:
        regions = lines[0].split(',')[1:]
        values = line.split(',')[1:-1]
        c = 0

        for index, value in enumerate(values[:-1]):
            print(index)
            print(regions[index])
            if float(value.strip()) >= 0.1 and (float(value)*int(regions[index].split('_')[1])) >= 3:
                c = c + 1
            else:
                pass

        if c >= 1:
            out.write(line)
            Mutations_recorded.append(line.split(',')[0])

    d = {}
    d2 = {}
    for mut in new_d:
        if mut not in Mutations_recorded:
            for seq in new_d[mut]:
                if SeqId_Date[seq] in d:
                    d[SeqId_Date[seq]].append(seq)
                else:
                    d[SeqId_Date[seq]] = [seq]

            for seq in new_d[mut]:
                d2[SeqId_Date[seq]] = float(len(d[SeqId_Date[seq]]))/float(len(set(Date_Seqs[SeqId_Date[seq]])))

    out.write('other, ')
    for date in sorted(Date_Seqs):
        out.write('%s, ' % (d2[date]))

    out.close()
