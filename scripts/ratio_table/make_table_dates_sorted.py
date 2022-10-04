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
    parser = argparse.ArgumentParser(
        usage, description=description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-j", "--json_path",
                        required=True,
                        help="Json file with mutation info (report)")
    parser.add_argument("-m", "--metadata_path",
                        required=True,
                        help="Metadata file")
    parser.add_argument("-d", "--date_column",
                        required=True,
                        type=int,
                        help="Sampling date column in metadata file")
    parser.add_argument("-o", "--output_path",
                        required=True,
                        help="Output table csv file")
    return parser.parse_args()


def date_to_week(date):
    if len(date.split('-')[0]) != 4:
        try:
            d, m, y = date.split('-')
            if int(y) <= 2019:
                week_number = 'NA'
            elif y == '2022' and datetime.date(int(y), int(m), int(d)).isocalendar()[1] == 52:
                week_number = '52-2021'
            else:
                week_number = datetime.date(int(y), int(m), int(d)).isocalendar()[1]
                if week_number in [1,2,3,4,5,6,7,8,9]:
                    week_number = '0%s-%s' % (week_number, y)
                else:
                    week_number = '%s-%s' % (week_number, y)
        except ValueError:
            print(date)
            week_number = 'NA'
    else:
        try:
            y, m, d = date.split('-')
            if int(y) <= 2019:
                week_number = 'NA'
            elif y == '2022' and datetime.date(int(y), int(m), int(d)).isocalendar()[1] == 52:
                week_number = '52-2021'
            else:
                week_number = datetime.date(int(y), int(m), int(d)).isocalendar()[1]
                if week_number in [1,2,3,4,5,6,7,8,9]:
                    week_number = '0%s-%s' % (week_number, y)
                else:
                    week_number = '%s-%s' % (week_number, y)
        except ValueError:
            print(date)
            week_number = 'NA'

    return(week_number)


if __name__ == '__main__':
    options = parse_options()
    infile = open(options.metadata_path)
    lines = infile.readlines()
    new_lines = []
    added_seqs = []

    for line in lines:
        if line.split('\t')[0] not in added_seqs:
            added_seqs.append(line.split('\t')[0])
            new_lines.append(line)

    lines = new_lines
    infile.close()
    SeqId_Date = {}
    Date_Seqs = {}

    for line in lines[1:]:
        Date = line.split('\t')[options.date_column - 1]
        SeqId_Date[line.split('\t')[0]] = date_to_week(Date)
        if date_to_week(Date) in Date_Seqs:
            Date_Seqs[date_to_week(Date)].append(line.split('\t')[0])
        else:
            Date_Seqs[date_to_week(Date)] = [line.split('\t')[0]]

    with open(options.json_path) as f:
        guys = json.loads(f.read())
        new_d1 = {}
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
                SeqId_Date[seq] = 'NA'
                Date_Seqs['NA'].append(seq)

    out = open(options.output_path, 'w')

    for date in sorted(Date_Seqs):
        out.write(',%s_%s' % (date, len(Date_Seqs[date])))

    out.write(',pocet')
    out.write('\n')

    # dates_d - number of sequences per day from the report.json
    # Date_Seqs - sequences per day { 'day': [seq, ...], ... }

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
                ratio = float(dates_d[date])/float(len(set(Date_Seqs[date])))
                out.write('%s, ' % (ratio))
            else:
                out.write('0,')

        out.write('%s, ' % (len(new_d[mutation])))
        out.write('\n')

    out.close()

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

        if c >= 1:
            out.write(line)
            Mutations_recorded.append(line.split(',')[0])

    d = {} # { 'date': ['seq_id', ...], ... }
    d2 = {} # { 'date': ratio, ... }

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

