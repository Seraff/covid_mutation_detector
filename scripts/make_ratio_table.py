#!/usr/bin/env python3
# @Author: Martin Kolisko, Serafim Nenarokov

import json
import argparse
import datetime
from tqdm import tqdm


def parse_options():
    usage = "./make_ratio_table.py"
    description = "Takes JSON file with per genome AA mutations"
    parser = argparse.ArgumentParser(
        usage, description=description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-r", "--report_path",
                        required=True,
                        help="Report file with mutation data in JSON format")
    parser.add_argument("-m", "--metadata_path",
                        required=True,
                        help="Metadata file")
    parser.add_argument("-d", "--data_column",
                        required=True,
                        type=int,
                        help="Sampling data column in metadata file (indexing starts from 1)")
    parser.add_argument("--to_week_format",
                        action='store_true',
                        help="Should the column be converted to the week format (i.e. `05-2022`)")
    parser.add_argument("-o", "--output_path",
                        required=True,
                        help="Output CSV file")
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


def format_report_seq_id(seq_id):
    result = seq_id.split('|')[0]
    result = result.replace('Czech_Republic', 'CzechRepublic')
    return result

if __name__ == '__main__':
    arguments = parse_options()
    # Reading metadata

    with open(arguments.metadata_path) as f:
        metadata_lines = f.readlines()

    new_lines = []
    added_seqs = []

    for line in metadata_lines:
        if line.split('\t')[0] not in added_seqs:
            added_seqs.append(line.split('\t')[0])
            new_lines.append(line)

    metadata_lines = new_lines
    dates_by_seq_ids = {}
    seq_ids_by_date = {}

    for line in metadata_lines[1:]:
        date = line.split('\t')[arguments.data_column - 1]

        if arguments.to_week_format:
            date = date_to_week(date)

        seq_id = line.split('\t')[0]

        dates_by_seq_ids[seq_id] = date

        if date in seq_ids_by_date:
            seq_ids_by_date[date].append(seq_id)
        else:
            seq_ids_by_date[date] = [seq_id]

    # Reading report

    with open(arguments.report_path) as f:
        report = json.load(f)

    for mut_id, data in report.items():
        new_seq_ids = []
        for seq_id in data['seq_ids']:
            new_seq_ids.append(format_report_seq_id(seq_id))
        data['seq_ids'] = new_seq_ids

        new_unknown_seq_ids = []
        for seq_id in data['unknown_seq_ids']:
            new_unknown_seq_ids.append(format_report_seq_id(seq_id))
        data['unknown_seq_ids'] = new_unknown_seq_ids

    for mut, data in report.items():
        seq_ids = set(data['seq_ids'])
        seq_ids = seq_ids.union(set(data['unknown_seq_ids']))

        for seq_id in seq_ids:
            if seq_id not in dates_by_seq_ids:
                dates_by_seq_ids[seq_id] = 'NA'
                if 'NA' not in seq_ids_by_date:
                    seq_ids_by_date['NA'] = []
                seq_ids_by_date['NA'].append(seq_id)

    # Making the main statistics

    out = open(arguments.output_path, 'w')

    for date in sorted(seq_ids_by_date):
        out.write(f",{ date }_{ len(seq_ids_by_date[date]) }")

    out.write(',pocet\n')

    for mut_id, data in tqdm(report.items(), desc="Processing the report"):
        # { date: { 'seq_ids': cnt, 'unknown_seq_ids': cnt }, ... }
        cnt_by_dates = {}

        all_seq_ids = set(data['seq_ids']).union(set(data['unknown_seq_ids']))

        for seq_id in all_seq_ids:
            if seq_id not in dates_by_seq_ids:
                raise ValueError(f'KURVA {seq_id} is still not in metadata')

            date = dates_by_seq_ids[seq_id]

            if date not in cnt_by_dates:
                cnt_by_dates[date] = {'seq_id_cnt': 0, 'unknown_seq_id_cnt': 0}

        for seq_id in data['seq_ids']:
            date = dates_by_seq_ids[seq_id]
            cnt_by_dates[date]['seq_id_cnt'] += 1

        for seq_id in data['unknown_seq_ids']:
            date = dates_by_seq_ids[seq_id]
            cnt_by_dates[date]['unknown_seq_id_cnt'] += 1

        out.write(f"{mut_id},")

        for date in sorted(seq_ids_by_date.keys()):
            if date in cnt_by_dates.keys():
                cnt = cnt_by_dates[date]['seq_id_cnt']
                unknown_cnt = cnt_by_dates[date]['unknown_seq_id_cnt']
                all_cnt = len(set(seq_ids_by_date[date]))

                if (all_cnt - unknown_cnt) == 0:
                    ratio = 0
                else:
                    ratio = float(cnt)/float(all_cnt - unknown_cnt)

                out.write(f"{ratio},")
            else:
                out.write('0,')

        out.write(f"{ len(report[mut_id]['seq_ids']) }\n")

    out.close()

    # Make the reduced statistics

    mutations_for_reduced = []

    with open(arguments.output_path) as infile:
        lines = infile.readlines()

    splitted_path = arguments.output_path.split('.')
    splitted_path.insert(-1, 'reduced')
    reduced_path = '.'.join(splitted_path)

    out = open(reduced_path, 'w')
    out.write(lines[0])

    for line in lines[1:]:
        regions = lines[0].split(',')[1:]
        values = line.split(',')[1:-1]

        cnt = 0
        for index, value in enumerate(values[:-1]):
            value = float(value.strip())
            cnt_in_region = int(regions[index].split('_')[1])

            if value >= 0.1 and (value*cnt_in_region) >= 3:
                cnt += 1

        if cnt >= 1:
            out.write(line)
            mutations_for_reduced.append(line.split(',')[0])


    reduced_seq_ids = {}
    for mut_id, data in report.items():
        if mut_id not in mutations_for_reduced:
            all_seq_ids = set(data['seq_ids']).union(set(data['unknown_seq_ids']))

            for seq_id in all_seq_ids:
                date = dates_by_seq_ids[seq_id]

                if date not in reduced_seq_ids:
                    reduced_seq_ids[date] = {'seq_ids': [], 'unknown_seq_ids': []}

            for seq_id in data['seq_ids']:
                date = dates_by_seq_ids[seq_id]
                reduced_seq_ids[date]['seq_ids'].append(seq_id)

            for seq_id in data['unknown_seq_ids']:
                date = dates_by_seq_ids[seq_id]
                reduced_seq_ids[date]['unknown_seq_ids'].append(seq_id)


    reduced_stats = {}
    for date in sorted(seq_ids_by_date.keys()):
        cnt = len(reduced_seq_ids[date]['seq_ids'])
        unknown_cnt = len(reduced_seq_ids[date]['unknown_seq_ids'])
        all_cnt = len(set(seq_ids_by_date[date]))

        reduced_stats[date] = float(cnt)/float(all_cnt)

    line = ','.join([str(reduced_stats[date]) for date in sorted(seq_ids_by_date)])
    out.write(f"other,{line}\n")

    out.close()

