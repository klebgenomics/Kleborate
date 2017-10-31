#!/usr/bin/env python3

"""
This script modifies the ARGannot_clustered80_r2.csv file (comes with SRST2) for use with
Kleborate. Specifically, it does two things:
  1) Add beta-lactamase description and class columns
  2) Removes some genes not suitable for inclusion in Kleborate

It takes two arguments:
  1) SRST2's copy of the ARGannot_clustered80_r2.csv file
  2) the table made by bla_info.py

Usage:
  ./make_argannot_csv.py path/to/srst2/data/ARGannot_clustered80_r2.csv bla_info_table > ARGannot_clustered80_r2.csv
"""


import sys

argannot_csv_filename = sys.argv[1]
bla_parse_output_filename = sys.argv[2]

descriptions = {}
bla_classes = {}

with open(bla_parse_output_filename, 'rt') as bla_parse_output_file:
    for line in bla_parse_output_file:
        allele, description, bla_class = line.strip().split('\t')
        descriptions[allele] = description
        bla_classes[allele] = bla_class

new_csv_lines = []
with open(argannot_csv_filename, 'rt') as argannot_csv_file:
    for line in argannot_csv_file:
        stripped_line = line.rstrip()

        if stripped_line.startswith('clusterid,queryID'):
            print(stripped_line + ',bla_description,bla_class')
            continue

        # Fix verbose gene name
        stripped_line = stripped_line.replace('NimB_Nitroimidazole_Gene', 'NimB')

        line_parts = stripped_line.split(',')
        res_class = line_parts[2]
        gene = line_parts[3]
        allele = line_parts[4]

        # Skip some genes not relevant for Kleborate.
        if gene == 'OqxB' or gene == 'OqxA' or gene == 'Far1' or gene == 'OptrA':
            continue

        # Fix Colistin class name to just three letters.
        if res_class == 'Colistin':
            stripped_line = stripped_line.replace(',Colistin,', ',Col,')

        if line_parts[2] == 'Bla' and allele in bla_classes:
            new_csv_lines.append(stripped_line + ',' + descriptions[allele] + ',' + bla_classes[allele])
        else:
            new_csv_lines.append(stripped_line + ',NA,NA')

# Sort so 'TEM-2' will come before 'TEM-10', etc.
new_csv_lines = sorted(new_csv_lines, key=lambda x: (int(x.split(',')[0]), int(x.split(',')[5])))

for line in new_csv_lines:
    print(line)

