#!/usr/bin/env python3
"""
Copyright 2018 Kat Holt
Copyright 2018 Ryan Wick (rrwick@gmail.com)
https://github.com/katholt/Kleborate/

This file is part of Kleborate. Kleborate is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Kleborate is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Kleborate. If
not, see <http://www.gnu.org/licenses/>.
"""

import sys
import argparse
import collections
import pandas as pd
from Bio import Phylo


def get_arguments():
    parser = argparse.ArgumentParser(description='A script for converting Kleborate output into a'
                                                 'format compatible with Microreact')
    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument('--kleborate_in', type=str, required=True,
                               help='Kleborate tab-delimited results file')
    required_args.add_argument('--tree_in', type=str, required=True,
                               help='Phylogenetic tree')
    required_args.add_argument('--csv_out', type=str, required=True,
                               help='Kleborate results in Microreact format')
    required_args.add_argument('--tree_out', type=str, required=True,
                               help='Tree in Microreact format')

    return parser.parse_args()


def main():
    args = get_arguments()
    name_subs = name_substitution(args.kleborate_in)
    check_for_unique_names(name_subs)
    save_tree_with_new_names(args.tree_in, args.tree_out, name_subs)
    colours = get_colours(args.kleborate_in)
    autocolour_columns = get_autocolour_columns(args.kleborate_in)

    csv_lines = []
    with open(args.kleborate_in, 'rt') as kleborate_results:
        header, columns = None, {}
        for line in kleborate_results:
            line = line.rstrip('\n')
            if header is None:
                header = get_header(line, columns, autocolour_columns)
                line_parts = header
            else:
                line_parts = get_data(line, columns, colours, name_subs)
            csv_lines.append((','.join(line_parts)))

    with open(args.csv_out, 'wt') as output_csv:
        for line in csv_lines:
            output_csv.write(line)
            output_csv.write('\n')

    print()


def get_colours(kleborate_in):
    table = pd.read_table(kleborate_in)
    try:
        max_classes = max(table['num_resistance_classes'])
    except (KeyError, ValueError):
        max_classes = 10
    try:
        max_genes = max(table['num_resistance_genes'])
    except (KeyError, ValueError):
        max_genes = 20
    colours = {'vir_score': colour_range("#CCCCCC", "#1414FF", 4),
               'res_score': colour_range("#CCCCCC", "#FF1414", 3),
               'res_classes': colour_range("#CCCCCC", "#FF1414", max_classes),
               'res_genes': colour_range("#CCCCCC", "#FF1414", max_genes)}
    return colours


def get_autocolour_columns(kleborate_in):
    autocolour_columns = []
    table = pd.read_table(kleborate_in)
    for col_name in ['species', 'ST', 'Yersiniabactin', 'YbST', 'Colibactin', 'CbST',
                     'Aerobactin', 'AbST', 'Salmochelin', 'SmST', 'hypermucoidy', 'wzi',
                     'K_locus', 'O_locus']:
        try:
            if len(set(table[col_name])) > 1:
                autocolour_columns.append(col_name)
        except KeyError:
            pass
    print()
    print('Using "__autocolour" on the following columns:')
    print('   ', ', '.join(autocolour_columns))
    return set(autocolour_columns)


def get_header(line, columns, autocolour_columns):
    header = line.split('\t')
    header[0] = 'id'  # Change 'strain' to 'id' for Microreact.

    # Add new colour columns.
    columns['vir_score'] = find_column_index(header, 'virulence_score')
    columns['res_score'] = find_column_index(header, 'resistance_score')
    columns['res_classes'] = find_column_index(header, 'num_resistance_classes')
    columns['res_genes'] = find_column_index(header, 'num_resistance_genes')

    assert (columns['vir_score'] < columns['res_score']
            < columns['res_classes'] < columns['res_genes'])

    header.insert(columns['res_genes'] + 1, 'num_resistance_genes__colour')
    header.insert(columns['res_classes'] + 1, 'num_resistance_classes__colour')
    header.insert(columns['res_score'] + 1, 'resistance_score__colour')
    header.insert(columns['vir_score'] + 1, 'virulence_score__colour')

    for autocolour_column in autocolour_columns:
        i = find_column_index(header, autocolour_column)
        header[i] = autocolour_column + '__autocolour'

    return header


def get_data(line, columns, colours, name_subs):
    line = line.replace(',', ';')
    line_parts = line.split('\t')

    line_parts[0] = name_subs[line_parts[0]]

    virulence_score = int(line_parts[columns['vir_score']])
    resistance_score = int(line_parts[columns['res_score']])
    num_resistance_classes = int(line_parts[columns['res_classes']])
    num_resistance_genes = int(line_parts[columns['res_genes']])

    virulence_score_colour = colours['vir_score'][virulence_score]
    resistance_score_colour = colours['res_score'][resistance_score]
    resistance_class_colour = colours['res_classes'][num_resistance_classes]
    resistance_gene_colour = colours['res_genes'][num_resistance_genes]

    line_parts.insert(columns['vir_score'] + 1, virulence_score_colour)
    line_parts.insert(columns['res_score'] + 1, resistance_score_colour)
    line_parts.insert(columns['res_classes'] + 1, resistance_class_colour)
    line_parts.insert(columns['res_genes'] + 1, resistance_gene_colour)

    # if '.' in line_parts[0]:
    #     line_parts[0] = line_parts[0].replace('.', '_')

    return line_parts


def name_substitution(kleborate_in):
    name_subs = {}
    with open(kleborate_in, 'rt') as kleborate_results:
        header = None
        for line in kleborate_results:
            if header is None:
                header = line.split('\t')
                if header[0] != 'strain':
                    sys.exit('Error: first column is not "strain" - is this Kleborate output?')
            else:
                line_parts = line.split('\t')
                if len(line_parts) != len(header):
                    sys.exit('Error: inconsistent number of columns')
                old_name = line_parts[0]
                if old_name in name_subs:
                    sys.exit('Error: duplicate sample ID: ' + old_name)
                new_name = old_name.replace('.', '_')
                new_name = new_name.replace(',', '_')
                new_name = new_name.replace("'", '_')
                new_name = new_name.replace('"', '_')
                name_subs[old_name] = new_name
    return name_subs


def check_for_unique_names(name_subs):
    names = list(name_subs.values())
    duplicate_names = [item for item, count in collections.Counter(names).items() if count > 1]
    if duplicate_names:
        sys.exit('Error: duplicate sample IDs: ' + ', '.join(duplicate_names))


def save_tree_with_new_names(tree_in, tree_out, name_subs):
    tree_format = None
    for try_tree_format in ['newick', 'nexus', 'nexml', 'phyloxml', 'cdao']:
        try:
            Phylo.read(tree_in, try_tree_format)
            tree_format = try_tree_format
            break
        except ValueError:
            pass
    if tree_format is None:
        sys.exit('Error: could not read input tree')

    tree = Phylo.read(tree_in, tree_format)
    for node in tree.get_terminals():
        name = str(node.name)
        try:
            node.name = name_subs[name]
        except IndexError:
            sys.exit('Error: sample name in tree not in Kleborate data: ' + name)
    Phylo.write(tree, tree_out, 'newick')


def scale_num(start, end, progress):
    return int(round(start * (1.0 - progress) + end * progress))


def colour_range(start, end, count):
    start, end = start.lower(), end.lower()
    if start.startswith('#'):
        start = start[1:]
    if end.startswith('#'):
        end = end[1:]
    start_r, start_g, start_b = int(start[0:2], 16), int(start[2:4], 16), int(start[4:6], 16)
    end_r, end_g, end_b = int(end[0:2], 16), int(end[2:4], 16), int(end[4:6], 16)
    colours = []
    for i in range(count):
        progress = i / (count - 1)
        r, g, b = scale_num(start_r, end_r, progress), scale_num(start_g, end_g, progress), \
            scale_num(start_b, end_b, progress)
        hex_colour = '"#' + ('0x%X' % r)[2:] + ('0x%X' % g)[2:] + ('0x%X' % b)[2:] + '"'
        colours.append(hex_colour)
    return colours


def find_column_index(header, col_name):
    try:
        return header.index(col_name)
    except ValueError:
        sys.exit('Error: could not find ' + col_name + ' column in Kleborate')


if __name__ == '__main__':
    main()
