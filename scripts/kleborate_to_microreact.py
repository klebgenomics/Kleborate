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
    autocolour_columns = get_autocolour_columns(args.kleborate_in)

    csv_lines = []
    with open(args.kleborate_in, 'rt') as kleborate_results:
        original_header, new_header = None, None
        for line in kleborate_results:
            line = line.rstrip('\n')
            if original_header is None:
                original_header = line.split('\t')
                new_header = get_new_header(original_header, autocolour_columns)
                line_parts = new_header
            else:
                line_parts = get_data(line, name_subs, original_header, new_header)
            csv_lines.append((','.join(line_parts)))

    print()
    print('Writing Microreact table to: {}'.format(args.csv_out))
    with open(args.csv_out, 'wt') as output_csv:
        for line in csv_lines:
            output_csv.write(line)
            output_csv.write('\n')

    print()




def get_autocolour_columns(kleborate_in):
    autocolour_columns = []
    table = pd.read_table(kleborate_in)
    for col_name in ['species', 'ST', 'YbST', 'CbST', 'AbST', 'SmST', 'wzi', 'K_locus',
                     'O_locus']:
        try:
            if len(set(table[col_name])) > 1:
                autocolour_columns.append(col_name)
        except KeyError:
            pass
    print()
    print('Using "__autocolour" on the following columns:')
    print('   ', ', '.join(autocolour_columns))
    return set(autocolour_columns)


def get_new_header(original_header, autocolour_columns):
    original_header[0] = 'id'  # Change 'strain' to 'id' for Microreact.
    for autocolour_column in autocolour_columns:
        i = find_column_index(original_header, autocolour_column)
        original_header[i] = autocolour_column + '__autocolour'

    header = list(original_header)
    for col in ['virulence_score', 'resistance_score', 'num_resistance_classes',
                'num_resistance_genes', 'Yersiniabactin', 'Colibactin', 'Aerobactin',
                'Salmochelin', 'rmpA', 'rmpA2']:
        header.insert(find_column_index(header, col) + 1, col + '__colour')

    for res in ['AGly', 'Col', 'Fcyn', 'Flq', 'Gly', 'MLS', 'Ntmdz', 'Phe', 'Rif', 'Sul', 'Tet',
                'Tgc', 'Tmt', 'Bla', 'Bla_Carb', 'Bla_ESBL', 'Bla_ESBL_inhR', 'Bla_broad',
                'Bla_broad_inhR', 'Omp']:
        header.insert(find_column_index(header, res) + 1, res + '__colour')
        header.remove(res)

    return header


def get_data(line, name_subs, original_header, new_header):
    line = line.replace(',', ';')
    line_parts = line.split('\t')

    line_parts[0] = name_subs[line_parts[0]]

    original_data = dict(zip(original_header, line_parts))
    new_data = {h: '' for h in new_header}

    for label, value in original_data.items():
        new_data[label] = value

    vir_score = int(original_data['virulence_score'])
    res_score = int(original_data['resistance_score'])
    res_classes = int(original_data['num_resistance_classes'])
    res_genes = int(original_data['num_resistance_genes'])

    new_data['virulence_score__colour'] = get_vir_score_colour(vir_score)
    new_data['resistance_score__colour'] = get_res_score_colour(res_score)
    new_data['num_resistance_classes__colour'] = get_res_classes_colour(res_classes)
    new_data['num_resistance_genes__colour'] = get_res_genes_colour(res_genes)
#    new_data['Yersiniabactin__colour'] = get_vir_lineage_colour(original_data['Yersiniabactin'])
    new_data['Colibactin__colour'] = get_vir_lineage_colour(original_data['Colibactin'])
    new_data['Aerobactin__colour'] = get_vir_lineage_colour(original_data['Aerobactin'])
    new_data['Salmochelin__colour'] = get_vir_lineage_colour(original_data['Salmochelin'])
    new_data['rmpA__colour'] = get_rmpA_colour(original_data['rmpA'])
    new_data['rmpA2__colour'] = get_rmpA2_colour(original_data['rmpA2'])
    for res_class in ['AGly', 'Col', 'Fcyn', 'Flq', 'Gly', 'MLS', 'Ntmdz', 'Phe', 'Rif', 'Sul',
                      'Tet', 'Tgc', 'Tmt', 'Bla', 'Bla_Carb', 'Bla_ESBL', 'Bla_ESBL_inhR', 'Bla_broad',
                      'Bla_broad_inhR' ,'Omp']:
        new_data[res_class + '__colour'] = get_res_class_colour(original_data[res_class])

    return [new_data[h] for h in new_header]


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
    print()
    print('Writing Microreact tree to: {}'.format(tree_out))
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


def get_vir_score_colour(vir_score):
    try:
        return ['#DEEBF7', '#9ECAE1', '#6BAED6', '#4292C6', '#2171B5', '#08306B'][vir_score]
    except IndexError:
        return '#BFBFBF'


def get_res_score_colour(res_score):
    try:
        return ['#FCBBA1', '#FC9272', '#FB6A4A', '#BE413D'][res_score]
    except IndexError:
        return '#BFBFBF'


def get_res_classes_colour(res_classes):
    try:
        return colour_range('#FCBBA1', '#BE413D', 11)[res_classes]
    except IndexError:
        return '#BE413D'


def get_res_genes_colour(res_genes):
    try:
        return colour_range('#FCBBA1', '#BE413D', 21)[res_genes]
    except IndexError:
        return '#BE413D'


def get_species_colour(species):
    try:
        return {'Klebsiella pneumoniae': '#875F9A',
                'Klebsiella variicola subsp. variicola': '#8CBDB2',
                'Klebsiella quasivariicola': '#F0B663',
                'Klebsiella quasipneumoniae subsp. quasipneumoniae': '#ED6060',
                'Klebsiella quasipneumoniae subsp. similipneumoniae': '#EDA483'}[species]
    except IndexError:
        return '#BFBFBF'


def get_vir_lineage_colour(vir_lineage):
    vir_lineage_colours = {'ybt 1': '#B27F91', 'ybt 2': '#CDA12C', 'ybt 3': '#56A354',
                           'ybt 4': '#F28FA2', 'ybt 5': '#DB7723', 'ybt 6': '#93539D',
                           'ybt 7': '#3A85A8', 'ybt 8': '#7B75CC', 'ybt 9': '#D9C5EF',
                           'ybt 10': '#449D72', 'ybt 11': '#EBD930', 'ybt 12': '#6AA3C6',
                           'ybt 13': '#A39F93', 'ybt 14': '#93539D', 'ybt 15': '#EDC59A',
                           'ybt 16': '#840639', 'ybt 17': '#E25065', 'clb 1': '#99BBE0',
                           'clb 2A': '#5972AF', 'clb 2B': '#242F69', 'clb 3': '#242F69',
                           'iro 1': '#B6D5EF', 'iro 2': '#DEC4E8', 'iro 3': '#E29771',
                           'iro 4': '#A4A4EA', 'iro 5': '#E0AAAA', 'iuc 1': '#B6D5EF',
                           'iuc 2': '#DEC4E8', 'iuc 2A': '#D8ABDD', 'iuc 3': '#C3EADB',
                           'iuc 4': '#9ACCBC', 'iuc 5': '#E0AAAA'}
    vir_lineage = vir_lineage.split(';')[0]
    if vir_lineage in vir_lineage_colours:
        return vir_lineage_colours[vir_lineage]
    elif vir_lineage == '-':
        return '#FFFFFF'
    else:
        return '#BFBFBF'


def get_rmpA_colour(rmpA):
    return '#FFFFFF' if rmpA == '-' else '#08306B'
    
def get_rmpA2_colour(rmpA2):
    return '#FFFFFF' if rmpA2 == '-' else '#08306B'


def get_res_class_colour(res_class):
    return '#FFFFFF' if res_class == '-' else '#BE413D'


if __name__ == '__main__':
    main()
