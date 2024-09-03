"""
Copyright 2024 Kat Holt
Copyright 2024 Ryan Wick (rrwick@gmail.com)
Copyright 2024 (gathonimaranga@gmail.com)
https://github.com/katholt/Kleborate/

This file is part of Kleborate. Kleborate is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Kleborate is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Kleborate. If
not, see <https://www.gnu.org/licenses/>.
"""
import os
import pathlib
import collections
from collections import defaultdict

from ...shared.resMinimap import read_class_file, get_res_headers

def description():
    return 'resistance score (0-3) for the Klebsiella pneumoniae species complex, based on the ' \
           'results of the klebsiella_pneumo_complex__amr module'


def prerequisite_modules():
    return ['klebsiella_pneumo_complex__amr']


def get_headers():
    full_headers = ['resistance_score']
    stdout_headers = ['resistance_score']
    return full_headers, stdout_headers


def add_cli_options(parser):
    pass


def check_cli_options(args):
    pass


def check_external_programs():
    return []

def data_dir():
    return pathlib.Path(__file__).parents[1] / 'klebsiella_pneumo_complex__amr' / 'data'


def get_results(assembly, minimap2_index, args, previous_results):
    _, res_classes, bla_classes = read_class_file(data_dir() / 'CARD_AMR_clustered.csv')
    res_headers = get_res_headers(res_classes, bla_classes)

    """
    Four possible resistance scores:
      * 0 = no ESBL, no carbapenemase (regardless of colistin resistance)
      * 1 = ESBL, no carbapenemase (regardless of colistin resistance)
      * 2 = Carbapenemase without colistin resistance
      * 3 = Carbapenemase and colistin resistance
    """

    #res_hits = {key.split('__')[1]: value for key, value in previous_results.items() if key.startswith('klebsiella_pneumo_complex__amr__')}
    res_hits = {key.split('__')[2]: value for key, value in previous_results.items() if key.startswith('klebsiella_pneumo_complex__amr__')}

    if not res_headers:
        return '-'
    
    # Look for a hit in any 'ESBL' column (e.g. 'Bla_ESBL' or 'Bla_ESBL_inhR').
    esbl_header_indices = [i for i, h in enumerate(res_headers) if 'esbl' in h.lower()]
    has_esbl = any(res_hits[res_headers[i]] != '-' for i in esbl_header_indices)

    # Look for a hit in any 'Carb' column (e.g. 'Bla_Carb').
    carb_header_indices = [i for i, h in enumerate(res_headers) if 'carb' in h.lower()]
    has_carb = any(res_hits[res_headers[i]] != '-' for i in carb_header_indices)

    # Look for a hit in the 'Col' columns.
    col_header_indices = [i for i, h in enumerate(res_headers)
                      if h.lower() == 'col_acquired' or h.lower() == 'col_mutations']
    has_col = any(res_hits[res_headers[i]] != '-' for i in col_header_indices)
    
    if has_carb and has_col:
        return {'resistance_score': '3'}
    elif has_carb:
        return {'resistance_score': '2'}
    elif has_esbl:
        return {'resistance_score': '1'}
    else:
        return {'resistance_score': '0'}
    







