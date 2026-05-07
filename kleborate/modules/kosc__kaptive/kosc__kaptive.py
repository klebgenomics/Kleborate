"""
This module contains classes for interacting with bacterial genome assemblies and contigs and a pipeline
to type them.

Copyright 2026 Mary Maranga, Tom Stanton
https://github.com/klebgenomics/Kleborate/
https://github.com/klebgenomics/Kaptive

This file is part of Kaptive. Kaptive is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Kaptive is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Kaptive.
If not, see <https://www.gnu.org/licenses/>.
"""

import os
import sys
import shutil
import pathlib
import warnings
from pathlib import Path
from Bio import SeqIO, BiopythonParserWarning
warnings.filterwarnings("ignore", category=BiopythonParserWarning)

from kaptive.database import load_database
from kaptive.utils import check_cpus
from kaptive.assembly import typing_pipeline


def description():
    return 'In silico serotyping of K and L locus for the Klebsiella oxytoca species complex'


def prerequisite_modules():
    return []


def get_headers():
    full_headers = [
        'K_locus', 'K_type', 'K_locus_confidence', 'K_locus_problems', 'K_locus_identity',
        'K_Missing_expected_genes',
        'O_locus', 'O_type', 'O_locus_confidence', 'O_locus_problems', 'O_locus_identity',
        'O_Missing_expected_genes'
    ]
    stdout_headers = []
    return full_headers, stdout_headers


def add_cli_options(parser):
    module_name = os.path.basename(__file__)[:-3]
    group = parser.add_argument_group(f'{module_name} module')
    return group


def check_cli_options(args):
    if hasattr(args, 'threads') and args.threads < 1:
        raise ValueError("The number of threads must be at least 1.")


def check_external_programs():
    if not shutil.which('Kaptive'):
        sys.exit('Error: could not find kaptive')
    return ['kaptive']


def data_dir():
    return pathlib.Path(__file__).parents[0] / 'data'

# define all headers

all_headers = [
    'Assembly', 'locus', 'type', 'locus confidence',
    'locus problems', 'locus identity', 'Coverage', 'Length discrepancy',
    'Expected genes in locus', 'Expected genes in locus, details',
    'Missing expected genes', 'Other genes in locus',
    'Other genes in locus, details', 'Expected genes outside locus',
    'Expected genes outside locus, details', 'Other genes outside locus',
    'Other genes outside locus, details', 'Truncated genes, details'
]

def get_results(assembly, minimap2_index, args, previous_results):
    full_headers, _ = get_headers()
    assembly_path = Path(assembly)

    results_dict = {}
    threads = getattr(args, 'threads', 1)

    k_db = load_database(data_dir() / 'KoSC_K_locus_database.gbk')
    o_db = load_database(data_dir() / 'KoSC_O_locus_database.gbk')

    k_results = typing_pipeline(assembly_path, k_db, threads=threads)
    if k_results is not None:
        k_result_table = k_results.format('tsv')
        for line in k_result_table.split('\n'):
            if line:
                parts = line.split('\t')
                for key, value in zip(all_headers, parts):
                    header = 'K_' + key.replace(' ', '_')
                    if header in full_headers:
                        results_dict[header] = value
    else:
        print("Warning: No gene alignments sufficient for typing. Skipping k_results processing.")

    o_results = typing_pipeline(assembly_path, o_db, threads=threads)
    if o_results is not None:
        o_result_table = o_results.format('tsv')
        for line in o_result_table.split('\n'):
            if line:
                parts = line.split('\t')
                for key, value in zip(all_headers, parts):
                    header = 'O_' + key.replace(' ', '_')
                    if header in full_headers:
                        results_dict[header] = value
    else:
        print("Warning: No gene alignments sufficient for typing. Skipping o_results processing.")

    for h in results_dict.keys():
        if h not in full_headers:
            sys.exit(f'Error: results contained a value ({h}) that is not covered by the full headers')

    results_dict = {k: (v if v else '-') for k, v in results_dict.items()}

    return results_dict
