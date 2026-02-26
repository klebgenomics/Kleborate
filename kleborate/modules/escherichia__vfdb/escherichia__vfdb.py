"""
Copyright 2025 Mary Maranga (gathonimaranga@gmail.com)
https://github.com/klebgenomics/Kleborate

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
import shutil
import sys
import subprocess
from pathlib import Path
from typing import Dict, Tuple
from typing import Dict

from .vfdb import map_virulence_factors, extract_fasta_headers



def description():
    return 'Typing of virulence factors in E. coli genomes'


def prerequisite_modules():
    return []


def data_dir():

    return pathlib.Path(__file__).parents[0] / 'data'


_NOTES_CACHE = None

def load_notes_mapping():
    global _NOTES_CACHE
    if _NOTES_CACHE is not None:
        return _NOTES_CACHE

    mapping = {}
    notes_path = data_dir() / 'notes.txt'
    if notes_path.exists():
        with notes_path.open('r', encoding='utf-8') as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith('#') or ':' not in line:
                    continue
                key, value = line.split(':', 1)
                key = key.strip()
                value = value.strip().rstrip(':').strip()
                if key and value:
                    mapping[key] = value

    _NOTES_CACHE = mapping
    return mapping


def display_name(short_id, notes):
    desc = notes.get(short_id)
    return f'{short_id} ({desc})' if desc else short_id


def get_headers():
    ref_file = data_dir() / 'virulence_ecoli.fsa'
    notes = load_notes_mapping()

    raw_headers = sorted(extract_fasta_headers(ref_file))
    full_headers = [display_name(h, notes) for h in raw_headers]
    stdout_headers = []
    return full_headers, stdout_headers

# def get_headers():
#     full_headers = list(extract_fasta_headers(data_dir() / 'virulence_ecoli.fsa'))
#     stdout_headers = []
#     return full_headers, stdout_headers


def add_cli_options(parser):
    module_name = os.path.basename(__file__)[:-3]
    group = parser.add_argument_group(f'{module_name} module')
    group.add_argument('--escherichia__vfdb_min_identity', type=float, default=90.0,
                       help='Minimum alignment percent identity for detecting virulence factors')
    group.add_argument('--escherichia__vfdb_min_coverage', type=float, default=80.0,
                       help='Minimum alignment percent coverage for detecting virulence factors')
    return group


def check_cli_options(args):
    if args.escherichia__pathovar_min_identity <= 50.0 or args.escherichia__pathovar_min_identity >= 100.0:
        sys.exit('Error: --escherichia__vfdb_min_identity must be between 50.0 and 100.0')
    if args.escherichia__pathovar_min_coverage <= 50.0 or args.escherichia__pathovar_min_coverage >= 100.0:
        sys.exit('Error: --escherichia__vfdb_min_coverage must be between 50.0 and 100.0')

def check_external_programs():
    if not shutil.which('minimap2'):
        sys.exit('Error: could not find minimap2')
    return ['minimap2']


def get_results(assembly, minimap2_index, args, previous_results):
    
    ref_file = data_dir() / 'virulence_ecoli.fsa'
    notes = load_notes_mapping()
    raw_headers = sorted(extract_fasta_headers(ref_file))

    # Run VF mapping
    virulence_markers = map_virulence_factors(
        assembly,
        minimap2_index,
        ref_file,
        args.escherichia__vfdb_min_identity,
        args.escherichia__vfdb_min_coverage
    )

    results = {}
    for short_id in raw_headers:
        key = display_name(short_id, notes)
        results[key] = virulence_markers.get(short_id, '-')

    return results

# def get_results(assembly, minimap2_index, args, previous_results):
#     full_headers, _ = get_headers()
#     ref_file = data_dir() / 'virulence_ecoli.fsa'
    
#     virulence_markers = map_virulence_factors(
#         assembly,
#         minimap2_index,
#         ref_file,
#         args.escherichia__vfdb_min_identity,
#         args.escherichia__vfdb_min_coverage
#     )

    
#     results = {header: virulence_markers.get(header, '-') for header in full_headers}

#     return results



