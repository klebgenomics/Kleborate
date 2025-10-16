
"""
Copyright 2025 Mary Maranga
https://github.com/klebgenomics/Kleborate/

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

from ...shared.alignment import align_query_to_ref, truncation_check
from ...shared.misc import load_fasta, reverse_complement



def description():
    return 'typing for the peg-344 gene'


def prerequisite_modules():
    return []


def get_headers():
    full_headers = ['peg-344']
    stdout_headers = []
    return full_headers, stdout_headers



def add_cli_options(parser):
    module_name = os.path.basename(__file__)[:-3]
    group = parser.add_argument_group(f'{module_name} module')

    group.add_argument('--klebsiella__peg_min_identity', type=float, default=90.0,
                       help='Minimum alignment percent identity for klebsiella__peg-344 results')
    group.add_argument('--klebsiella__peg_min_coverage', type=float, default=80.0,
                       help='Minimum alignment percent coverage for klebsiella__peg-344 results')

    return group


def check_cli_options(args):
    if args.klebsiella__peg_min_identity <= 50.0 or args.klebsiella__peg_min_identity >= 100.0:
        sys.exit('Error: --klebsiella__peg-344_min_identity must be between 50.0 and 100.0')
    if args.klebsiella__peg_min_coverage <= 50.0 or args.klebsiella__peg_min_coverage >= 100.0:
        sys.exit('Error: --klebsiella__peg-344_min_coverage must be between 50.0 and 100.0')


def check_external_programs():
    if not shutil.which('minimap2'):
        sys.exit('Error: could not find minimap2')
    return ['minimap2']


def data_dir():
    return pathlib.Path(__file__).parents[0] / 'data'


def check_peg(ref_file, assembly, minimap2_index, min_coverage, min_identity):
    """
    Checks presence of peg-344 gene.
    """
    alignment_hits = align_query_to_ref(
        ref_file,  
        assembly, 
        ref_index=minimap2_index,
        min_identity=min_identity,
        min_query_coverage=min_coverage
    )

    if not alignment_hits:
        return '-'
    
    for hit in alignment_hits:
        _, coverage, _ = truncation_check(hit)
        if coverage >= 100.0:
            return 'present'
        else:
            return f'truncated ({coverage}%)'


def get_results(assembly, minimap2_index, args, previous_results):
    ref_file = data_dir() / 'peg-344.fasta'
    
    result = check_peg(
        ref_file, 
        assembly,
        minimap2_index, 
        args.klebsiella__peg_min_coverage,
        args.klebsiella__peg_min_identity)

    if not result:
        return {'peg-344': '-'}

    else:
        return {'peg-344': result}
   
