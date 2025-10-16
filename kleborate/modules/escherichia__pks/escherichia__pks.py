"""
Copyright 2025 Mary Maranga (gathonimaranga@gmailcom)
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

from ...shared.alignment import align_query_to_ref


def description():
    return 'Typing for colibactin polyketide biosynthesis gene cluster'


def prerequisite_modules():
    return []


def get_headers():
    full_headers = ['clbB']
    stdout_headers = []
    return full_headers, stdout_headers


def add_cli_options(parser):
    module_name = os.path.basename(__file__)[:-3]
    group = parser.add_argument_group(f'{module_name} module')
    group.add_argument('--escherichia__pks_min_identity', type=float, default=90.0,
                       help='Minimum alignment percent identity for detecting clbB')
    group.add_argument('--escherichia__pks_min_coverage', type=float, default=80.0,
                       help='Minimum alignment percent coverage for detecting clbB')
    return group


def check_cli_options(args):
    if args.escherichia__pks_min_identity <= 50.0 or args.escherichia__pks_min_identity >= 100.0:
        sys.exit('Error: --escherichia__pks_min_identity must be between 50.0 and 100.0')
    if args.escherichia__pks_min_coverage <= 50.0 or args.escherichia__pks_min_coverage >= 100.0:
        sys.exit('Error: --escherichia__pks_min_coverage must be between 50.0 and 100.0')



def check_external_programs():
    if not shutil.which('minimap2'):
        sys.exit('Error: could not find minimap2')
    return ['minimap2']


def data_dir():
    return pathlib.Path(__file__).parents[0] / 'data'


def pks_minimap(assembly, minimap2_index, ref_file, min_identity, min_coverage):
    alignment_hits = align_query_to_ref(
        ref_file,
        assembly,
        ref_index=minimap2_index,
        min_identity=min_identity,
        min_query_coverage=min_coverage
    )

    if alignment_hits:
        for hit in alignment_hits:
            alignment_length = hit.query_end - hit.query_start
            coverage = (alignment_length / hit.query_length) * 100
            allele = hit.query_name
            if coverage < 100.0 or hit.percent_identity < 100.0:
                return f"{allele}?"
            else:
                return allele
    else:
        return ''

# def pks_minimap(assembly, minimap2_index, ref_file, min_identity, min_coverage):
#     alignment_hits = align_query_to_ref(
#         ref_file,
#         assembly,
#         ref_index=minimap2_index,
#         min_identity=min_identity,
#         min_query_coverage=min_coverage
#     )
#     if alignment_hits:
#         allele = alignment_hits[0].query_name
#         return allele
#     else:
#         return ''



def get_results(assembly, minimap2_index, args, previous_results):
    clbB_ref = data_dir() / 'clbB.fasta'

    result = pks_minimap(
        assembly,
        minimap2_index,
        clbB_ref,
        args.escherichia__pks_min_identity,
        args.escherichia__pks_min_coverage
    )

    if not result:
        return {'clbB': '-'}
    else:
        return {'clbB': result}

