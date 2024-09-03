"""
Copyright 2023 Kat Holt
Copyright 2023 Ryan Wick (rrwick@gmail.com)
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
import shutil
import sys

from ...shared.multi_mlst import multi_mlst
from ...shared.alignment import truncation_check


def description():
    return 'MLST on the KpSC yersiniabactin locus (ybt and irp genes)'


def prerequisite_modules():
    return []


def get_headers():
    full_headers = ['YbST', 'Yersiniabactin', 'ybtS', 'ybtX', 'ybtQ', 'ybtP', 'ybtA', 'irp2', 'irp1',
                    'ybtU', 'ybtT', 'ybtE', 'fyuA', 'spurious_ybt_hits']
    stdout_headers = []
    return full_headers, stdout_headers


def add_cli_options(parser):
    module_name = os.path.basename(__file__)[:-3]
    group = parser.add_argument_group(f'{module_name} module')
    group.add_argument('--klebsiella__ybst_min_identity', type=float, default=90.0,
                       help='Minimum alignment percent identity for yersiniabactin MLST')
    group.add_argument('--klebsiella__ybst_min_coverage', type=float, default=80.0,
                       help='Minimum alignment percent coverage for yersiniabactin MLST')
    group.add_argument('--klebsiella__ybst_min_spurious_identity', type=float, default=80.0,
                       help='Minimum alignment percent identity for klebsiella__ybst spurious results')
    group.add_argument('--klebsiella__ybst_min_spurious_coverage', type=float, default=40.0,
                       help='Minimum alignment percent coverage for klebsiella__ybst spurious results')
    group.add_argument('--klebsiella__ybst_required_exact_matches', type=int, default=6,
                       help='At least this many exact matches are required to call an ST')
    return group


def check_cli_options(args):
    if args.klebsiella__ybst_min_identity <= 50.0 or args.klebsiella__ybst_min_identity >= 100.0:
        sys.exit('Error: --klebsiella__ybst_min_identity must be between 50.0 and 100.0')
    if args.klebsiella__ybst_min_coverage <= 50.0 or args.klebsiella__ybst_min_coverage >= 100.0:
        sys.exit('Error: --klebsiella__ybst_min_coverage must be between 50.0 and 100.0')
    if args.klebsiella__ybst_min_spurious_identity <= 50.0 or args.klebsiella__ybst_min_spurious_identity >= 100.0:
        sys.exit('Error: --klebsiella__ybst_min_spurious_identity must be between 50.0 and 100.0')
    if args.klebsiella__ybst_min_spurious_coverage <= 30.0 or args.klebsiella__ybst_min_spurious_coverage >= 100.0:
        sys.exit('Error: --klebsiella__ybst_min_spurious_coverage must be between 30.0 and 100.0')
    if args.klebsiella__ybst_required_exact_matches < 0:
        sys.exit('Error: --klebsiella__ybst_required_exact_matches must be a positive integer')


def check_external_programs():
    if not shutil.which('minimap2'):
        sys.exit('Error: could not find minimap2')
    return ['minimap2']


def data_dir():
    return pathlib.Path(__file__).parents[0] / 'data'


def get_results(assembly, minimap2_index, args, previous_results):
    genes = ['ybtS', 'ybtX', 'ybtQ', 'ybtP', 'ybtA', 'irp2', 'irp1', 'ybtU', 'ybtT', 'ybtE', 'fyuA']
    profiles = data_dir() / 'profiles.tsv'
    alleles = {gene: data_dir() / f'{gene}.fasta' for gene in genes}

    results, spurious_hits = multi_mlst(assembly, minimap2_index, profiles, alleles, genes,
                                      'lineage_ICE', args.klebsiella__ybst_min_identity,
                                      args.klebsiella__ybst_min_coverage, args.klebsiella__ybst_required_exact_matches,
                                      check_for_truncation=True, report_incomplete=True, 
                                      min_spurious_identity=args.klebsiella__ybst_min_spurious_identity,
                                      min_spurious_coverage=args.klebsiella__ybst_min_spurious_coverage)
    st, lineage, alleles = results

    if st == 'NA':
        st = 0
    else:
        st = st[2:]

    # spurious hits
    spurious_hits = [item for h in spurious_hits.values() for item in h]

    spurious_virulence_hits = ';'.join(spurious_hits )if spurious_hits else '-'

    return {'YbST': st, 'Yersiniabactin': lineage,
            'ybtS': alleles['ybtS'], 'ybtX': alleles['ybtX'], 'ybtQ': alleles['ybtQ'],
            'ybtP': alleles['ybtP'], 'ybtA': alleles['ybtA'], 'irp2': alleles['irp2'],
            'irp1': alleles['irp1'], 'ybtU': alleles['ybtU'], 'ybtT': alleles['ybtT'],
            'ybtE': alleles['ybtE'], 'fyuA': alleles['fyuA'],
            'spurious_ybt_hits':spurious_virulence_hits}
