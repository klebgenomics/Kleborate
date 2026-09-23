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
import pandas as pd
import ast

from .shigella_classify import classify_shigella


def description():
    return 'classification of Shigella based on serotype markers'


def prerequisite_modules():
    return []


def get_headers():
    full_headers = ['Serotype']
    stdout_headers = []
    return full_headers, stdout_headers


def add_cli_options(parser):
    module_name = os.path.basename(__file__)[:-3]
    group = parser.add_argument_group(f'{module_name} module')
    group.add_argument('--shigellatyping__min_identity', type=float, default=90.0,
                       help='Minimum alignment percent identity for detecting shigella')
    group.add_argument('--shigellatyping__min_coverage', type=float, default=80.0,
                       help='Minimum alignment percent coverage for detecting shigella')
    return 



def check_cli_options(args):
    if args.shigellatyping__min_identity <= 50.0 or args.shigellatyping__min_identity >= 100.0:
        sys.exit('Error: --shigellatyping__min_identity must be between 50.0 and 100.0')
    if args.shigellatyping__min_coverage <= 50.0 or args.shigellatyping__min_coverage >= 100.0:
        sys.exit('Error: --shigellatyping__min_coverage must be between 50.0 and 100.0')
    

def check_external_programs():
    if not shutil.which('minimap2'):
        sys.exit('Error: could not find minimap2')
    return ['minimap2']


def data_dir():
    return pathlib.Path(__file__).parents[0] / 'data'



def get_results(assembly, minimap2_index, args, previous_results):

    ShigellaRef = data_dir() / 'ShigellaRef5.fasta'

    # Load the shigella_serotype_markers dictionary from the file
    with open(data_dir() / 'shigella_serotype_markers.txt', 'r') as file:
        file_content = file.read()
        shigella_serotype_markers = ast.literal_eval(file_content)


    Serotype = classify_shigella(
        assembly,
        minimap2_index,
        ShigellaRef,
        args.shigellatyping__min_identity,
        args.shigellatyping__min_coverage,
        shigella_serotype_markers
    )

    results = {
        'Serotype': Serotype
    }


    return results 


