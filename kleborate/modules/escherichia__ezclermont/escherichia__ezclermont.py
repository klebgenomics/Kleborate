
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
import io
from contextlib import redirect_stdout, redirect_stderr
import itertools
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from . import run


def description():
    return ' Clermont PCR typing method for in silico analysis of E. coli whole genomes or assembled contigs'


def prerequisite_modules():
    return []


def get_headers():
    full_headers = ['Clermont_type', 'Clermont_profile']
    stdout_headers = []
    return full_headers, stdout_headers



def add_cli_options(parser):
    module_name = os.path.basename(__file__)[:-3]
    group = parser.add_argument_group(f'{module_name} module')
    group.add_argument('--escherichia__ezclermont_min_length', type=float, default=500.0,
                       help='Minimum contig length to consider in the analysis')

    return


def check_cli_options(args):
    if args.escherichia__ezclermont_min_length <= 400.0:
        sys.exit('Error: --escherichia__ezclermont_min_length must be between 400.0 and 500.0')
   

def check_external_programs():
    if not shutil.which('minimap2'):
        sys.exit('Error: could not find minimap2')
    return ['minimap2']


def data_dir():
    return pathlib.Path(__file__).parents[0] / 'data'


def get_results(assembly, minimap2_index, args, previous_results):

    args.contigs = assembly
    args.min_length = args.escherichia__ezclermont_min_length if hasattr(args, 'min_length') else 500

    for attr, default in [('logfile', None), ('no_partial', False), ('experiment_name', None)]:
        setattr(args, attr, getattr(args, attr, default))


    # Run the PCR analysis from the imported `run` module
    f_stdout = io.StringIO()
    f_stderr = io.StringIO()
    with redirect_stdout(f_stdout), redirect_stderr(f_stderr):
        clermont_type, profile = run.main(args)

    profile = profile.replace('\n', '; ').strip('; ')

    # Create a results dictionary
    results = {
        "Clermont_type": clermont_type,
        "Clermont_profile": profile
    }

    return results
