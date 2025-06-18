
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
import subprocess
import re

# from . import run


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
    group.add_argument('--escherichia__ezclermont_min_length', type=int, default=500,
                       help='Minimum contig length to consider in the analysis (integer, 400–500)')
    return group

def check_cli_options(args):
    if not (400 <= args.escherichia__ezclermont_min_length <= 500):
        sys.exit('Error: --escherichia__ezclermont_min_length must be between 400 and 500')


def check_external_programs():
    if not shutil.which('ezclermont'):
        sys.exit('Error: could not find ezclermont')
    return ['ezclermont']


def data_dir():
    return pathlib.Path(__file__).parents[0] / 'data'


def run_ezclermont(input_fasta, min_length):
    """
    Run ezClermont.

    Parameters:
        input_fasta (str): Path to the input FASTA file.
        min_length (int): Minimum contig length.

    Returns:
        str: Output from ezClermont.
    """
    command = [
        "ezclermont",
        "-m", str(min_length),
        input_fasta
    ]
    
    try:
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True)
        output = result.stdout or result.stderr
        return output
    except subprocess.CalledProcessError as e:
        output = e.stderr
        if "Clermont type:" in output or "Results" in output:
            return output
        print(f"Error occurred: {e}")
        return None



def get_results(assembly, minimap2_index, args, previous_results):
    min_length = args.escherichia__ezclermont_min_length
    output = run_ezclermont(assembly, min_length)
    if not output:
        return {"Clermont_type": '', "Clermont_profile": ''}

    # Extract Clermont type
    type_match = re.search(r'Clermont type:\s*([A-Za-z0-9]+)', output)
    clermont_type = type_match.group(1) if type_match else ''

    # Extract Clermont profile
    markers = ['TspE4', 'arpA', 'chu', 'yjaA']
    profile_lines = []
    for marker in markers:
        m = re.search(rf'{marker}:\s*[+-]', output)
        if m:
            profile_lines.append(m.group(0))
    clermont_profile = '; '.join(profile_lines)

    results = {
        "Clermont_type": clermont_type,
        "Clermont_profile": clermont_profile
    }
    return results