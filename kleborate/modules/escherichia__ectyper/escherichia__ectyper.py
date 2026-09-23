"""
Copyright 2026 Mary Maranga (gathonimaranga@gmail.com)
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
import sys
import shutil
from pathlib import Path
import pandas as pd
import tempfile
from collections import defaultdict
import subprocess
import csv
import re


def description():
    return 'E. coli serotyping of O and H antigens'


def prerequisite_modules():
    return []


def get_headers():

    full_headers = [
        'O-type','H-type','Serotype','QC','Evidence','GeneScores',
        'AlleleKeys','GeneIdentities(%)','GeneCoverages(%)','GeneLengths','Warnings'
    ]
    stdout_headers = []
    
    return full_headers, stdout_headers


def add_cli_options(parser):
    module_name = os.path.basename(__file__)[:-3]
    group = parser.add_argument_group(f'{module_name} module')
    group.add_argument('--escherichia__serotyping_cores', type=int, default=4,
                       help='Number of CPU cores to use for ectyper (default: 4)')
    return group


    
def check_cli_options(args):
    """
    Validate the command-line arguments.
    """

    if args.escherichia__serotyping_cores < 4:
        raise ValueError("The number of threads must be at least 4")

    if not shutil.which('ectyper'):
        sys.exit('Error: ectyper is not installed or not in PATH.')


def check_external_programs():
    """
    Ensure the required external programs are available.
    """
    if not shutil.which('ectyper'):
        sys.exit('Error: could not find ectyper executable.')
    return ['ectyper']


def data_dir():
    return pathlib.Path(__file__).parents[0] / 'data'


def run_ectyper(input_fasta, output_dir, quiet,cores) :
    """
    Run ectyper to serotype E. coli from a FASTA file.

    Parameters:
        input_fasta (str): Path to the input FASTA file.
        output_dir (str): Directory where ectyper will write its output.
    """


    os.makedirs(output_dir, exist_ok=True)
    command = [
        "ectyper",
        "-i", input_fasta,
        "-o", output_dir,
        "-c", str(cores)
    ]
    if quiet:
        stdout_dest = subprocess.DEVNULL
        stderr_dest = subprocess.DEVNULL
    else:
        stdout_dest = subprocess.PIPE
        stderr_dest = subprocess.PIPE

    subprocess.run(
        command,
        stdout=stdout_dest,
        stderr=stderr_dest,
        check=True,
        text=not quiet
    )

def get_results(assembly, ref_index, args, previous_results):
    results = {}
    quiet = getattr(args, 'quiet', True)

    with tempfile.TemporaryDirectory() as tmpdir:
        run_ectyper(assembly, tmpdir, quiet=quiet, cores=args.escherichia__serotyping_cores)

        tsv_path = os.path.join(tmpdir, "output.tsv")
        full_headers, _ = get_headers()

        if not os.path.exists(tsv_path):
            if not quiet:
                print(f"[get_results] output.tsv not found in temporary directory")
            return {}

        with open(tsv_path, newline='') as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            try:
                row = next(reader)
            except StopIteration:
                return {}

            for col in full_headers:
                raw = row.get(col, "") or ""

                raw = re.sub(r"^-:", "", raw)
                clean = raw.rstrip(";").strip()

                results[col] = clean

    return results


