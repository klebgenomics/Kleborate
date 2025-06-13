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
import sys
import shutil
from pathlib import Path
import subprocess


def description():
    return 'Determining stx type from E.coli genomes using STXTyper'


def prerequisite_modules():
    return []


def get_headers():
    """
    Define the headers for STXTyper results.
    """
    full_headers = [
        'stx_type', 'operon', 'identity', 'target_start', 'target_stop', 
        'A_reference', 'A_identity', 'A_reference_subtype', 'A_coverage',
        'B_reference', 'B_reference_subtype', 'B_identity', 'B_coverage'
    ]
    stdout_headers = []
    return full_headers, stdout_headers


def add_cli_options(parser):
    """
    Add command-line options for this module.
    """
    module_name = os.path.basename(__file__)[:-3]
    group = parser.add_argument_group(f'{module_name} module')
    
    # group.add_argument(
    #     '-n', '--nucleotide', required=True,
    #     help="Assembled nucleotide sequence to search in FASTA format."

    # )

    return group


def check_cli_options(args):
    """
    Validate the command-line arguments.
    """
    # Check if stxtyper is in PATH
    if not shutil.which('stxtyper'):
        sys.exit('Error: STXTyper is not installed or not in PATH.')

def check_external_programs():
    """
    Ensure the required external programs are available.
    """
    if not shutil.which('stxtyper'):
        sys.exit('Error: could not find STXTyper executable.')
    return ['stxtyper']



def run_stxtyper(nucleotide_file):
    """
    Run STXTyper.

    Parameters:
        nucleotide_file (str): Path to the nucleotide FASTA file.

    Returns:
        str: Output from STXTyper.
    """
    command = [
        "stxtyper",
        "-n", nucleotide_file
    ]
    
    try:
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e}")
        return None


# all headers for stxtyper
all_headers = [
        '#target_contig', 'stx_type', 'operon', 'identity', 
        'target_start', 'target_stop', 'target_strand', 
        'A_reference',  'A_reference_subtype', 'A_identity','A_coverage',
        'B_reference', 'B_reference_subtype', 'B_identity', 'B_coverage'
    ]


def parse_stxtyper_output(output, headers, all_headers):
    """
    Parse STXTyper output and concatenate values for each selected column.
    Returns a dict mapping each header to a semicolon-separated string.
    """
    lines = [line.strip() for line in output.strip().split('\n') if line.strip() and not line.startswith('#')]
    if not lines:
        return {header: '' for header in headers}

    # Split each line into columns
    rows = [line.split('\t') for line in lines]

    # Ensure that all rows have the correct number of columns
    expected_col_count = len(all_headers)
    filtered_rows = [row for row in rows if len(row) >= expected_col_count]
    if not filtered_rows:
        return {header: '' for header in headers}

    # Create a mapping of header to column index
    header_to_index = {header: idx for idx, header in enumerate(all_headers)}

    result_dict = {}
    for header in headers:
        idx = header_to_index.get(header)
        if idx is not None:
            # Collect the relevant column from each row, join with ;
            values = [row[idx] for row in filtered_rows]
            result_dict[header] = ';'.join(values)
        else:
            result_dict[header] = ''
    return result_dict



def get_results(assembly, index, previous_results, args):
    """
    Runs STXTyper for a single assembly and parses the output.
    Returns a dict mapping each header to a semicolon-separated string.
    """
    full_headers, _ = get_headers()

    raw_output = run_stxtyper(
        nucleotide_file=assembly
    )

    # print(raw_output)

    if raw_output:
        results = parse_stxtyper_output(raw_output, full_headers, all_headers)
        return results

    return {header: '' for header in full_headers}


