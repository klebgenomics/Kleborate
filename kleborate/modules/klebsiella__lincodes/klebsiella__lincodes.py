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
import os
import re
import tempfile
from pathlib import Path

def description():
    return 'Klebsiella cgMLST and LINcodes with MIST'


def prerequisite_modules():
    return []


def get_headers():
    full_headers = ['cgST', 'LINcodes']
    stdout_headers = []
    return full_headers, stdout_headers


def add_cli_options(parser):
    pass


def check_cli_options(args):
    pass


def check_external_programs():
    if not shutil.which('mist'):
        sys.exit('Error: could not find mist')
    return ['mist']


def data_dir():
    return pathlib.Path(__file__).parents[0] / 'data'


def extract_lincode_from_stdout(stdout):
    """
    Extracts scgST and full or partial LINcode from mist_to_partial_lincode.py
    """
    results = {}

    st_match = re.search(r'Best matching: (scgST-\d+)', stdout)
    lin_match = re.search(r'LINcode for scgMST-\d+: ([\d_]+|-|n/a)', stdout, re.IGNORECASE)
    partial_lin_match = re.search(
        r'Partial LINcode for input strain: ([\d_*]+)',
        stdout
    )

    # Extract scgST
    if st_match:
        try:
            st_value = st_match.group(1).split('-')[1]
            results['cgST'] = int(st_value)
        except (IndexError, ValueError):
            pass

    # Extract LINcode
    if lin_match:
        lin_value = lin_match.group(1)

        # If full LINcode is missing or invalid, fall back to partial
        if lin_value.lower() in {'n/a', '-'}:
            if partial_lin_match:
                results['LINcodes'] = partial_lin_match.group(1)
        else:
            results['LINcodes'] = lin_value

    else:
        # No full LINcode try partial
        if partial_lin_match:
            results['LINcodes'] = partial_lin_match.group(1)

    return results
    

def run_mist_and_extract_lincode(assembly, db_path, mist_script_path):
    """
    Runs MiST and mist_to_partial_lincode.py to extract LINcodes for an assembly, extracts the cgST and Lincode.

    """
    assembly_id = assembly.stem
    with tempfile.TemporaryDirectory() as tempdir:
        json_path = os.path.join(tempdir, f"{assembly_id}.json")
        mist_cmd = [
            "mist", "call",
            "--fasta", str(assembly),
            "--db", str(db_path),
            "--out-json", json_path
        ]
        mist_run = subprocess.run(
            mist_cmd,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        if mist_run.returncode != 0:
            raise subprocess.CalledProcessError(
                mist_run.returncode, mist_cmd,
                output=mist_run.stdout, stderr=mist_run.stderr
            )
        python_cmd = [
            "python", str(mist_script_path), json_path
        ]
        completed = subprocess.run(
            python_cmd,
            capture_output=True,
            text=True,
            check=False
        )
        if completed.returncode != 0:
            raise subprocess.CalledProcessError(
                completed.returncode, python_cmd,
                output=completed.stdout, stderr=completed.stderr
            )
        results = extract_lincode_from_stdout(completed.stdout)
        return results
        

def get_results(assembly, minimap2_index, args, previous_results):
    """
    Returns:
        - "cgST": best matching scgST
        - "LINcodes": full LINcode string
    """
    if isinstance(assembly, str):
        assembly = pathlib.Path(assembly)

    #  Retrieve the DB path from the CLI arguments
    if not args.kpsc_cgmlst_db:
        print("Error: --kpsc_cgmlst_db path is required for kpsc__cgmlst module", file=sys.stderr)
        return {"cgST": "-", "LINcodes": "-"}

    db_path = pathlib.Path(args.kpsc_cgmlst_db)
    mist_script_path = data_dir() / "mist_to_partial_lincode.py"
    if not db_path.exists():
        print(f"Error: The provided database path does not exist: {db_path}", file=sys.stderr)
        return {"cgST": "-", "LINcodes": "-"}

    if not mist_script_path.exists():
        print(f"Error: Helper script not found at {mist_script_path}", file=sys.stderr)
        return {"cgST": "-", "LINcodes": "-"}

    try:
        extracted_data = run_mist_and_extract_lincode(assembly, db_path, mist_script_path)
        
        return {
            "cgST": extracted_data.get('cgST', '-'),
            "LINcodes": extracted_data.get('LINcodes', '-')
        }
    except Exception:
        return {
            "cgST": "-",
            "LINcodes": "-"
        }

# def get_results(assembly, minimap2_index, args, previous_results):
#     """
#     Returns:
#         - "cgST": best matching scgST
#         - "LINcodes": full LINcode string
#     """
#     if isinstance(assembly, str):
#         assembly = pathlib.Path(assembly)
#     db_path = data_dir() / "kleb_scgmlst_s-index"
#     mist_script_path = data_dir() / "mist_to_partial_lincode.py"
#     if not db_path.exists():
#         pass
#     if not mist_script_path.exists():
#         pass
#     try:
#         extracted_data = run_mist_and_extract_lincode(assembly, db_path, mist_script_path)
#         return {
#             "cgST": extracted_data['cgST'],
#             "LINcodes": extracted_data['LINcodes']
#         }
#     except Exception:
#         return {
#             "cgST": "-",
#             "LINcodes": "-"
#         }

