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
import pathlib
import shutil
import sys
import subprocess
import os
import re
import tempfile
from pathlib import Path

def description():
    return 'Klebsiella cgMLST and LINcodes typing with MIST'


def prerequisite_modules():
    return []


def get_headers():
    full_headers = ['cgST', 'LIN code','Sublineage', 'Clonal group']
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

    results = {
        'cgST': '-',
        'LINcodes': '-',
        'Clonal group': '-',
        'Sublineage': '-'
    }

    st_matches = re.findall(r'(?:Best matching: scgST-|Profile: scgST-)(\d+)', stdout)
    lin_match = re.search(r'LINcode for scgMST-\d+: ([\d_]+|-|n/a)', stdout, re.IGNORECASE)
    partial_lin_match = re.search(r'Partial LINcode for input strain: ([\d_*]+)', stdout)
    
    cg_matches = re.findall(r'Clonal group:\s*(.+)', stdout)
    sl_matches = re.findall(r'Sublineage:\s*(.+)', stdout)

    if st_matches:
        results['cgST'] = "; ".join(dict.fromkeys(st_matches))

    lin_value = lin_match.group(1) if lin_match else 'n/a'
    
    if lin_value.lower() in {'n/a', '-'} or not lin_match:
        if partial_lin_match:
            results['LINcodes'] = partial_lin_match.group(1)
    else:
        results['LINcodes'] = lin_value

    if cg_matches:
        vals = [v.strip() for v in cg_matches if v.strip().lower() != 'n/a']
        results['Clonal group'] = "; ".join(dict.fromkeys(vals)) if vals else '-'
        
    if sl_matches:
        vals = [v.strip() for v in sl_matches if v.strip().lower() != 'n/a']
        results['Sublineage'] = "; ".join(dict.fromkeys(vals)) if vals else '-'

    return results


def run_mist_and_extract_lincode(assembly, db_path, mist_script_path):
    """
    Runs MiST and mist_to_partial_lincode.py to extract LINcodes for an assembly, extracts the cgST and Lincode

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

    if isinstance(assembly, str):
        assembly = pathlib.Path(assembly)
        
    db_path = data_dir() / "kleb_scgmlst_s-index"
    mist_script_path = data_dir() / "mist_to_partial_lincode.py"
    
    try:
        extracted_data = run_mist_and_extract_lincode(assembly, db_path, mist_script_path)
        
        raw_st = extracted_data.get('cgST', '-')
        
        if raw_st != '-':
            formatted_cgst = "; ".join([f"cgST{st.strip()}" for st in raw_st.split(";")])
        else:
            formatted_cgst = "-"
        
        return {
            "cgST": formatted_cgst,
            "LIN code": extracted_data.get('LINcodes', '-'),
            "Sublineage": extracted_data.get('Sublineage', '-'),
            "Clonal group": extracted_data.get('Clonal group', '-')
        }
    except Exception as e:
        return {
            "cgST": "-",
            "LIN code": "-",
            "Sublineage": "-",
            "Clonal group": "-"
        }

