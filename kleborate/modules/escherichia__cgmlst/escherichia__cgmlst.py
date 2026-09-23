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

import json
import os
import pathlib
import shutil
import sys
import subprocess
import tempfile
from pathlib import Path


def description():
    return 'E.coli cgMLST and LINcodes typing with MIST (LIN codes retrieved from EnteroBase)'


def prerequisite_modules():
    return []


def get_headers():

    full_headers = ['cgST', 'LIN code']
    stdout_headers = []
    return full_headers, stdout_headers


ENTEROBASE_TOKEN_ENV_VAR = 'KLEBORATE_ENTEROBASE_TOKEN'
DEFAULT_TOKEN_PATH = Path.home() / '.enterobase_token'


def add_cli_options(parser):
    group = parser.add_argument_group('E. coli cgMLST/LIN code (EnteroBase)')
    group.add_argument(
        '--ecoli_entero_token', type=str, default=None,
        help='Path to an EnteroBase API token file, required to retrieve LIN codes for '
             'E. coli cgSTs '
    )
    group.add_argument(
        '--ecoli_entero_preset', type=str, default='ecoli',
        help="EnteroBase preset passed to 'mist lincode' "
    )


def resolve_entero_token_path(cli_value):
    if cli_value:
        return Path(cli_value).expanduser()
    env_value = os.environ.get(ENTEROBASE_TOKEN_ENV_VAR)
    if env_value:
        return Path(env_value).expanduser()
    return DEFAULT_TOKEN_PATH



def check_cli_options(args):
    token_path = resolve_entero_token_path(args.ecoli_entero_token)

    if not token_path.exists():
        sys.exit(f'Error: EnteroBase API token file not found at {token_path}')

    args.ecoli_entero_token = str(token_path)



def check_external_programs():
    if not shutil.which('mist'):
        sys.exit('Error: could not find mist')

    db_path = data_dir() / "ecoli_cgmlst_v1-index"
    if not db_path.exists() or not any(db_path.iterdir()):
        sys.exit(f'Error: MiST cgMLST database not found at {db_path}')

    return ['mist']



def data_dir():
    return pathlib.Path(__file__).parents[0] / 'data'


def parse_lincode_json(lincode_json_path):
    """
    Parses the JSON file written by `mist lincode` into this module's results dict.
    """
    results = {'cgST': '-', 'LIN code': '-'}

    with open(lincode_json_path) as f:
        data = json.load(f)

    st = data.get('st')
    if st:
        results['cgST'] = f'cgST{st}'

    lincode_full = data.get('lincode_full')
    if lincode_full:
        results['LIN code'] = '-'.join(str(v) for v in lincode_full)

    return results


def run_mist_and_get_lincode(assembly, db_path, entero_token, entero_preset, timeout=300):
    """
    Runs `mist call` for the cgST, then `mist lincode` for the matching LIN code.
    """
    if not db_path.exists() or not any(db_path.iterdir()):
        raise FileNotFoundError(f'MiST database not found at {db_path}')

    assembly_id = assembly.stem
    with tempfile.TemporaryDirectory() as tempdir:
        call_json_path = os.path.join(tempdir, f"{assembly_id}.json")
        mist_call_cmd = [
            "mist", "call",
            "--fasta", str(assembly),
            "--db", str(db_path),
            "--out-json", call_json_path
        ]
        call_run = subprocess.run(
            mist_call_cmd,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=timeout
        )
        if call_run.returncode != 0:
            raise subprocess.CalledProcessError(
                call_run.returncode, mist_call_cmd,
                output=call_run.stdout, stderr=call_run.stderr
            )

        lincode_json_path = os.path.join(tempdir, f"{assembly_id}_lincode.json")
        mist_lincode_cmd = [
            "mist", "lincode",
            call_json_path,
            "--db", str(db_path),
            "--entero-token", str(entero_token),
            "--entero-preset", entero_preset,
            "--output", lincode_json_path
        ]
        lincode_run = subprocess.run(
            mist_lincode_cmd,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=timeout
        )
        if lincode_run.returncode != 0:
            raise subprocess.CalledProcessError(
                lincode_run.returncode, mist_lincode_cmd,
                output=lincode_run.stdout, stderr=lincode_run.stderr
            )

        return parse_lincode_json(lincode_json_path)



def get_results(assembly, ref_index, args, previous_results):

    if isinstance(assembly, str):
        assembly = pathlib.Path(assembly)

    db_path = data_dir() / "ecoli_cgmlst_v1-index"

    try:
        return run_mist_and_get_lincode(
            assembly, db_path, args.ecoli_entero_token, args.ecoli_entero_preset)
    except Exception:
        return {
            "cgST": "-",
            "LIN code": "-",
        }