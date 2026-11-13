"""
Copyright 2026 Mary Maranga (gathonimaranga@gmail.com)
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
import re
import sys
import shutil
import tempfile
import subprocess
import warnings
from pathlib import Path
from collections import defaultdict
import pandas as pd
from functools import lru_cache
warnings.filterwarnings('ignore', message=r'Skipping \d+ records with INTERNAL_STOP', category=UserWarning)


def description():
    return 'Antimicrobial resistance (AMR) gene detection in Klebsiella oxytoca using AMRFinderPlus and hAMRonization'


def prerequisite_modules():
    return []


def get_headers():
    """
    Define the headers for AMRFinderPlus results.
    """
    full_headers = [
        'Aminoglycosides', 'Quinolones', 'Fosfomycin', 'Sulfonamides', 'Tetracyclines',
        'Colistin', 'Phenicols', 'Macrolides', 'Rifamycin', 'Tigecycline',
        'Trimethoprim', 'Penicillins', 'Bla_acquired', 'Beta-lactamase inhibitor',
        'Carbapenems', 'ESBL', 'Other Classes',
    ]
    stdout_headers = []
    return full_headers, stdout_headers


def add_cli_options(parser):
    """
    command-line options for kosc__amr module
    """
    module_name = os.path.basename(__file__)[:-3]
    group = parser.add_argument_group(f'{module_name} module')
    return group


def check_cli_options(args):
    if not shutil.which('amrfinder'):
        sys.exit('Error: AMRFinderPlus is not installed or not in PATH.')
    if not getattr(args, 'no_hamronization', False):
        try:
            import hAMRonization
        except ImportError:
            sys.exit('Error: hAMRonization Python package is not installed.')


def check_external_programs():
    """
    Ensure the required external programs are available.
    """
    programs = ['amrfinder']
    if not shutil.which('amrfinder'):
        sys.exit('Error: could not find AMRFinderPlus executable.')
    if shutil.which('hamronize'):
        programs.append('hamronize')
    return programs


@lru_cache(maxsize=1)
def get_amrfinder_versions():
    """Return (software_version, database_version) from `amrfinder -V`."""
    out = subprocess.run(
        ["amrfinder", "-V"],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, check=True,
    ).stdout

    sw = re.search(r'Software version:\s*(\S+)', out)
    db = re.search(r'Database version:\s*(\S+)', out)
    if not (sw and db):
        raise RuntimeError(f"Could not parse AMRFinderPlus versions from:\n{out}")
    return sw.group(1), db.group(1)



def categorize_class(class_name):
    """
    Categorize the AMR determinants into classes.
    Betalactams are classified based on sub-class
    """
    class_name = class_name.strip().upper()
    class_map = {
        'AMINOGLYCOSIDE': 'Aminoglycosides',
        'BETA-LACTAM': 'Penicillins',
        'CARBAPENEM': 'Carbapenems',
        'CEPHALOSPORIN': 'ESBL',
        'MACROLIDE': 'Macrolides',
        'PHENICOL': 'Phenicols',
        'QUINOLONE': 'Quinolones',
        'FLUOROQUINOLONE': 'Quinolones',
        'SULFONAMIDE': 'Sulfonamides',
        'TETRACYCLINE': 'Tetracyclines',
        'TIGECYCLINE': 'Tigecycline',
        'TRIMETHOPRIM': 'Trimethoprim',
        'RIFAMYCIN': 'Rifamycin',
        'COLISTIN': 'Colistin',
        'CEPHALOTHIN': 'Penicillins',
        'FOSFOMYCIN': 'Fosfomycin',
        'INHIBITOR': 'Beta-lactamase inhibitor',
    }
    return class_map.get(class_name, 'Other Classes')


def run_amrfinder(input_fasta, organism):
    """
    run AMRFinder
    """
    command = [
        "amrfinder",
        "-n", input_fasta,
        "-O", organism,
        "--plus",
        "-q"
    ]
    try:
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error occurred running AMRFinderPlus: {e}")
        return None


def filter_amr_only(raw_amrfinder_output):
    """
    Keep the header plus rows where Type == AMR, dropping the STRESS, and VIRULENCE factors .
    """
    lines = [line for line in raw_amrfinder_output.splitlines() if line.strip()]
    if not lines:
        return ""

    header = lines[0].split("\t")
    if "Type" not in header:
        return raw_amrfinder_output

    type_idx = header.index("Type")
    kept = [lines[0]]
    for line in lines[1:]:
        columns = line.split("\t")
        if len(columns) > type_idx and columns[type_idx].strip() == "AMR":
            kept.append(line)

    if len(kept) == 1:
        return ""
    return "\n".join(kept) + "\n"


def run_hamronize(raw_amrfinder_output, input_fasta_path, input_name=None):
    """
    Run hAMRonization (Python API) on the raw AMRFinderPlus output.
    """
    if not raw_amrfinder_output or not raw_amrfinder_output.strip():
        return ""

    # Filter out non-AMR rows (stress response, virulence, etc.) before hAMRonization
    raw_amrfinder_output = filter_amr_only(raw_amrfinder_output)
    if not raw_amrfinder_output:
        return ""

    import hAMRonization

    sw_ver, db_ver = get_amrfinder_versions()
    input_file_name = input_name or Path(input_fasta_path).name

    metadata = {
        "analysis_software_version": sw_ver,
        "reference_database_version": db_ver,
        "input_file_name": input_file_name,
    }

    # hAMRonization.parse needs a file path, so write the AMRFinderPlus output to a temp file
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.tsv') as temp_amr_tsv:
        temp_amr_tsv.write(raw_amrfinder_output)
        temp_amr_path = temp_amr_tsv.name
    temp_out_path = temp_amr_path + '.hamronized.tsv'

    try:
        parsed_report = hAMRonization.parse(temp_amr_path, metadata, "amrfinderplus")
        parsed_report.write(output_location=temp_out_path)
        if os.path.exists(temp_out_path):
            with open(temp_out_path) as f:
                return f.read()
        return ""
    except Exception as e:
        print(f"Error occurred running hAMRonization: {e}")
        return ""
    finally:
        for p in (temp_amr_path, temp_out_path):
            if os.path.exists(p):
                os.remove(p)


def parse_amrfinder_results(output, split_map=None):
    full_headers, _ = get_headers()
    results = {header: '-' for header in full_headers if header != 'Other Classes'}
    results['Other Classes'] = []

    lines = [line.strip() for line in output.splitlines() if line.strip()]
    if not lines:
        print("AMRFinder output is empty")
        return results

    file_headers = lines[0].split("\t")
    for line in lines[1:]:
        columns = line.strip().split("\t")
        result_dict = dict(zip(file_headers, columns))
        
        if result_dict.get('Type') != 'AMR':
            continue

        raw_class = result_dict.get('Class', 'Other Classes').strip()
        subclass_name = result_dict.get('Subclass', '').strip()
        element_symbol = result_dict.get('Element symbol', '').strip()
        subtype = result_dict.get('Subtype', '').strip().upper()
        identity = result_dict.get('% Identity to reference', '').strip()
        coverage = result_dict.get('% Coverage of reference', '').strip()

        # annotate genes based on identity and coverage thresholds
        if subtype == 'AMR':
            if coverage and float(coverage) < 100.0:
                element_symbol += '?'
            if identity and float(identity) < 100.0:
                element_symbol += '*'

        # Split class string by "/" to handle multi-class entries (e.g. AMINOGLYCOSIDE/QUINOLONE)
        class_parts = [part.strip() for part in raw_class.split('/') if part.strip()]

        for class_part in class_parts:
            class_upper = class_part.upper()
            categories = []
            if class_upper == 'BETA-LACTAM':
                subclass_map = {
                    'BETA-LACTAM': 'Penicillins',
                    'CARBAPENEM': 'Carbapenems',
                    'CEPHALOSPORIN': 'ESBL',
                    'TANIBORBACTAM': 'Beta-lactamase inhibitor',
                    'AMOXICILLIN-CLAVULANIC ACID': 'Beta-lactamase inhibitor',
                    'PIPERACILLIN-TAZOBACTAM': 'Beta-lactamase inhibitor',
                    'TICARCILLIN-CLAVULANIC ACID': 'Beta-lactamase inhibitor',
                    'SULBACTAM-DURLOBACTAM': 'Beta-lactamase inhibitor',
                    'CEFTAZIDIME-AVIBACTAM': 'Beta-lactamase inhibitor',
                    'TAZOBACTAM': 'Beta-lactamase inhibitor'
                }

                subclass_parts = [part.strip().upper() for part in subclass_name.split('/') if part.strip()]
                for sc in subclass_parts:
                    if sc in subclass_map:
                        cat = subclass_map[sc]
                        if cat not in categories:
                            categories.append(cat)

                if not categories:
                    categories.append('Penicillins')
            else:
                categories.append(categorize_class(class_upper))

            for category in categories:
                if category == 'Penicillins' and subtype != 'POINT':
                    category = 'Bla_acquired'

                if category == 'Other Classes':
                    results['Other Classes'].append((class_part.title(), element_symbol))
                elif category in results:
                    if results[category] == '-':
                        results[category] = element_symbol
                    else:
                        existing_genes = results[category].split(',')
                        if element_symbol not in existing_genes:
                            results[category] += ',' + element_symbol

    # Normalise the output
    for category, value in results.items():
        if category == 'Other Classes':
            if value:
                grouped_classes = defaultdict(list)
                for cls, gene in value:
                    if gene not in grouped_classes[cls]:
                        grouped_classes[cls].append(gene)
                results[category] = ";".join(f"{cls}:{';'.join(genes)}" for cls, genes in grouped_classes.items())
            else:
                results[category] = '-'
        else:
            if value != '-':
                value = value.replace(" ", "").replace(";", "").replace(",", ";")
            results[category] = value

    return results


def get_results(assembly, ref_index, args, previous_results=None):
    organism = "Klebsiella_oxytoca"
    raw_output = run_amrfinder(assembly, organism)
    if raw_output:
        results = parse_amrfinder_results(raw_output)
        if not getattr(args, 'no_hamronization', False):
            # run_hamronize filters to Type == AMR internally
            strain_name = (previous_results or {}).get('strain')
            results['hamronization_output'] = run_hamronize(raw_output, assembly, strain_name)
        return results
    return {}
