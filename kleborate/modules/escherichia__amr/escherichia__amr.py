
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
import sys
import shutil
from pathlib import Path
import pandas as pd
from collections import defaultdict
import subprocess


def description():
    return 'Antimicrobial resistance (AMR) gene detection using AMRFinderPlus'


def prerequisite_modules():
    return []


def get_headers():
    """
    Define the headers for AMRFinderPlus results.
    """
    full_headers = [
        'Aminoglycosides', 'Quinolones', 'Fosfomycin', 'Sulfonamides', 'Tetracyclines',
        'Colistin', 'Phenicols', 'Macrolides', 'Rifamycin', 'Tigecycline',
        'Trimethoprim', 'Penicillins', 'Beta-lactamase inhibitor', 'Carbapenems', 'ESBL',
        'Bla+Inhibitor', 'Carb+Inhibitor', 'ESBL+Inhibitor',
        'Other Classes'
    ]
    stdout_headers = []
    return full_headers, stdout_headers


def add_cli_options(parser):
    """
    command-line options for amr module
    """
    module_name = os.path.basename(__file__)[:-3]
    group = parser.add_argument_group(f'{module_name} module')
    group.add_argument('--plus', action='store_true', default=False,
                       help="Use the --plus option in AMRFinderPlus (default: %(default)s).")
    group.add_argument('-q', '--quiet', action='store_true', default=False,
                       help="Suppress additional AMRFinderPlus output (default: %(default)s).")
    return group


def check_cli_options(args):
    if not shutil.which('amrfinder'):
        sys.exit('Error: AMRFinderPlus is not installed or not in PATH.')


def check_external_programs():
    """
    Ensure the required external programs are available.
    """
    if not shutil.which('amrfinder'):
        sys.exit('Error: could not find AMRFinderPlus executable.')
    return ['amrfinder']


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


combined_category = {
    'Carbapenems': 'Carb+Inhibitor',
    'ESBL': 'ESBL+Inhibitor',
    'Penicillins': 'Bla+Inhibitor',
}


def run_amrfinder(input_fasta, organism):
    """
    run AMRFinder
    Parameters:
        input_fasta (str): Path to the input FASTA file.
        organ organism (str): The organism name.
    Returns:
        str: output from the AMRFinderPlus.
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
        print(f"Error occurred: {e}")
        return None


# header columns for the AMRFinder
headers = [
    'Protein id', 'Contig id', 'Start', 'Stop', 'Strand', 'Element symbol', 'Element name', 
    'Scope', 'Type', 'Subtype', 'Class', 'Subclass', 'Method', 'Target length', 'Reference sequence length', 
    '% Coverage of reference', '% Identity to reference', 'Alignment length', 'Closest reference accession', 
    'Closest reference name', 'HMM accession', 'HMM description', 'Hierarchy node'
]


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

        if element_symbol == 'blaEC':
            continue

        # Split class string by "/" to handle multi-class entries (e.g. AMINOGLYCOSIDE/QUINOLONE)
        class_parts = [part.strip() for part in raw_class.split('/') if part.strip()]

        for class_part in class_parts:
            class_upper = class_part.upper()

            categories = []
            if class_upper == 'BETA-LACTAM':
                subclass_map = {
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
                
                # Split subclass string by "/" to handle combined subclasses (e.g. CARBAPENEM/TANIBORBACTAM)
                subclass_parts = [part.strip().upper() for part in subclass_name.split('/') if part.strip()]
                
                for sc in subclass_parts:
                    if sc in subclass_map:
                        cat = subclass_map[sc]
                        if cat not in categories:
                            categories.append(cat)
                
                # Default to Penicillins
                if not categories:
                    categories.append('Penicillins')

                
                # combination column (Carb+Inhibitor / ESBL+Inhibitor / Bla+Inhibitor)
                if 'Beta-lactamase inhibitor' in categories and len(categories) > 1:
                    other_cats = [c for c in categories if c != 'Beta-lactamase inhibitor']
                    categories = [combined_category.get(oc, oc) for oc in other_cats]
            else:
                categories.append(categorize_class(class_upper))

            for category in categories:
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


def get_results(assembly, ref_index, previous_results, args):
    organism = "Escherichia"
    raw_output = run_amrfinder(assembly, organism)
    if raw_output:
        results = parse_amrfinder_results(raw_output)
        return results
    return {}
