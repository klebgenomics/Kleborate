
"""Copyright 2025 Mary Maranga
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

    full_headers = ['Aminoglycoside', 'Fluoroquinolone','Fosfomycin','Sulfonamide','Tetracycline',
                    'Glycopeptide','Colistin','Phenicol','Macrolide','Rifamycin',
                    'Trimethoprim','Betalactam','Carbapenem','Cephalosporin',
                    'Methicillin', 'Other Classes']               

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
    Betalacm are classified based on sub-class
    """
    class_name = class_name.strip().upper()

    class_map = {
        'AMINOGLYCOSIDE': 'Aminoglycoside',
        'AMINOGLYCOSIDE/QUINOLONE' : 'Aminoglycoside/Fluoroquinolone', # genes such as aac(6')-Ib-cr5 
        'BETA-LACTAM': 'Betalactam', 
        'CARBAPENEM': 'Carbapenem',
        'CEPHALOSPORIN': 'Cephalosporin',
        'MACROLIDE': 'Macrolide',
        'PHENICOL': 'Phenicol',
        'QUINOLONE': 'Fluoroquinolone',
        'FLUOROQUINOLONE': 'Fluoroquinolone',
        'SULFONAMIDE': 'Sulfonamide',
        'TETRACYCLINE': 'Tetracycline',
        # 'TIGECYCLINE': 'Tigecycline',
        'TRIMETHOPRIM': 'Trimethoprim',
        'RIFAMYCIN':'Rifamycin',
        'COLISTIN': 'Colistin',
        'CEPHALOTHIN': 'Betalactam',
        'QUINOLONE/TRICLOSAN':'Fluoroquinolone',
        'METHICILLIN': 'Methicillin',
        'FOSFOMYCIN': 'Fosfomycin',
        'GLYCOPEPTIDE':'Glycopeptide',
        'PHENICOL/QUINOLONE':'Phenicol/Quinolone',
        'PHENICOL/OXAZOLIDINONE':'Phenicol/Oxazolidinone',
        'MACROLIDE/PHENICOL':'Macrolide/Phenicol',
        'MACROLIDE/PLEUROMUTILIN':'Macrolide/Pleuromutilin'
    }

    # Return the mapped class or 'Other Classes'
    return class_map.get(class_name, 'Other Classes')


def run_amrfinder(input_fasta, organism):
    """
    run AMRFinder
    
    Parameters:
        input_fasta (str): Path to the input FASTA file.
        organism (str): The organism name.
    
    Returns:
        str: output from the AMRFinderPlus.
    """
    command = [
        "amrfinder",
        "-n", input_fasta,
        "-O", organism,
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

def parse_amrfinder_results(output, split_map):
    full_headers, _ = get_headers()

    results = {
        'Aminoglycoside': '-',
        'Fluoroquinolone': '-',
        'Fosfomycin': '-',
        'Sulfonamide': '-',
        'Tetracycline': '-',
        'Glycopeptide': '-',
        'Colistin': '-',
        'Phenicol': '-',
        'Macrolide': '-',
        # 'Tigecycline': '-',
        'Rifamycin': '-',
        'Trimethoprim': '-',
        'Betalactam': '-',  
        'Carbapenem': '-',  
        'Cephalosporin': '-',
        'Methicillin': '-',
        'Other Classes': []   # Change from '-' to an empty list
    }

    genes_to_remove = ['blaEC-5', 'pmrB', 'glpT', 'mdf(A)']
    lines = [line.strip() for line in output.splitlines() if line.strip()]
    if not lines:
        print("AMRFinder output is empty")
        return results

    headers = lines[0].split("\t")

    for line in lines[1:]:
        columns = line.strip().split("\t")
        result_dict = dict(zip(headers, columns))
        if result_dict.get('Type') != 'AMR':
            continue

        class_name = result_dict.get('Class', 'Other Classes')
        subclass_name = result_dict.get('Subclass', '')
        element_symbol = result_dict.get('Element symbol', '')

        if element_symbol in genes_to_remove:
            continue

        if class_name == 'BETA-LACTAM':
            category = 'Betalactam'
        else:
            category = categorize_class(class_name)

        if category in split_map:
            for cat in split_map[category]:
                if results[cat] == '-':
                    results[cat] = element_symbol
                else:
                    results[cat] += ',' + element_symbol
        elif category == 'Other Classes':
            # Store (class_name, element_symbol) pairs for later formatting
            results['Other Classes'].append((class_name.title(), element_symbol))
        else:
            if category not in results:
                continue  
            if results[category] == '-':
                results[category] = element_symbol
            else:
                results[category] += ',' + element_symbol

    # Normalise the output
    for category, value in results.items():
        if category == 'Other Classes':
            if value:
                # Format as Class:gene;Class:gene
                results[category] = ";".join(f"{cls}:{gene}" for cls, gene in value)
            else:
                results[category] = '-'
        else:
            if value != '-':
                value = value.replace(" ", "").replace(";", "").replace(",", ";")
            results[category] = value

    return results



def get_results(assembly, index, previous_results, args):
    organism = "Escherichia"
    raw_output = run_amrfinder(assembly, organism)
    
    if raw_output:

        # genes reported in multiple classes
        split_map = {
            'Aminoglycoside/Fluoroquinolone': ('Aminoglycoside', 'Fluoroquinolone'),
            'Phenicol/Quinolone': ('Phenicol', 'Fluoroquinolone'), 
            'QUINOLONE/TRICLOSAN':('Fluoroquinolone'),
            'Phenicol/Oxazolidinone': ('Phenicol'),
            'Macrolide/Phenicol': ('Macrolide', 'Phenicol'),
            'Macrolide/Pleuromutilin': ('Macrolide',),
        }
        results = parse_amrfinder_results(raw_output, split_map)
        return results
    return {}

