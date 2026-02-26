
"""
Copyright 2025 Mary Maranga
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

    amrfinder_headers = ['Aminoglycoside', 'Fluoroquinolone','Fosfomycin',
                        'Sulfonamide','Tetracycline',
                        'Glycopeptide','Colistin','Phenicol','Macrolide','Rifamycin',
                        'Trimethoprim','Betalactam','Carbapenem','Cephalosporin',
                        'Methicillin', 'Other Classes'] 

    hamronize_headers=['input_file_name', 'gene_symbol', 'gene_name', 'reference_database_name', 
                        'reference_database_version', 'reference_accession', 'analysis_software_name', 
                        'analysis_software_version', 'genetic_variation_type', 'antimicrobial_agent', 
                        'coverage_percentage', 'coverage_depth', 'coverage_ratio', 'drug_class', 
                        'input_gene_length', 'input_gene_start', 'input_gene_stop', 'input_protein_length', 
                        'input_protein_start', 'input_protein_stop', 'input_sequence_id', 'nucleotide_mutation', 
                        'nucleotide_mutation_interpretation', 'predicted_phenotype', 
                        'predicted_phenotype_confidence_level', 'amino_acid_mutation', 
                        'amino_acid_mutation_interpretation', 'reference_gene_length', 
                        'reference_gene_start', 'reference_gene_stop', 'reference_protein_length', 
                        'reference_protein_start', 'reference_protein_stop', 'resistance_mechanism', 
                        'strand_orientation', 'sequence_identity']  

    full_headers = amrfinder_headers + hamronize_headers       

    stdout_headers = []
    return full_headers, stdout_headers

# def get_headers():
#     """
#     Define the headers for AMRFinderPlus results.
#     """

#     full_headers = ['Aminoglycoside', 'Fluoroquinolone','Fosfomycin','Sulfonamide','Tetracycline',
#                     'Glycopeptide','Colistin','Phenicol','Macrolide','Rifamycin',
#                     'Trimethoprim','Betalactam','Carbapenem','Cephalosporin',
#                     'Methicillin', 'Other Classes']               

#     stdout_headers = []
#     return full_headers, stdout_headers



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
        'AMINOGLYCOSIDE/QUINOLONE' : 'Aminoglycoside/Fluoroquinolone',
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


def run_hamronize(raw_amrfinder_output: str, input_fasta_name: str):
    """
    Run hamronize on the raw AMRFinderPlus output and return parsed dicts.

    Parameters:
        raw_amrfinder_output (str): stdout from run_amrfinder(...)
        input_fasta_name (str): label for this input (used in --input_file_name)

    Returns:
        List[Dict[str, str]]: one dict per hamronize row
    """
    if not raw_amrfinder_output:
        return []

    input_tsv = Path("amrfinder_input.tsv")

    try:
        # Write AMRFinderPlus raw TSV to a temporary file for hamronize to read
        input_tsv.write_text(raw_amrfinder_output, encoding='utf-8')

        command = [
            "hamronize",
            "amrfinderplus",
            str(input_tsv),
            "--analysis_software_version", "4.0.23",
            "--reference_database_version", "2025-06-03.1",
            "--format", "tsv",
            "--input_file_name", input_fasta_name
        ]

        result = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
            text=True
        )

        hamronize_stdout = result.stdout
        parsed_rows = parse_hamronize_results(hamronize_stdout)
        return parsed_rows

    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running hamronize: {e}")
        print(f"Stderr: {e.stderr}")
        return []

    except FileNotFoundError:
        print("Error: hamronize command not found. Please ensure it is installed and in your PATH.")
        return []

    finally:
        if input_tsv.exists():
            try:
                os.remove(input_tsv)
            except OSError:
                pass


def parse_hamronize_results(hamronize_output: str):
    """
    Parse hamronize TSV stdout into a list of dicts
    (one dict per AMR determinant / gene_symbol row).

    Parameters:
        hamronize_output (str): TSV text from hamronize stdout.
    """
    if not hamronize_output:
        return []

    # Drop empty lines
    lines = [ln.strip() for ln in hamronize_output.splitlines() if ln.strip()]
    if not lines:
        return []

    # First line = header
    header_cols = lines[0].split('\t')

    parsed_rows = []
    for ln in lines[1:]:
        cols = ln.split('\t')
        row_dict = {
            header_cols[i]: cols[i] if i < len(cols) else ''
            for i in range(len(header_cols))
        }
        parsed_rows.append(row_dict)

    return parsed_rows
    
# def run_hamronize(raw_amrfinder_output, input_fasta_name):
#     """
#     Runs hamronize on the raw output from AMRFinderPlus.

#     Parameters:
#         raw_amrfinder_output (str): The stdout from the run_amrfinder function.
#         input_fasta_name (str): Name of the input FASTA file for Hamronize metadata.

#     Returns:
#         str: The harmonized TSV output from hamronize
#     """
#     if not raw_amrfinder_output:
#         return None
    
#     # Define temporary file names
#     input_tsv = Path("amrfinder_input.tsv")

#     try:
#         # Write the raw AMRFinderPlus output to a temporary file
#         input_tsv.write_text(raw_amrfinder_output, encoding='utf-8')

#         # Define the hamronize command
#         command = [
#             "hamronize", 
#             "amrfinderplus",
#             str(input_tsv),
#             "--analysis_software_version", "4.0.23", 
#             "--reference_database_version", "2025-06-03.1", 
#             "--format", "tsv",
#             "--input_file_name"
#         ]
        
#         # Execute hamronize
#         result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True)
        
#         # The hamronize output is written to stdout by default, capture it
#         return result.stdout
        
#     except subprocess.CalledProcessError as e:
#         print(f"Error occurred while running hamronize: {e}")
#         print(f"Stderr: {e.stderr}")
#         return None
#     except FileNotFoundError:
#         print("Error: hamronize command not found. Please ensure it is installed and in your PATH.")
#         return None
#     finally:
        
#         if input_tsv.exists():
#             os.remove(input_tsv)
#         if output_tsv.exists():
#             os.remove(output_tsv) 


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
        'Other Classes': []
    }
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
        
        if class_name == 'BETA-LACTAM':
            subclass_upper = subclass_name.strip().upper()
            subclass_map = {
                'CARBAPENEM': 'Carbapenem',
                'CEPHALOSPORIN': 'Cephalosporin',
                'METHICILLIN': 'Methicillin'
            }
            if subclass_upper in subclass_map:
                category = subclass_map[subclass_upper]
            else:
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
    
    results = {}
    hamronize_output = None
    
    if raw_output:
        split_map = {
            'Aminoglycoside/Fluoroquinolone': ('Aminoglycoside', 'Fluoroquinolone'),
            'Phenicol/Quinolone': ('Phenicol', 'Fluoroquinolone'), 
            'QUINOLONE/TRICLOSAN':('Fluoroquinolone',),
            'Phenicol/Oxazolidinone': ('Phenicol',),
            'Macrolide/Phenicol': ('Macrolide', 'Phenicol'),
            'Macrolide/Pleuromutilin': ('Macrolide',),
        }
        results = parse_amrfinder_results(raw_output, split_map)
        
        input_fasta_name = Path(assembly).name
        hamronize_output = run_hamronize(raw_output, input_fasta_name)

    results['hamronize_tsv_output'] = hamronize_output if hamronize_output else ''
    
    return results


