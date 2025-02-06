
"""
Copyright 2025 Mary Maranga,  ebenezer foster nyarko, Kat Holt
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

    full_headers = ['Aminoglycoside', 'Carbapenemase', 'BETA-LACTAM', 'Cephalosporin','Macrolide', 'Colistin','Phenicol', 'Quinolone', 'Sulfonamide', 'Tetracycline', 'Trimethoprim', 'Other Classes']
    stdout_headers = []
    return full_headers, stdout_headers



def add_cli_options(parser):
    """
    Add command-line options for this module.
    """
    module_name = os.path.basename(__file__)[:-3]
    group = parser.add_argument_group(f'{module_name} module')
    # group.add_argument('-t', '--threads', type=int, default=8, metavar='',
                       # help="Number of threads for AMRFinderPlus (default: %(default)s).")
    group.add_argument('--plus', action='store_true', default=False,
                       help="Use the --plus option in AMRFinderPlus (default: %(default)s).")
    group.add_argument('-q', '--quiet', action='store_true', default=False,
                       help="Suppress additional AMRFinderPlus output (default: %(default)s).")
    return group


def check_cli_options(args):
    """
    Validate the command-line arguments.
    """
    # if args.threads < 1:
    #     raise ValueError("The number of threads must be at least 1.")
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
    Categorize the AMR determinants into sub-class.
    """
    class_name = class_name.strip().upper()

    class_map = {
        'AMINOGLYCOSIDE': 'Aminoglycoside',
        'BETA-LACTAM': 'BETA-LACTAM', 
        'CARBAPENEM': 'Carbapenemase',
        'CEPHALOSPORIN': 'Cephalosporin',
        'MACROLIDE': 'Macrolide',
        'PHENICOL': 'Phenicol',
        'QUINOLONE': 'Quinolone',
        'SULFONAMIDE': 'Sulfonamide',
        'TETRACYCLINE': 'Tetracycline',
        'TRIMETHOPRIM': 'Trimethoprim',
        'COLISTIN': 'Colistin'
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
        "--plus",
        "-q"
    ]
    
    try:
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e}")
        return None


# headers for the AMRFINDER output
headers = [
        'Protein id', 'Contig id', 'Start', 'Stop', 'Strand', 'Element symbol', 'Element name', 
        'Scope', 'Type', 'Subtype', 'Class', 'Subclass', 'Method', 'Target length', 'Reference sequence length', 
        '% Coverage of reference', '% Identity to reference', 'Alignment length', 'Closest reference accession', 
        'Closest reference name', 'HMM accession', 'HMM description', 'Hierarchy node'
    ]


def parse_amrfinder_results(output):
    full_headers, _ = get_headers()  
    
    results = {
        'Aminoglycoside': '-',
        'Carbapenemase': '-',
        'BETA-LACTAM': '-',
        'Cephalosporin': '-',
        'Colistin': '-',
        'Macrolide': '-',
        'Phenicol': '-',
        'Quinolone': '-',
        'Sulfonamide': '-',
        'Tetracycline': '-',
        'Trimethoprim': '-',
        'Other Classes': '-'
    }

    
    lines = output.splitlines()
    lines = [line.strip() for line in lines if line.strip()]

    if not lines:
        print("AMRFinder output is empty")
        return results  

    
    headers = lines[0].split("\t")
    
    for line in lines[1:]:
        columns = line.strip().split("\t")
        result_dict = dict(zip(headers, columns))  
    
        if result_dict.get('Type') != 'AMR':
            continue

        class_name = result_dict.get('Subclass', 'Other Classes')
        element_symbol = result_dict.get('Element symbol', '')

        
        category = categorize_class(class_name)

        if element_symbol:
            if results.get(category) == '-':
                results[category] = element_symbol
            else:
                results[category] += ',' + element_symbol

    for category in results:
        if results[category] != '-':
            results[category] = results[category].replace(" ", "").replace(";", "").replace(",", ";")

    return results



def get_results(assembly, index, previous_results, args):
    
    organism = "Escherichia"
    
    assembly_name = previous_results.assemblies[0].split('/')[-1].replace('.fasta', '')
    
    # Run AMRFinder
    raw_output = run_amrfinder(assembly, organism)

    if raw_output:
        results = parse_amrfinder_results(raw_output)
        return results
    
    return {}



# def parse_amrfinder_results(output, assembly_name):
#     full_headers, _ = get_headers() 
#     results = {}

#     # Split output into lines
#     lines = output.splitlines()

#     lines = [line.strip() for line in lines if line.strip()]

#     if not lines:
#         print("AMRFinder output is empty")
#         return results  

#     results[assembly_name] = {
#         'Aminoglycoside': '-',
#         'Carbapenemase': '-',
#         'BETA-LACTAM': '-',
#         'Macrolide': '-',
#         'Phenicol': '-',
#         'Quinolone': '-',
#         'Sulfonamide': '-',
#         'Tetracycline': '-',
#         'Trimethoprim': '-',
#         'Other Classes': '-'
#     }

#     # Loop through the lines (skip header line)
#     for line in lines[1:]:
#         columns = line.strip().split("\t")
#         result_dict = dict(zip(headers, columns))  # Map column values to headers

#         if result_dict.get('Type') != 'AMR':
#             continue

#         class_name = result_dict.get('Class', 'Other Classes')
#         element_symbol = result_dict.get('Element symbol', '')

#         # Categorize the class
#         category = categorize_class(class_name)

#         if element_symbol:
#             if results[assembly_name].get(category) == '-':
#                 results[assembly_name][category] = element_symbol
#             else:
#                 results[assembly_name][category] += element_symbol

#     for category in results[assembly_name]:
#         if results[assembly_name][category] != '-':
#             # Remove any unwanted characters like semicolons or extra spaces
#             results[assembly_name][category] = results[assembly_name][category].replace(" ", "").replace(";", "")

#     return results



# def get_results(assembly, minimap2_index, previous_results, args):
#     """
#     Execute AMRFinderPlus
#     """
#     organism = "Escherichia"
    
#     # Extract assembly name from previous_results
#     assembly_name = previous_results.assemblies[0].split('/')[-1].replace('.fasta', '')
    
#     # Run AMRFinder
#     raw_output = run_amrfinder(assembly, organism)

#     if raw_output:
#         # Pass the assembly name to parse_amrfinder_results
#         results = parse_amrfinder_results(raw_output, assembly_name)

#         # Prepare the final output format
#         final_result = {assembly_name: {}}
        
#         if assembly_name in results and results[assembly_name]:
#             for category in results[assembly_name]:
#                 # Only add categories that have elements
#                 if results[assembly_name][category]:
#                     final_result[assembly_name][category] = results[assembly_name][category]
#                 else:
#                     final_result[assembly_name][category] = '-'
#         else:
#             print(f"Error: {assembly_name} not found or empty in results")

#         print(f"Final result for {assembly_name}: {final_result}")
        
#         return final_result
#     else:
#         print(f"AMRFinder analysis failed for assembly: {assembly_name}")
#         return {}



    


