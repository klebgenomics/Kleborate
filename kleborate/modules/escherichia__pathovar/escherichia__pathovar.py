"""
Copyright 2024 Mary Maranga (gathonimaranga@gmail.com), Kat Holt
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
import ast

from .pathovar import minimap_pathovar
from .shigella_classify import classify_shigella


def description():
    return 'Pathotyping of E. coli genomes'


def prerequisite_modules():
    return []


def get_headers():
    full_headers = ['Pathotype','Stx1', 'Stx2','Shiga toxin','ST', 'LT', 'eae', 'ipaH']
    stdout_headers = []
    return full_headers, stdout_headers


def add_cli_options(parser):
    module_name = os.path.basename(__file__)[:-3]
    group = parser.add_argument_group(f'{module_name} module')
    group.add_argument('--escherichia__pathovar_min_identity', type=float, default=90.0,
                       help='Minimum alignment percent identity for detecting virulence factors')
    group.add_argument('--escherichia__pathovar_min_coverage', type=float, default=80.0,
                       help='Minimum alignment percent coverage for detecting virulence factors')
    return 



def check_cli_options(args):
    if args.escherichia__pathovar_min_identity <= 50.0 or args.escherichia__pathovar_min_identity >= 100.0:
        sys.exit('Error: --escherichia__pathovar_min_identity must be between 50.0 and 100.0')
    if args.escherichia__pathovar_min_coverage <= 50.0 or args.escherichia__pathovar_min_coverage >= 100.0:
        sys.exit('Error: --escherichia__pathovar_min_coverage must be between 50.0 and 100.0')
    

def check_external_programs():
    if not shutil.which('minimap2'):
        sys.exit('Error: could not find minimap2')
    return ['minimap2']


def data_dir():
    return pathlib.Path(__file__).parents[0] / 'data'


def get_results(assembly, minimap2_index, args, previous_results):
    full_headers, _ = get_headers()
    
    # Load the virulence factors database
    ref_file = data_dir() / 'virulence_ecoli.fsa'
    
    # Load the Shigella reference file and serotype marker dictionary
    ShigellaRef = data_dir() / 'ShigellaRef5.fasta'
    with open(data_dir() / 'shigella_serotype_markers.txt', 'r') as file:
        file_content = file.read()
        shigella_serotype_markers = ast.literal_eval(file_content)

    # pathotyping
    pathovar, virulence_markers = minimap_pathovar(
        assembly,
        minimap2_index,
        ref_file,
        args.escherichia__pathovar_min_identity,
        args.escherichia__pathovar_min_coverage
    )

      # print(virulence_markers)

    # Initialize the result dictionary with headers as keys
    result_dict = {header: '-' for header in full_headers}

    # Set the virulence markers
    for marker, marker_hits in virulence_markers.items():
        if marker in result_dict:
            if isinstance(marker_hits, list):
                result_dict[marker] = ";".join(marker_hits) if marker_hits else '-'
            else:  # If it's a string (e.g., 'Shiga toxin')
                result_dict[marker] = marker_hits

    # Call classify_shigella
    Serotype = classify_shigella(
        assembly,
        minimap2_index,
        ShigellaRef,
        args.escherichia__pathovar_min_identity,
        args.escherichia__pathovar_min_coverage,
        shigella_serotype_markers
    )

    # Combine the pathotype and Shigella Serotype results
    combined_pathovar = []
    if pathovar != '-':
        combined_pathovar.append(pathovar)
    if Serotype != '-':
        combined_pathovar.append(Serotype)

    # Set the combined Pathovar value in the result dictionary
    result_dict['Pathotype'] = ', '.join(combined_pathovar) if combined_pathovar else '-'

    # Return the result dictionary
    return result_dict
    

# def get_results(assembly, minimap2_index, args, previous_results):
    
#     full_headers, _ = get_headers()
    
#     # Load the virulence factors database
#     ref_file = data_dir() / 'virulence_ecoli.fsa'
    
#     # Load the Shigella reference file and serotype marker dictionary
#     ShigellaRef = data_dir() / 'ShigellaRef5.fasta'
#     with open(data_dir() / 'shigella_serotype_markers.txt', 'r') as file:
#         file_content = file.read()
#         shigella_serotype_markers = ast.literal_eval(file_content)

#     # pathotyping
#     pathovar, virulence_markers = minimap_pathovar(
#         assembly,
#         minimap2_index,
#         ref_file,
#         args.escherichia__pathovar_min_identity,
#         args.escherichia__pathovar_min_coverage
#     )

#     print(virulence_markers)

#     # Initialize the result dictionary with headers as keys
#     result_dict = {header: '-' for header in full_headers}

#     # Set the virulence markers
#     for marker, marker_hits in virulence_markers.items():
#         if marker in result_dict:
#             result_dict[marker] = ";".join(marker_hits) if marker_hits else '-'

#     # Call classify_shigella
#     Serotype = classify_shigella(
#         assembly,
#         minimap2_index,
#         ShigellaRef,
#         args.escherichia__pathovar_min_identity,
#         args.escherichia__pathovar_min_coverage,
#         shigella_serotype_markers
#     )

#     # Combine the pathotype and Shigella Serotype results
#     combined_pathovar = []
#     if pathovar != '-':
#         combined_pathovar.append(pathovar)
#     if Serotype != '-':
#         combined_pathovar.append(Serotype)

#     # Set the combined Pathovar value in the result dictionary
#     result_dict['Pathotype'] = ', '.join(combined_pathovar) if combined_pathovar else '-'

#     # Return the result dictionary
#     return result_dict


# def get_results(assembly, minimap2_index, args, previous_results):
#     full_headers, _ = get_headers() 
#     # Load the virulence factors file
#     ref_file = data_dir() / 'virulence_ecoli.fsa'
    
#     # Call minimap_pathovar to get the pathovar and virulence factors data
#     pathovar, virulence_markers = minimap_pathovar(
#         assembly,
#         minimap2_index,
#         ref_file,
#         args.escherichia__pathovar_min_identity,
#         args.escherichia__pathovar_min_coverage
#     )

#     # Initialize the result dictionary with headers as keys and default value as None
#     result_dict = {header: None for header in full_headers}

#     # Set the pathovar value
#     result_dict['Pathovar'] = pathovar

#     # Set the virulence factors if they are present in virulence_markers
#     for marker in virulence_markers:
#         if marker in result_dict:
#             result_dict[marker] = ", ".join(virulence_markers[marker])

#     # Replace any None values with '-'
#     result_dict = {key: (value if value is not None else '-') for key, value in result_dict.items()}

#     # Return the result dictionary
#     return result_dict


