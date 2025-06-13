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

from ...shared.alignment import align_query_to_ref, cull_redundant_hits


def minimap_pathovar(assembly, minimap2_index, ref_file, min_identity, min_coverage):
    """
    Aligns assembled genomes to the virulence alleles and classifies the pathotype.

    Parameters:
    - assembly: Assembly in FASTA format.
    - ref_file: Virulence factors in FASTA format.
    - minimap2_index: Path to the assembly's minimap2 index for faster alignment (optional).
    - min_coverage: Minimum query coverage for alignment.
    - min_identity: Minimum identity percentage for alignment.

    Returns:
    - pathotype and virulence markers.
    """

    alignment_hits = align_query_to_ref(
        ref_file,
        assembly,
        ref_index=minimap2_index,
        min_identity=min_identity,
        min_query_coverage=min_coverage
    )
    alignment_hits = cull_redundant_hits(alignment_hits)

    # Load virulence factors mapping file
    virulence_factors_map = load_virulence_factors(ref_file)
    
    # Identify virulence factors from the genome
    virulence_factors, virulence_markers = identify_virulence_factors(
        alignment_hits, 
        virulence_factors_map
    )
    
    # Classify the pathotype
    pathovar = classify_pathovar(virulence_factors)

    return pathovar, virulence_markers



# Define the virulence factors map
virulence_factors_map = {
    'ltcA': {'name': 'LT', 'headers': []},
    'sta1': {'name': 'ST', 'headers': []},
    'stx1A': {'name': 'Stx1', 'headers': []},
    'stx1B': {'name': 'Stx1', 'headers': []},
    'stx2A': {'name': 'Stx2', 'headers': []},
    'stx2B': {'name': 'Stx2', 'headers': []},
    'eae': {'name': 'eae', 'headers': []},
    'ipaH': {'name': 'ipaH', 'headers': []},
 }


def load_virulence_factors(fasta_file_path):
    """
    Loads virulence factor headers from a FASTA file.

    Parameters:
    - fasta_file_path (str): Path to the FASTA file.

    Returns:
    - Dictionary with virulence factors mapped to the headers.
    """
    with open(fasta_file_path, 'r') as fasta_file:
        for line in fasta_file:
            if line.startswith('>'):
                header = line[1:].strip()
                parts = header.split(':')
                virulence_factor = parts[0]

                # Check if the virulence factor exists in the map
                if virulence_factor in virulence_factors_map:
                    virulence_factors_map[virulence_factor]['headers'].append(header)
                elif virulence_factor.startswith('ipaH') and 'ipaH' in virulence_factors_map:
                    # Handle cases like 'ipaH9.8'
                    virulence_factors_map['ipaH']['headers'].append(header)

    return virulence_factors_map


def identify_virulence_factors(alignment_hits, virulence_factors_map):
    """
    Identifies presence of virulence factors based on alignment hits.

    Parameters:
    - alignment_hits (list): List of alignment hit objects with query_name attributes.
    - virulence_factors_map (dict): Dictionary with virulence factors mapped to headers.

    Returns:
    - virulence_factors dictionary
    - Dictionary mapping detected markers to virulence factor names.
    """
    virulence_factors = {
        'ipaH': '-',
        'ST': '-',
        'LT': '-',
        'Stx1': '-',
        'Stx2': '-',
        'eae': '-',
    }
    virulence_markers = {}

    for hit in alignment_hits:
        header = hit.query_name
        for key, value in virulence_factors_map.items():
            if key in header:
                virulence_factors[value['name']] = '+'
                
                # For Stx1 and Stx2 
                if value['name'] in ['Stx1', 'Stx2']:
                    # take the first part
                    first_part = header.split(':')[0]
                    # If header ends with a variant append it
                    if header.endswith(':d'):
                        combined_name = f"{first_part}d"
                    else:
                        combined_name = first_part
                else:
                    combined_name = header.split(':')[0]
                
                factor_name = value['name']
                if factor_name not in virulence_markers:
                    virulence_markers[factor_name] = []
                virulence_markers[factor_name].append(combined_name)

    return virulence_factors, virulence_markers


def classify_pathovar(virulence_factors):
    """
    Classifies E. coli pathotype based on detected virulence factors.

    Parameters:
    - virulence_factors (dict): Dictionary with detected virulence factors.

    Returns:
    - Pathotype
    """
    pathovar = ''

    if virulence_factors['ipaH'] == '+':
        pathovar = 'EIEC'
    if virulence_factors['ST'] == '+' or virulence_factors['LT'] == '+':
        # Ensure ipaH is not in combination with ST or LT
        pathovar = 'ERROR' if pathovar else 'ETEC'
    if virulence_factors['Stx1'] == '+' or virulence_factors['Stx2'] == '+':
        # Stx factors classify as STEC unless eae is also present, which makes it EHEC
        if virulence_factors['eae'] == '+':
            pathovar = pathovar + ('/' if len(pathovar) else '') + 'EHEC'
        else:
            pathovar = pathovar + ('/' if len(pathovar) else '') + 'STEC'
    elif virulence_factors['eae'] == '+':
        # eae without Stx is classified as EPEC
        pathovar = pathovar + ('/' if len(pathovar) else '') + 'EPEC'

    if len(pathovar):
        pathovar = pathovar
    else:
        pathovar = '-'

    return pathovar  


    