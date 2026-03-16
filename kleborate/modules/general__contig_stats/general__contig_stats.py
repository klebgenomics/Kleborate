"""
Copyright 2023 Kat Holt, Ryan Wick (rrwick@gmail.com), Mary Maranga (gathonimaranga@gmail.com)
https://github.com/klebgenomics/KleborateModular/

This file is part of Kleborate. Kleborate is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Kleborate is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Kleborate. If
not, see <https://www.gnu.org/licenses/>.
"""

import collections
import json
import pathlib
from pathlib import Path
import ast

from ...shared.misc import load_fasta


def description():
    return 'basic stats on the assembly\'s contigs'


def prerequisite_modules():
    return ['enterobacterales__species']


def get_headers():
    full_headers = ['contig_count', 'N50', 'largest_contig', 'total_size', 'GC_content','ambiguous_bases',
                    'QC_warnings']
    stdout_headers = ['N50']
    return full_headers, stdout_headers


def add_cli_options(parser):
    pass


def check_cli_options(args):
    pass


def check_external_programs():
    return []


def data_dir():
    return pathlib.Path(__file__).parents[0] / 'data'

def get_results(assembly, minimap2_index, args, previous_results):
    species_file = data_dir() / 'species_specification.txt'
    species_specification_dict = load_species_specifications(species_file)
    species = previous_results['enterobacterales__species__species']
    
    # Calculate stats
    contig_count, n50, longest_contig, total_size, ambig, gc = get_contig_stats(assembly)
    
    QC_warnings = get_qc_warnings(
        total_size, 
        n50, 
        contig_count, 
        gc, 
        ambig, 
        species, 
        species_specification_dict
    )
    
    return {
        'contig_count': str(contig_count),
        'N50': str(n50),
        'largest_contig': str(longest_contig),
        'total_size': str(total_size),
        'GC_content': str(gc),
        'ambiguous_bases': ambig,
        'QC_warnings': QC_warnings
    }



def get_contig_stats(assembly):
    fasta = load_fasta(assembly)

    base_counts = collections.defaultdict(int)
    total_gc = 0
    total_acgt = 0

    for _, seq in fasta:
        for b in seq.upper():
            base_counts[b] += 1
            if b in ("A", "C", "G", "T"):
                total_acgt += 1
                if b in ("G", "C"):
                    total_gc += 1

    # ambiguous bases = anything not A/C/G/T
    base_counts.pop('A', None)
    base_counts.pop('C', None)
    base_counts.pop('G', None)
    base_counts.pop('T', None)
    ambiguous_base_count = sum(base_counts.values())

    if ambiguous_base_count:
        ambiguous_bases = f"yes ({ambiguous_base_count})"
    else:
        ambiguous_bases = "no"

    contig_lengths = sorted(len(seq) for _, seq in fasta)

    contig_lengths = sorted([len(x[1]) for x in fasta])
    if not contig_lengths:
        return 0, 0, 0, 0, 'no'
    longest = contig_lengths[-1]

    total_size = sum(contig_lengths)

    # N50
    half_total_length = total_size / 2
    total_so_far = 0
    segment_lengths = contig_lengths[::-1]
    N50 = 0
    for length in segment_lengths:
        total_so_far += length
        if total_so_far >= half_total_length:
            N50 = length
            break


    # GC content as percentage of A/C/G/T bases
    if total_acgt > 0:
        gc_content = round((total_gc / total_acgt) * 100.0, 2)
    else:
        gc_content = 0.0

    # number_of_contigs, N50, longest, total_size, ambiguous_bases, GC_Content
    return len(contig_lengths), N50, longest, total_size, ambiguous_bases, gc_content


def load_species_specifications(file_path):
    """Loads the species thresholds from a JSON-formatted text file."""
    with open(file_path, 'r') as file:
        return json.load(file)



def get_qc_warnings(total_size, N50, contig_count, gc_content, ambiguous_bases, species, species_spec_dict):
    warnings = []
    
    # Logic to handle subspecies/versions using startswith
    spec = None
    for ref_species, ref_data in species_spec_dict.items():
        if species.startswith(ref_species):
            spec = ref_data
            break
    
    if spec is None:
        return '-' # No matching species found in dictionary

    # Genome Size Check (Min and Max)
    min_size = spec.get('min_genome_size')
    max_size = spec.get('max_genome_size')
    if min_size is not None and total_size < min_size:
        warnings.append('total_size')
    if max_size is not None and total_size > max_size:
        warnings.append('total_size')

    # N50 Check (Min only)
    min_n50 = spec.get('min_N50')
    if min_n50 is not None and N50 < min_n50:
        warnings.append('N50')

    # Contig Count Check (Min and Max)
    # min_contigs = spec.get('min_no_of_contigs')
    # max_contigs = spec.get('max_no_of_contigs')
    # if min_contigs is not None and contig_count < min_contigs:
    #     warnings.append('contig_count')
    # if max_contigs is not None and contig_count > max_contigs:
    #     warnings.append('contig_count')

    # GC Content Check (Min and Max)
    min_gc = spec.get('min_GC_Content')
    max_gc = spec.get('max_GC_Content')
    if min_gc is not None and gc_content < min_gc:
        warnings.append('GC_content')
    if max_gc is not None and gc_content > max_gc:
        warnings.append('GC_content')

    # Ambiguous Bases Check
    if 'yes' in ambiguous_bases:
        warnings.append('ambiguous_bases')

    return ','.join(warnings) if warnings else '-'


# def get_qc_warnings(total_size, N50, ambiguous_bases, species, species_specification_dict):
#     warnings = []
#     if species in species_specification_dict:
#         species_spec = species_specification_dict[species]
#         min_size, max_size = species_spec['min_genome_size'], species_spec['max_genome_size']
#         if total_size < min_size:
#             warnings.append('total_size')
#         elif total_size > max_size:
#             warnings.append('total_size')
#     else:
#         return '-'  # Skip QC for species not in the dictionary

#     if N50 < 10000:
#         warnings.append('N50')
#     if 'yes' in ambiguous_bases:
#         warnings.append('ambiguous_bases')
#     return ','.join(warnings) if warnings else '-'
