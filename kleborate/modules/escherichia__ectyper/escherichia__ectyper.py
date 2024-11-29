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
import json
import pandas as pd
import logging

from .ectyper import load_json_file, alignment, minimap2_output_to_df, ectyper_dict_to_df, setAlleleMeta, get_prediction, predict_serotype, process_sample_data

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)


def description():
    return 'E. coli serotyping of O and H antigens'


def prerequisite_modules():
    return []


def get_headers():
    full_headers = ['Otype', 'Htype', 'Serotype', 'Evidence', 'AlleleKeys',
    'GeneIdentities'] 
    stdout_headers = []
    return full_headers, stdout_headers


def add_cli_options(parser):
    module_name = os.path.basename(__file__)[:-3]
    group = parser.add_argument_group(f'{module_name} module')
    group.add_argument('--escherichia__serotyping_min_identity', type=float, default=90.0,
                       help='Minimum alignment percent identity for detecting serotypes')
    group.add_argument('--escherichia__serotyping_min_coverage', type=float, default=80.0,
                       help='Minimum alignment percent coverage for detecting serotypes')
    # New arguments
    group.add_argument('--escherichia__serotyping_percentIdentityOtype', type=float, default=90.0,
                       help='Minimum percent identity for O-type detection')
    group.add_argument('--escherichia__serotyping_percentIdentityHtype', type=float, default=95.0,
                       help='Minimum percent identity for H-type detection')
    group.add_argument('--escherichia__serotyping_percentCoverageOtype', type=float, default=90.0,
                       help='Minimum percent coverage for O-type detection')
    group.add_argument('--escherichia__serotyping_percentCoverageHtype', type=float, default=50.0,
                       help='Minimum percent coverage for H-type detection')
    group.add_argument('--escherichia__serotyping_HIGH_SIMILARITY_THRESHOLD_O', type=float, default=0.00771,
                       help='High similarity threshold for O-type')
    group.add_argument('--escherichia__serotyping_MIN_O_IDENTITY_LS', type=float, default=95.0,
                       help='Minimum percent identity for O-type (LS threshold)')
    group.add_argument('--escherichia__serotyping_MIN_O_COVERAGE_LS', type=float, default=48.0,
                       help='Minimum percent coverage for O-type (LS threshold)')
    return


def check_cli_options(args):
    if not (50.0 <= args.escherichia__serotyping_min_identity <= 100.0):
        sys.exit('Error: --escherichia__serotyping_min_identity must be between 50.0 and 100.0')
    if not (50.0 <= args.escherichia__serotyping_min_coverage <= 100.0):
        sys.exit('Error: --escherichia__serotyping_min_coverage must be between 50.0 and 100.0')

    # Additional checks for new arguments
    if not (50.0 <= args.escherichia__serotyping_percentIdentityOtype <= 100.0):
        sys.exit('Error: --percentIdentityOtype must be between 50.0 and 100.0')
    if not (50.0 <= args.escherichia__serotyping_percentIdentityHtype <= 100.0):
        sys.exit('Error: --percentIdentityHtype must be between 50.0 and 100.0')
    if not (50.0 <= args.escherichia__serotyping_percentCoverageOtype <= 100.0):
        sys.exit('Error: --percentCoverageOtype must be between 50.0 and 100.0')
    if not (50.0 <= args.escherichia__serotyping_percentCoverageHtype <= 100.0):
        sys.exit('Error: --percentCoverageHtype must be between 50.0 and 100.0')
    if not (50.0 <= args.escherichia__serotyping_MIN_O_IDENTITY_LS <= 100.0):
        sys.exit('Error: --MIN_O_IDENTITY_LS must be between 50.0 and 100.0')
    if not (10.0 <= args.escherichia__serotyping_MIN_O_COVERAGE_LS <= 100.0):
        sys.exit('Error: --MIN_O_COVERAGE_LS must be between 10.0 and 100.0')

    return

    

def check_external_programs():
    if not shutil.which('minimap2'):
        sys.exit('Error: could not find minimap2')
    return ['minimap2']



def data_dir():
    return pathlib.Path(__file__).parents[0] / 'data'



def get_results(assembly, minimap2_index, args, previous_results):
    # Step 1: Load the JSON database
    ref_db = data_dir() / 'ectyper_alleles_db.json'
    ref_file = data_dir() / 'ectyper_alleles_db.fasta'
    ectyper_dict = load_json_file(ref_db)
    if ectyper_dict is None:
        sys.exit(f"Error: Could not load the JSON reference database from {ref_db}")

    # Step 2: Perform the alignment
    alignment_hits = alignment(
        ref_file,
        minimap2_index,
        assembly,
        args.escherichia__serotyping_min_identity,
        args.escherichia__serotyping_min_coverage
    )

    # Step 3: Predict the serotype
    predictions_dict, output_df = predict_serotype(alignment_hits, ectyper_dict, args, previous_results)

    # Step 4: Process each data for final serotype prediction
    final_results_dict = {}
    for genome_name in predictions_dict:
        final_results_dict[genome_name] = process_sample_data(genome_name, predictions_dict, ectyper_dict)

    # Step 5: Generate the output dictionary including the full headers
    combined_results = {}  # Initialize combined_results dictionary
    output_results = []
    for genome_name, result in final_results_dict.items():
        # Extract the necessary values based on headers
        otype = result.get('Otype', '-')
        htype = result.get('Htype', '-')
        serotype = result.get('Serotype', '-')
        # gene_scores = result.get('gene_scores', '-')
        evidence = result.get('O_allele_No', ['-'])[0]
        allele_keys = result.get('alleles', '-')
        # qc = result.get('QC', '-')

        # Gene identities, contig names, coordinates, and lengths
        gene_identities = []
        gene_contig_names = []
        gene_ranges = []
        gene_lengths = []

        # Iterate over all genes
        for gene, gene_data in result.items():
            if isinstance(gene_data, dict):
                gene_identities.append(str(gene_data.get('identity', '-')))
                gene_contig_names.append(gene_data.get('contigname', '-'))
                gene_ranges.append(gene_data.get('coordinates', '-'))
                gene_lengths.append(str(gene_data.get('length', '-')))

        # If no gene data is found, append '-'
        if not gene_identities:
            gene_identities.append('-')
        if not gene_contig_names:
            gene_contig_names.append('-')
        if not gene_ranges:
            gene_ranges.append('-')
        if not gene_lengths:
            gene_lengths.append('-')

        # Join the lists with semicolons
        gene_identities_str = ';'.join(gene_identities)
        gene_contig_names_str = ';'.join(gene_contig_names)
        gene_ranges_str = ';'.join(gene_ranges)
        gene_lengths_str = ';'.join(gene_lengths)

        results = {
            'Otype': otype,
            'Htype': htype,
            'Serotype': serotype,
            # 'gene_scores': gene_scores,
            'Evidence': evidence,
            'AlleleKeys': allele_keys,
            'GeneIdentities': gene_identities_str,
            # 'GeneContigNames': gene_contig_names_str,
            # 'Gene ranges': gene_ranges_str,
            # 'GeneLengths': gene_lengths_str
            # 'QC':qc
        }

        full_headers, _ = get_headers()
        for h in results.keys():
            if h not in full_headers:
                sys.exit(f'Error: results contained a value ({h}) that is not covered by the full headers')

        combined_results.update({
            r: ';'.join(sorted(results[r].split(';'))) if isinstance(results[r], str) and r in [
            'GeneIdentities' 
            ] else results[r] for r in full_headers
        })


        # combined_results.update({
        #     r: ';'.join(sorted(results[r].split(';'))) if isinstance(results[r], str) and r in [
        #         'GeneIdentities', 'GeneContigNames', 'Gene ranges', 'GeneLengths'
        #     ] else results[r] for r in full_headers
        # })

    return combined_results




# def get_results(assembly, minimap2_index, args, previous_results):
#     # Step 1: Load the JSON database
#     ref_db = data_dir() / 'ectyper_alleles_db.json'
#     ref_file = data_dir() / 'ectyper_alleles_db.fasta'
#     ectyper_dict = load_json_file(ref_db)
#     if ectyper_dict is None:
#         sys.exit(f"Error: Could not load the JSON reference database from {ref_db}")

#     # Step 2: Perform the alignment
#     alignment_hits = alignment(ref_file, 
#         minimap2_index, 
#         assembly, 
#         args.escherichia__serotyping_min_identity,
#         args.escherichia__serotyping_min_coverage)
#     if not alignment_hits:
#         sys.exit("Error: Alignment failed or no hits found.")

#     # Step 3: Predict the serotype of all genomes
#     predictions_dict, output_df = predict_serotype(alignment_hits, ectyper_dict, args)

#     # Step 4: Process each sample's data for final serotype prediction
#     final_results_dict = {}
#     for genome_name in predictions_dict:
#         final_results_dict[genome_name] = process_sample_data(genome_name, predictions_dict, ectyper_dict)

#     return final_results_dict


