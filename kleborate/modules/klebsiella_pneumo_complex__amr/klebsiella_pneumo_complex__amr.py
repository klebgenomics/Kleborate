"""
Copyright 2025 Mary Maranga (gathonimaranga@gmail.com)
https://github.com/klebgenomics/KleborateModular/

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
import pandas as pd
from pathlib import Path
import csv
import re

from ...shared.resMinimap import read_class_file, get_res_headers, resminimap_assembly


def description():
    return 'Genotyping acquired genes and mutations for the Klebsiella pneumoniae species complex'


def prerequisite_modules():
    return []



def get_headers():
    _, res_classes, bla_classes = read_class_file(data_dir() / 'CARD_AMR_clustered.csv')
    res_headers = get_res_headers(res_classes, bla_classes)
    res_headers += ['truncated_resistance_hits', 'spurious_resistance_hits']
    

    full_headers = res_headers + [
                    'Input_sequence_ID','Input_gene_length', 'Input_gene_start', 'Input_gene_stop', 'Reference_gene_length', 
                    'Reference_gene_start', 'Reference_gene_stop', 'Sequence_identity', 'Coverage',
                    'Software_name', 'Software_version', 'Reference_database_name', 
                    'Reference_database_version', 'Reference_accession', 'Genetic_variation_type',
                    'Antimicrobial_agent', 'Coverage_depth', 'Coverage_ratio', 'Drug_class',
                    'Input_protein_length', 'Input_protein_start', 'Input_protein_stop',
                    'Predicted_phenotype','predicted_phenotype_confidence_level',
                    'Reference_protein_length', 'Reference_protein_start', 'Reference_protein_stop',
                    'Resistance_mechanism', 'Strand_orientation']

    stdout_headers = []
    return full_headers, stdout_headers


def add_cli_options(parser):
    module_name = os.path.basename(__file__)[:-3]
    group = parser.add_argument_group(f'{module_name} module')

    group.add_argument('--klebsiella_pneumo_complex__amr_min_identity', type=float, default=90.0,
                       help='Minimum alignment percent identity for klebsiella_pneumo_complex Amr results')
    group.add_argument('--klebsiella_pneumo_complex__amr_min_coverage', type=float, default=80.0,
                       help='Minimum alignment percent coverage for klebsiella_pneumo_complex Amr  results')
    group.add_argument('--klebsiella_pneumo_complex__amr_min_spurious_identity', type=float, default=80.0,
                       help='Minimum alignment percent identity for klebsiella_pneumo_complex Amr spurious results')
    group.add_argument('--klebsiella_pneumo_complex__amr_min_spurious_coverage', type=float, default=40.0,
                       help='Minimum alignment percent coverage for klebsiella_pneumo_complex Amr spurious results')
    
    return group


def check_cli_options(args):
    if args.klebsiella_pneumo_complex__amr_min_identity <= 50.0 or args.klebsiella_pneumo_complex__amr_min_identity >= 100.0:
        sys.exit('Error: --klebsiella_pneumo_complex__amr_min_identity must be between 50.0 and 100.0')
    if args.klebsiella_pneumo_complex__amr_min_coverage <= 50.0 or args.klebsiella_pneumo_complex__amr_min_coverage >= 100.0:
        sys.exit('Error: --klebsiella_pneumo_complex__amr_min_coverage must be between 50.0 and 100.0')
    if args.klebsiella_pneumo_complex__amr_min_spurious_identity <= 50.0 or args.klebsiella_pneumo_complex__amr_min_spurious_identity >= 100.0:
        sys.exit('Error: --klebsiella_pneumo_complex__amr_min_spurious_identity must be between 50.0 and 100.0')
    if args.klebsiella_pneumo_complex__amr_min_spurious_coverage <= 30.0 or args.klebsiella_pneumo_complex__amr_min_spurious_coverage >= 100.0:
        sys.exit('Error: --klebsiella_pneumo_complex__amr_min_spurious_coverage must be between 30.0 and 100.0')


def check_external_programs():
    if not shutil.which('minimap2'):
        sys.exit('Error: could not find minimap2')
    return ['minimap2']


def data_dir():
    return pathlib.Path(__file__).parents[0] / 'data'



def format_res_hits(res_hits, full_headers):
    specific_keys = [
        'Input_sequence_ID','Input_gene_length', 'Input_gene_start', 'Input_gene_stop', 'Reference_gene_length',
        'Reference_gene_start', 'Reference_gene_stop', 'Sequence_identity', 'Coverage',
        'Reference_accession', 'Genetic_variation_type', 'Antimicrobial_agent', 'Coverage_depth',
        'Coverage_ratio', 'Drug_class', 'Input_protein_length', 'Input_protein_start',
        'Input_protein_stop','Predicted_phenotype', 'predicted_phenotype_confidence_level', 
        'Reference_protein_length', 'Reference_protein_start',
        'Reference_protein_stop', 'Resistance_mechanism', 'strand_orientation'
    ]
    software_fields = [
        'Software_name', 'Software_version', 'Reference_database_name', 'Reference_database_version'
    ]
    extra_categories = [
        'truncated_resistance_hits', 'spurious_resistance_hits', 'Col_mutations', 'Omp_mutations', 'Flq_mutations'
    ]

    formatted_dict = {}

    # 1. Process all categories in res_hits
    for category in set(res_hits) | set(extra_categories):
        if category in software_fields:
            continue
        values = []
        data = res_hits.get(category)
        if data is not None:
            if isinstance(data, str):
                values.append(data)
            elif isinstance(data, list):
                for hit in data:
                    if isinstance(hit, str):
                        values.append(hit)
                    elif isinstance(hit, list) and len(hit) > 0:
                        values.append(str(hit[0]))
                    else:
                        values.append(str(hit))
            else:
                values.append(str(data))
        if values:
            formatted_dict[category] = [";".join(map(str, values))]

    # 2. Extract values for specific fields from nested dicts in lists
    for field in full_headers:
        if field in software_fields or field in formatted_dict:
            continue
        values = []
        for data in res_hits.values():
            if isinstance(data, list):
                for hit in data:
                    if isinstance(hit, list) and len(hit) > 1:
                        for idx in (1, 2):  # Check second and third element if they exist
                            if len(hit) > idx and isinstance(hit[idx], dict) and field in hit[idx]:
                                values.append(f"{hit[0]}:{hit[idx][field]}")
        if values:
            formatted_dict[field] = [";".join(values)]

    # 3. Add software fields
    for field in software_fields:
        if field in res_hits:
            formatted_dict[field] = [res_hits[field]]

    return formatted_dict



def get_results(assembly, minimap2_index, args, previous_results):
    # Read gene info and headers
    gene_info, _, _ = read_class_file(data_dir() / 'CARD_AMR_clustered.csv')
    full_headers, _ = get_headers()
    qrdr = data_dir() / 'QRDR_120.fasta'
    trunc = data_dir() / 'MgrB_and_PmrB.fasta'
    omp = data_dir() / 'OmpK.fasta'
    ref_file = data_dir() / 'CARD_v3.2.9.fasta'

    # Run minimap and get results
    res_hits = resminimap_assembly(
        assembly,
        minimap2_index,
        ref_file,
        gene_info,
        qrdr,
        trunc,
        omp,
        args.klebsiella_pneumo_complex__amr_min_coverage,
        args.klebsiella_pneumo_complex__amr_min_identity,
        args.klebsiella_pneumo_complex__amr_min_spurious_coverage,
        args.klebsiella_pneumo_complex__amr_min_spurious_identity
    )

    # --- report in aac(6`)-Ib-cr* allele in both AGly_acquired and Flq_acquired---
    aac6_pattern = re.compile(r"aac\(6[’']\)-Ib-cr", re.IGNORECASE)
    agly_hits = res_hits.get('AGly_acquired', [])
    flq_hits = res_hits['Flq_acquired']

    for hit in agly_hits:
        if not isinstance(hit, list) or len(hit) < 2:
            continue
        gene_string = hit[0]
        # Split gene_string
        for gene in re.split(r'[;^]', gene_string):
            if aac6_pattern.search(gene.strip()):
                flq_hits.append(hit)
                break

    # --- software and database metadata ---
    res_hits['Software_name'] = 'Kleborate'
    res_hits['Software_version'] = '3.1.3'
    res_hits['Reference_database_name'] = 'CARD'
    res_hits['Reference_database_version'] = '3.2.9'

    # --- map alleles to CARD ARO_accessions and CARD_class---
    allele_to_accession = {}
    allele_to_drug_class = {}
    drug_class_to_accession = {}

    with open(data_dir() / 'CARD_AMR_clustered.csv', newline='', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',', skipinitialspace=True)
        for row in reader:
            allele = row['allele'].strip()
            accession = row['ARO_class'].strip()
            drug_class = row.get('CARD_class', '-').strip()
            if allele:
                allele_to_accession[allele] = accession
                allele_to_drug_class[allele] = drug_class
            # Build reverse mapping for mutation lookup
            if drug_class and accession:
                if drug_class not in drug_class_to_accession:
                    drug_class_to_accession[drug_class] = accession

    sorted_alleles = sorted(allele_to_accession.keys(), key=len, reverse=True)

    # --- Define mutation drug class mapping ---
    mutation_drug_class = {
        'SHV_mutations': 'penicillin beta-lactam',
        'Col_mutations': 'peptide antibiotic',
        'Omp_mutations': '',  # OMP is not class specific
        'Flq_mutations': 'fluoroquinolone antibiotic'
    }

    # --- Reference accession and drug class assignment ---
    for category in res_hits:
        if not isinstance(res_hits[category], list):
            continue

        #  retrieve drug class for mutation category
        mutation_class = mutation_drug_class.get(category, None)
        mutation_accession = None

        if mutation_class:
            mutation_accession = drug_class_to_accession.get(mutation_class, '-')

        # Iterate through individual hits
        for hit in res_hits[category]:
            # process only lists with metadata
            if not isinstance(hit, list) or len(hit) < 2:
                continue
            # extract gene name
            gene_name = hit[0]
            merged = {}
            # merge all metadata
            for elem in hit[1:]:
                if isinstance(elem, dict):
                    merged.update(elem)
            
            matched_accession = None
            matched_drug_class = None

            # Mutation fields
            if mutation_class is not None:
                # use the mapped class, and accession 
                matched_drug_class = mutation_class
                matched_accession = mutation_accession
            else:
                # map gene to accession
                if gene_name in allele_to_accession:
                    matched_accession = allele_to_accession[gene_name]
                    matched_drug_class = allele_to_drug_class.get(gene_name, '-')
                else:
                    # find the longest matching prefix
                    for allele in sorted_alleles:
                        if gene_name.startswith(allele):
                            matched_accession = allele_to_accession[allele]
                            matched_drug_class = allele_to_drug_class.get(allele, '-')
                            break
            # Add accession and drug class
            merged['Reference_accession'] = matched_accession if matched_accession else '-'
            merged['Drug_class'] = matched_drug_class if matched_drug_class else '-'

            for elem in hit[1:]:
                if isinstance(elem, dict):
                    elem['Reference_accession'] = merged['Reference_accession']
                    elem['Drug_class'] = merged['Drug_class']

    # Format the res_hits dictionary
    res_hits = format_res_hits(res_hits, full_headers)

    # Double check that all results correspond to a header in full_headers
    for h in res_hits.keys():
        if h not in full_headers:
            sys.exit(f'Error: results contained a value ({h}) that is not covered by the full_headers')

    # Return a dictionary with values from res_hits or '-'
    return {r: res_hits[r] if r in res_hits else '-' for r in full_headers}



# def get_results(assembly, minimap2_index, args, previous_results):
#     gene_info, _, _ = read_class_file(data_dir() / 'CARD_AMR_clustered.csv')
#     full_headers, _ = get_headers() 
#     qrdr = data_dir() / 'QRDR_120.fasta'
#     trunc = data_dir() / 'MgrB_and_PmrB.fasta'
#     omp = data_dir() / 'OmpK.fasta'

#     ref_file = data_dir() / 'CARD_v3.2.9.fasta'

#     res_hits = resminimap_assembly(
#         assembly,
#         minimap2_index, 
#         ref_file, 
#         gene_info, 
#         qrdr, 
#         trunc, 
#         omp,   
#         args.klebsiella_pneumo_complex__amr_min_coverage, 
#         args.klebsiella_pneumo_complex__amr_min_identity,
#         args.klebsiella_pneumo_complex__amr_min_spurious_coverage,
#         args.klebsiella_pneumo_complex__amr_min_spurious_identity
#     )

#     # Double check that there weren't any results without a corresponding full headers.
#     for h in res_hits.keys():
#         if h not in full_headers:
#             sys.exit( f'Error: results contained a value ({h}) that is not covered by the '
#                       f'full headers')

#     return {r: ';'.join(sorted(res_hits[r])) if r in res_hits else '-' for r in full_headers}



