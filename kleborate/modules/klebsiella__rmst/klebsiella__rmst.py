"""
Copyright 2025 Kat Holt, Mary Maranga, Ryan Wick
https://github.com/katholt/Kleborate/

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
import re


from ...shared.multi_mlst import multi_mlst, poly_G_variation, check_argR_box, check_argR_status,check_polyT_tract, translate_nucl_to_prot, allele_type
from ...shared.alignment import truncation_check
from ...shared.misc import load_fasta, reverse_complement

def description():
    return 'MLST on the KpSC Rmp locus (rmp genes)'


def prerequisite_modules():
    return []


def get_headers():
    full_headers = ['RmST', 'RmpADC', 'RmpADC_status','rmpA', 'rmpD', 'rmpC', 'rmpA_promoter', 'argR','spurious_rmst_hits']
    stdout_headers = []
    return full_headers, stdout_headers


def add_cli_options(parser):
    module_name = os.path.basename(__file__)[:-3]
    group = parser.add_argument_group(f'{module_name} module')
    group.add_argument('--klebsiella__rmst_min_identity', type=float, default=90.0,
                       help='Minimum alignment percent identity for Rmp MLST')
    group.add_argument('--klebsiella__rmst_min_coverage', type=float, default=80.0,
                       help='Minimum alignment percent coverage for Rmp MLST')
    group.add_argument('--klebsiella__rmst_min_spurious_identity', type=float, default=80.0,
                       help='Minimum alignment percent identity for klebsiella__rmst spurious results')
    group.add_argument('--klebsiella__rmst_min_spurious_coverage', type=float, default=40.0,
                       help='Minimum alignment percent coverage for klebsiella__rmst spurious results')
    group.add_argument('--klebsiella__rmst_required_exact_matches', type=int, default=1,
                       help='At least this many exact matches are required to call an ST')
    group.add_argument('--klebsiella__rmst_min_gene_count', type=int, default=2,
                       help='At least this many exact alleles required to report a novel allele')
    return group


def check_cli_options(args):
    if args.klebsiella__rmst_min_identity <= 50.0 or args.klebsiella__rmst_min_identity >= 100.0:
        sys.exit('Error: --klebsiella__rmst_min_identity must be between 50.0 and 100.0')
    if args.klebsiella__rmst_min_coverage <= 50.0 or args.klebsiella__rmst_min_coverage >= 100.0:
        sys.exit('Error: --klebsiella__rmst_min_coverage must be between 50.0 and 100.0')
    if args.klebsiella__rmst_min_spurious_identity <= 50.0 or args.klebsiella__rmst_min_spurious_identity >= 100.0:
        sys.exit('Error: --klebsiella__rmst_min_spurious_identity must be between 50.0 and 100.0')
    if args.klebsiella__rmst_min_spurious_coverage <= 30.0 or args.klebsiella__rmst_min_spurious_coverage >= 100.0:
        sys.exit('Error: --klebsiella__rmst_min_spurious_coverage must be between 30.0 and 100.0')
    if args.klebsiella__rmst_min_gene_count < 0:
        sys.exit('Error: --klebsiella__rmst_min_gene_count must be a positive integer')


def check_external_programs():
    if not shutil.which('minimap2'):
        sys.exit('Error: could not find minimap2')
    return ['minimap2']


def data_dir():
    return pathlib.Path(__file__).parents[0] / 'data'


def get_results(assembly, minimap2_index, args, previous_results):
    argR_ref = data_dir() / 'argR.fasta'
    genes = ['rmpA', 'rmpC', 'rmpD']
    profiles = data_dir() / 'profiles.tsv'
    alleles = {gene: data_dir() / f'{gene}.fasta' for gene in genes}
    rmpA_status_dict = data_dir()/'rampA_polyG_status.txt'


    with open(rmpA_status_dict) as f:
        rmpA_dict_raw = ast.literal_eval(f.read())
    rmpA_dict = {
        re.sub(r'^rmpA_', '', k).rstrip('*'): v
        for k, v in rmpA_dict_raw.items()
    }

    results, spurious_hits, hits_per_gene  = multi_mlst(assembly, minimap2_index, profiles, alleles, genes,
                                      'rmp_lineage', args.klebsiella__rmst_min_identity,
                                      args.klebsiella__rmst_min_coverage, args.klebsiella__rmst_required_exact_matches,
                                      check_for_truncation=True, report_incomplete=True,
                                      min_spurious_identity=args.klebsiella__rmst_min_spurious_identity,
                                      min_spurious_coverage=args.klebsiella__rmst_min_spurious_coverage,
                                      unknown_group_name='rmp unknown',
                                      min_gene_count=args.klebsiella__rmst_min_gene_count)
    st, lineage, alleles = results

    if st == 'NA':
        st = 0
    else:
        st = st[2:]

    # spurious hits
    spurious_hits = [item for h in spurious_hits.values() for item in h]

    spurious_virulence_hits = ';'.join(spurious_hits )if spurious_hits else '-'

    rmpA_allele = alleles.get('rmpA', None)
    has_rmpA = rmpA_allele and rmpA_allele != '-' and rmpA_allele.strip() != ''

    if has_rmpA:
        promoter_polyT = check_polyT_tract(hits_per_gene, assembly)
        promoter_argR = check_argR_box(hits_per_gene, assembly)
        if promoter_argR:
            rmpA_promoter = f"{promoter_polyT}, {promoter_argR}"
        else:
            rmpA_promoter = f"{promoter_polyT}"
    else:
        
        promoter_polyT = "-"
        promoter_argR = "-"
        rmpA_promoter = "-"

    # rmpA_status logic
    allele_value = alleles['rmpA']
    allele_key = re.sub(r'^rmpA_', '', allele_value).rstrip('*')

    if allele_key == '-' or allele_value == '-':
        RmpADC_status = "-"
    else:
        match_allele_type = allele_type(allele_value)
        if match_allele_type == "truncated":
            RmpADC_status = "truncated"
        elif match_allele_type == "inexact":
            RmpADC_status = poly_G_variation(hits_per_gene)
        elif match_allele_type == "exact":
            RmpADC_status = rmpA_dict[allele_value][0]
        else:
            RmpADC_status = "-"

    # add rmpA_promoter
    def append_status_annotation(status, annotation):
        status = status.strip()
        annotation = annotation.strip()
        if not annotation:
            return status
        # Find existing parenthetical annotation
        if status.endswith(")"):
            # Insert new annotation inside the existing parentheses
            return status[:-1] + ", " + annotation + ")"
        else:
            # Add a new parenthetical annotation
            return f"{status} ({annotation})"

    # Add (promoter 10T) if polyT is 10T (OFF)
    if str(promoter_polyT).strip() == "10T (OFF)":
        RmpADC_status = append_status_annotation(RmpADC_status, "promoter 10T")

    # Add (ARG box lost) if reported in promoter_argR
    if promoter_argR and "ARG box lost" in promoter_argR:
        RmpADC_status = append_status_annotation(RmpADC_status, "ARG box lost")

    argR_status = check_argR_status(hits_per_gene, assembly, argR_ref, args.klebsiella__rmst_min_identity,args.klebsiella__rmst_min_coverage)                           
                          
    return {'RmST': st, 
            'RmpADC': lineage,
            'rmpA': alleles['rmpA'], 'rmpD': alleles['rmpD'], 'rmpC': alleles['rmpC'],
            'RmpADC_status': RmpADC_status,
            'rmpA_promoter': rmpA_promoter,
            'argR': argR_status,
            'spurious_rmst_hits':spurious_virulence_hits}



# def get_results(assembly, minimap2_index, args, previous_results):
#     genes = ['rmpA', 'rmpC', 'rmpD']
#     profiles = data_dir() / 'profiles.tsv'
#     alleles = {gene: data_dir() / f'{gene}.fasta' for gene in genes}

#     results, spurious_hits  = multi_mlst(assembly, minimap2_index, profiles, alleles, genes,
#                                       'rmp_lineage', args.klebsiella__rmst_min_identity,
#                                       args.klebsiella__rmst_min_coverage, args.klebsiella__rmst_required_exact_matches,
#                                       check_for_truncation=True, report_incomplete=True,
#                                       min_spurious_identity=args.klebsiella__rmst_min_spurious_identity,
#                                       min_spurious_coverage=args.klebsiella__rmst_min_spurious_coverage,
#                                       unknown_group_name='rmp unknown',
#                                       min_gene_count=args.klebsiella__rmst_min_gene_count)
#     st, lineage, alleles = results

#     if st == 'NA':
#         st = 0
#     else:
#         st = st[2:]

#     # spurious hits
#     spurious_hits = [item for h in spurious_hits.values() for item in h]

#     spurious_virulence_hits = ';'.join(spurious_hits )if spurious_hits else '-'

#     return {'RmST': st, 'RmpADC': lineage,
#             'rmpA': alleles['rmpA'], 'rmpD': alleles['rmpD'], 'rmpC': alleles['rmpC'],
#             'spurious_rmst_hits':spurious_virulence_hits}

