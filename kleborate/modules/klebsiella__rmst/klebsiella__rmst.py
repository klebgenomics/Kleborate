"""
Copyright 2025  Mary Maranga
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


from ...shared.multi_mlst import multi_mlst, poly_G_variation, poly_G_rmpC_variation, poly_A_variation, check_argR_box, check_argR_status,check_polyT_tract, translate_nucl_to_prot, allele_type, process_status_dict, get_gene_status
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
    
    rmpA_status_dict = data_dir()/'rmpA_polyG_status.txt'
    rmpC_status_dict = data_dir()/'rmpC_polyG_status.txt' 
    rmpD_status_dict = data_dir()/'rmpD_polyA_status.txt'

    rmpA_dict = process_status_dict(rmpA_status_dict, 'rmpA')
    rmpC_dict = process_status_dict(rmpC_status_dict, 'rmpC')
    rmpD_dict = process_status_dict(rmpD_status_dict, 'rmpD')

    results, spurious_hits, hits_per_gene = multi_mlst(
        assembly, minimap2_index, profiles, alleles, genes,
        'rmp_lineage', 
        args.klebsiella__rmst_min_identity,
        args.klebsiella__rmst_min_coverage, 
        args.klebsiella__rmst_required_exact_matches,
        check_for_truncation=True, 
        report_incomplete=True,
        min_spurious_identity=args.klebsiella__rmst_min_spurious_identity,
        min_spurious_coverage=args.klebsiella__rmst_min_spurious_coverage,
        unknown_group_name='rmp unknown',
        min_gene_count=args.klebsiella__rmst_min_gene_count
    )
    
    st, lineage, alleles = results
    # print(alleles)

    if st == 'NA':
        st = 0
    else:
        st = st[2:]

    spurious_hits = [item for h in spurious_hits.values() for item in h]
    spurious_virulence_hits = ';'.join(spurious_hits) if spurious_hits else '-'

    # Promoter Checks (PolyT and ArgR box)
    # checks if the rmpA gene is found
    rmpA_allele = alleles.get('rmpA', None)
    has_rmpA = rmpA_allele and rmpA_allele != '-' and rmpA_allele.strip() != ''

    rmpA_promoter = "-"
    promoter_argR = "-" 

    if has_rmpA:
        promoter_polyT = check_polyT_tract(hits_per_gene, assembly)
        promoter_argR = check_argR_box(hits_per_gene, assembly)
        
        promoter_string = [x for x in [promoter_polyT, promoter_argR] if x != "-"]
        if promoter_string:
            rmpA_promoter = ", ".join(promoter_string)
        else:
            rmpA_promoter = "-"

    # Check argR gene status
    argR_status = check_argR_status( 
        assembly, 
        argR_ref, 
        args.klebsiella__rmst_min_identity,
        args.klebsiella__rmst_min_coverage
    ) 

    rmpA_status = get_gene_status(alleles['rmpA'], rmpA_dict, hits_per_gene, poly_G_variation)
    rmpD_status = get_gene_status(alleles['rmpD'], rmpD_dict, hits_per_gene, poly_A_variation)
    rmpC_status = get_gene_status(alleles['rmpC'], rmpC_dict, hits_per_gene, poly_G_rmpC_variation)
    
    if "OFF" in str(rmpA_status):
        alleles['rmpA'] = str(rmpA_status)
    
    if "OFF" in str(rmpD_status):
        alleles['rmpD'] = str(rmpD_status)

    if "OFF" in str(rmpC_status):
        alleles['rmpC'] = str(rmpC_status)

    combined_status = [rmpA_status, rmpD_status, rmpC_status]
    
    has_reversible_off = any("OFF" in str(s) for s in combined_status)
    has_true_truncation = any("%" in str(s) for s in combined_status)

    if has_reversible_off and not has_true_truncation:
        if isinstance(lineage, str):
            lineage = lineage.replace(" (truncated)", "").replace(" (incomplete)", "").replace(" (unknown)", "")

    if isinstance(lineage, str) and ("truncated" in lineage or "incomplete" in lineage or "unknown" in lineage):
        lineage = lineage.replace(" (truncated)", "").replace(" (incomplete)", "").replace(" (unknown)", "")
        lineage = f"{lineage} (partial)"
        RmpADC_status = "-"
    else:
        if any("OFF" in str(status) for status in combined_status):
            RmpADC_status = "OFF"
        elif any(str(status) == "-" for status in combined_status):
            RmpADC_status = "-"
        else:
            RmpADC_status = "ON"
    
    # RmpADC_status annotation
    if RmpADC_status != "-":
        ann = []
        # check rmpA promoter
        if promoter_polyT:
            if "reduced expression" in promoter_polyT:
                ann.append("reduced expression")
            elif "untypable" in promoter_polyT:
                ann.append("untypable")

        # Check ARG-box lost
        if promoter_argR and "ARG-box lost" in promoter_argR:
            ann.append("ARG box lost")
        
        # Check argR gene status
        argR_str = str(argR_status) 
        
        if "truncated" in argR_str:
            clean_status = argR_status if isinstance(argR_status, str) else "".join(argR_status)
            ann.append(f"argR {clean_status}")
            
        elif "-" in argR_str and "truncated" not in argR_str:
            ann.append("argR missing")
            
        if ann:
            RmpADC_status = f"{RmpADC_status} ({', '.join(ann)})"
            
    return {'RmST': st, 
            'RmpADC': lineage,
            'rmpA': alleles['rmpA'], 
            'rmpD': alleles['rmpD'], 
            'rmpC': alleles['rmpC'],
            'RmpADC_status': RmpADC_status,
            'rmpA_promoter': rmpA_promoter,
            'argR': argR_status,
            'spurious_rmst_hits': spurious_virulence_hits}



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

