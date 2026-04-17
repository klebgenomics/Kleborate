"""
Copyright 2026 Mary Maranga
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


from ...shared.multi_mlst import multi_mlst
from ...shared.alignment import truncation_check
from ...shared.misc import load_fasta, reverse_complement
from .rmpADC_calling import poly_G_variation, poly_G_rmpC_variation, poly_A_variation, check_argR_box, check_argR_status,check_polyT_tract, translate_nucl_to_prot, process_status_dict, get_gene_status


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
    group.add_argument('--klebsiella__rmst_required_exact_matches', type=int, default=2,
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
    alleles_files = {gene: data_dir() / f'{gene}.fasta' for gene in genes}
    
    # Load status dictionaries
    rmpA_dict = process_status_dict(data_dir()/'rmpA_polyG_status.txt', 'rmpA')
    rmpC_dict = process_status_dict(data_dir()/'rmpC_polyG_status.txt', 'rmpC')
    rmpD_dict = process_status_dict(data_dir()/'rmpD_polyA_status.txt', 'rmpD')

    results, spurious_hits, hits_per_gene = multi_mlst(
        assembly, minimap2_index, profiles, alleles_files, genes,
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
    if st == 'NA': st = 0
    else: st = st[2:]

    spurious_hits_list = [item for h in spurious_hits.values() for item in h]
    spurious_virulence_hits = ';'.join(spurious_hits_list) if spurious_hits_list else '-'

    # rmpA Promoter Checks
    rmpA_allele = alleles.get('rmpA', None)
    has_rmpA = rmpA_allele and rmpA_allele != '-' and rmpA_allele.strip() != ''

    rmpA_promoter = "-"
    promoter_argR = "-" 
    promoter_polyT = "-"

    if has_rmpA:
        promoter_polyT = check_polyT_tract(hits_per_gene, assembly)
        promoter_argR = check_argR_box(hits_per_gene, assembly)
        
        if promoter_polyT != "-" and promoter_argR != "-":
            rmpA_promoter = f"{promoter_polyT} ({promoter_argR})"
        elif promoter_polyT != "-":
            rmpA_promoter = promoter_polyT
        elif promoter_argR != "-":
            rmpA_promoter = promoter_argR
        else:
            rmpA_promoter = "-"


    # checks for argR 
    argR_status = check_argR_status(assembly, argR_ref, args.klebsiella__rmst_min_identity, args.klebsiella__rmst_min_coverage) 
    
    # rmpADC status
    rmpA_status = get_gene_status(alleles.get('rmpA'), rmpA_dict, hits_per_gene, poly_G_variation)
    rmpD_status = get_gene_status(alleles.get('rmpD'), rmpD_dict, hits_per_gene, poly_A_variation)
    rmpC_status = get_gene_status(alleles.get('rmpC'), rmpC_dict, hits_per_gene, poly_G_rmpC_variation)
    
    current_statuses = {'rmpA': rmpA_status, 'rmpD': rmpD_status, 'rmpC': rmpC_status}
    for gene in genes:
        status_str = str(current_statuses[gene])
        if "OFF" in status_str:
            alleles[gene] = status_str
            
    combined_status_f = [str(rmpA_status), str(rmpD_status), str(rmpC_status)]
    
    # Determine overall Phase Status 
    existing_statuses = [s for s in combined_status_f if s != "-"]
    
    has_reversible_off = any("OFF" in s for s in existing_statuses)
    has_true_truncation = any(("%" in s and "OFF" not in s) for s in existing_statuses)

    # annotate lineage with partial if any gene is truncated or incomplete
    if hasattr(lineage, "replace"):
        is_partial = any(x in lineage for x in ["(truncated)", "(incomplete)", "(unknown)"])
        lineage = lineage.replace(" (truncated)", "").replace(" (incomplete)", "").replace(" (unknown)", "")
        
        if is_partial:
            lineage = f"{lineage} (partial)"

    # Determine RmpADC_status
    # if any gene is missing the rmpADC_status is set as '-'
    if any(s == "-" for s in combined_status_f):
        RmpADC_status = "-"
    elif has_true_truncation:
        RmpADC_status = "-"
    elif has_reversible_off:
        RmpADC_status = "Phase OFF"
    else:
        RmpADC_status = "Phase ON"
    
    
    if RmpADC_status != "-":
        annotation_groups = []
        
        # annotate the rmpADC status column with rmpA_promoter annotations
        promoter_anns = []
        if promoter_polyT != "-":
            if "reduced expression" in promoter_polyT: promoter_anns.append("reduced expression")
            elif "untypable" in promoter_polyT: promoter_anns.append("untypable")
        if promoter_argR != "-" and "ARG-box lost" in promoter_argR:
            promoter_anns.append("ARG box lost")
        
        if promoter_anns:
            annotation_groups.append(f"({', '.join(promoter_anns)})")

        # annotate the rmpADC status column with argR status annotations
        argR_str = str(argR_status)
        argR_ann = ""
        if "truncated" in argR_str:
            argR_f = argR_status if isinstance(argR_status, str) else "".join(argR_status)
            argR_ann = f"argR {argR_f}"
        elif "-" in argR_str and "truncated" not in argR_str:
            argR_ann = "argR missing"
        
        if argR_ann:
            annotation_groups.append(f"({argR_ann})")
            
        if annotation_groups:
            RmpADC_status = f"{RmpADC_status} {', '.join(annotation_groups)}"
        
    if isinstance(lineage, str):
        lineage = lineage.lstrip('- ').strip()
        if lineage == "(partial)":
            lineage = "partial"
        if not lineage:
            lineage = "-"
            
    return {
        'RmST': st, 
        'RmpADC': lineage,
        'rmpA': alleles.get('rmpA'), 
        'rmpD': alleles.get('rmpD'), 
        'rmpC': alleles.get('rmpC'),
        'RmpADC_status': RmpADC_status,
        'rmpA_promoter': rmpA_promoter,
        'argR': argR_status,
        'spurious_rmst_hits': spurious_virulence_hits
    }




# def get_results(assembly, minimap2_index, args, previous_results):
#     argR_ref = data_dir() / 'argR.fasta'
#     genes = ['rmpA', 'rmpC', 'rmpD']
#     profiles = data_dir() / 'profiles.tsv'
#     alleles_files = {gene: data_dir() / f'{gene}.fasta' for gene in genes}
    
#     # Load status dictionaries
#     rmpA_dict = process_status_dict(data_dir()/'rmpA_polyG_status.txt', 'rmpA')
#     rmpC_dict = process_status_dict(data_dir()/'rmpC_polyG_status.txt', 'rmpC')
#     rmpD_dict = process_status_dict(data_dir()/'rmpD_polyA_status.txt', 'rmpD')

#     results, spurious_hits, hits_per_gene = multi_mlst(
#         assembly, minimap2_index, profiles, alleles_files, genes,
#         'rmp_lineage', 
#         args.klebsiella__rmst_min_identity,
#         args.klebsiella__rmst_min_coverage, 
#         args.klebsiella__rmst_required_exact_matches,
#         check_for_truncation=True, 
#         report_incomplete=True,
#         min_spurious_identity=args.klebsiella__rmst_min_spurious_identity,
#         min_spurious_coverage=args.klebsiella__rmst_min_spurious_coverage,
#         unknown_group_name='rmp unknown',
#         min_gene_count=args.klebsiella__rmst_min_gene_count
#     )
    
#     st, lineage, alleles = results
#     if st == 'NA': st = 0
#     else: st = st[2:]

#     spurious_hits_list = [item for h in spurious_hits.values() for item in h]
#     spurious_virulence_hits = ';'.join(spurious_hits_list) if spurious_hits_list else '-'

#     # rmpA Promoter Checks
#     rmpA_allele = alleles.get('rmpA', None)
#     has_rmpA = rmpA_allele and rmpA_allele != '-' and rmpA_allele.strip() != ''

#     rmpA_promoter = "-"
#     promoter_argR = "-" 
#     promoter_polyT = "-"

#     if has_rmpA:
#         promoter_polyT = check_polyT_tract(hits_per_gene, assembly)
#         promoter_argR = check_argR_box(hits_per_gene, assembly)
        
#         if promoter_polyT != "-" and promoter_argR != "-":
#             rmpA_promoter = f"{promoter_polyT} ({promoter_argR})"
#         elif promoter_polyT != "-":
#             rmpA_promoter = promoter_polyT
#         elif promoter_argR != "-":
#             rmpA_promoter = promoter_argR
#         else:
#             rmpA_promoter = "-"


#     # checks for argR 
#     argR_status = check_argR_status(assembly, argR_ref, args.klebsiella__rmst_min_identity, args.klebsiella__rmst_min_coverage) 
#     # rmpADC status
#     rmpA_status = get_gene_status(alleles.get('rmpA'), rmpA_dict, hits_per_gene, poly_G_variation)
#     print(rmpA_status)
#     rmpD_status = get_gene_status(alleles.get('rmpD'), rmpD_dict, hits_per_gene, poly_A_variation)
#     print(rmpD_status)
#     rmpC_status = get_gene_status(alleles.get('rmpC'), rmpC_dict, hits_per_gene, poly_G_rmpC_variation)
#     print(rmpC_status)
    
#     current_statuses = {'rmpA': rmpA_status, 'rmpD': rmpD_status, 'rmpC': rmpC_status}
#     for gene in genes:
#         status_str = str(current_statuses[gene])
#         if "OFF" in status_str:
#             alleles[gene] = status_str
#     combined_status_f = [str(rmpA_status), str(rmpD_status), str(rmpC_status)]
#     print(combined_status_f)
    
#     # Determine overall Phase Status
#     has_reversible_off = any("OFF" in s for s in combined_status_f)
#     has_true_truncation = any(("%" in s and "OFF" not in s) for s in combined_status_f)

#     if isinstance(lineage, str):
#         lineage = lineage.replace(" (truncated)", "").replace(" (incomplete)", "").replace(" (unknown)", "")

#     if has_true_truncation:
#         RmpADC_status = "-"
#         lineage = f"{lineage} (partial)"
#     elif has_reversible_off:
#         RmpADC_status = "Phase OFF"
#     elif any(s == "-" for s in combined_status_f):
#         RmpADC_status = "-"
#     else:
#         RmpADC_status = "Phase ON"
    

#     if RmpADC_status != "-":
#         annotation_groups = []
        
#         # annotate the rmpADC status
#         promoter_anns = []
#         if promoter_polyT:
#             if "reduced expression" in promoter_polyT: promoter_anns.append("reduced expression")
#             elif "untypable" in promoter_polyT: promoter_anns.append("untypable")
#         if promoter_argR and "ARG-box lost" in promoter_argR:
#             promoter_anns.append("ARG box lost")
        
#         if promoter_anns:
#             annotation_groups.append(f"({', '.join(promoter_anns)})")

#         # argR checks
#         argR_str = str(argR_status)
#         argR_ann = ""
#         if "truncated" in argR_str:
#             argR_f = argR_status if isinstance(argR_status, str) else "".join(argR_status)
#             argR_ann = f"argR {argR_f}"
#         elif "-" in argR_str and "truncated" not in argR_str:
#             argR_ann = "argR missing"
        
#         if argR_ann:
#             annotation_groups.append(f"({argR_ann})")
            
#         if annotation_groups:
#             RmpADC_status = f"{RmpADC_status} {', '.join(annotation_groups)}"
        
#     if isinstance(lineage, str):
#         lineage = lineage.lstrip('- ').strip()
#         if lineage == "(partial)":
#             lineage = "partial"
#         if not lineage:
#             lineage = "-"
    
#     print(RmpADC_status)     
#     return {
#         'RmST': st, 
#         'RmpADC': lineage,
#         'rmpA': alleles.get('rmpA'), 
#         'rmpD': alleles.get('rmpD'), 
#         'rmpC': alleles.get('rmpC'),
#         'RmpADC_status': RmpADC_status,
#         'rmpA_promoter': rmpA_promoter,
#         'argR': argR_status,
#         'spurious_rmst_hits': spurious_virulence_hits
#     }

