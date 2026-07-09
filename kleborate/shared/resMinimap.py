"""
Copyright 2025 Mary Maranga
https://github.com/klebgenomics/Kleborate

This file is part of Kleborate. Kleborate is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Kleborate is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Kleborate. If
not, see <http://www.gnu.org/licenses/>.
"""

import collections
from collections import defaultdict
from Bio.Seq import Seq
from Bio import Align
from Bio.Align import substitution_matrices
from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError
import traceback
import re
 
from .alignment import align_query_to_ref, cull_redundant_hits, is_exact_aa_match, translate_nucl_to_prot, check_for_exact_aa_match, truncation_check,get_bases_per_ref_pos, translate_nucl_to_prot
from .misc import load_fasta, reverse_complement
from kleborate.modules.kpsc__amr.shv_mutations import*
from kleborate.modules.kpsc__amr.qrdr_mutations import*
from kleborate.modules.kpsc__amr.omp_mutations import*
from kleborate.modules.kpsc__amr.col_mutations import*


def resminimap_assembly(assembly, minimap2_index, ref_file, gene_info, qrdr, trunc, omp,  min_coverage, min_identity,
                          min_spurious_coverage, min_spurious_identity):
    hits_dict = minimap_against_all(assembly, minimap2_index, ref_file, gene_info, min_coverage, min_identity, min_spurious_coverage, min_spurious_identity)
    
    if qrdr:
        check_for_qrdr_mutations(hits_dict, assembly, qrdr, min_identity, 90.0)
        
    if trunc:
        check_for_mgrb_pmrb_gene_truncations(hits_dict, assembly, trunc, min_identity)
    if omp:
        check_omp_genes(hits_dict, assembly, omp, min_identity, 90.0)
    return hits_dict



def read_class_file(res_class_file):
    gene_info = {}  # key = sequence id (fasta header in ref file), value = (allele,class,Bla_Class)
    res_classes = []
    bla_classes = ['Bla', 'Bla_inhR', 'Bla_ESBL', 'Bla_ESBL_inhR', 'Bla_Carb', 'Bla_chr']

    with open(res_class_file, 'r') as f:
        header = 0
        for line in f:
            if header == 0:
                header = 1
                # clusterid,queryID,class,gene,allele,seqID,accession,positions,size,
                # cluster_contains_multiple_genes,gene_found_in_multiple_clusters,bla_description,
                # bla_class
            else:
                fields = line.rstrip().split(',')

                cluster_id, res_class, gene, allele_symbol, seq_id, bla_class = \
                    fields[0], fields[2], fields[3], fields[4], fields[5], fields[12]
                seq_header = '__'.join([cluster_id, gene + '_' + res_class, allele_symbol, seq_id])

                if res_class == 'Bla' and bla_class == 'NA':
                    bla_class = 'Bla'
                gene_info[seq_header] = (allele_symbol, res_class, bla_class)
                if res_class not in res_classes:
                    res_classes.append(res_class)
                if bla_class not in bla_classes:
                    bla_classes.append(bla_class)

    res_classes.sort()
    if 'Bla' in res_classes:
        res_classes.remove('Bla')
    if 'NA' in bla_classes:
        bla_classes.remove('NA')

    if 'SHV_mutations' not in res_classes:
        res_classes.append('SHV_mutations')
    if 'Omp_mutations' not in res_classes:
        res_classes.append('Omp_mutations')
    if 'Col_mutations' not in res_classes:
        res_classes.append('Col_mutations')
    if 'Flq_mutations' not in res_classes:
        res_classes.append('Flq_mutations')

    return gene_info, res_classes, bla_classes



def get_res_headers(res_classes, bla_classes):
    res_headers = res_classes + bla_classes

    # Rearrange the headers a bit. First move Bla_chr to the end:
    res_headers = ([h for h in res_headers if h != 'Bla_chr'] +
                   [h for h in res_headers if h == 'Bla_chr'])

    # Then move mutation columns to the end:
    res_headers = ([h for h in res_headers if '_mutations' not in h] +
                   [h for h in res_headers if '_mutations' in h])

    # Add '_acquired' to the end of the rest of the columns:
    res_headers = [h if h.endswith('_chr') or h.endswith('_mutations') else h + '_acquired'
                   for h in res_headers]

    return res_headers



def minimap_against_all(assembly, minimap2_index, ref_file, gene_info, min_coverage, min_identity, min_spurious_coverage, min_spurious_identity):
    """
    This function takes:
    * assembly:  assembly in FASTA format
    * ref_file: a path for a CARD reference in FASTA format
    * minimap2_index: a path for the assembly's minimap2 index (for faster alignment) (optional)
    * min_identity: hits with a lower percent identity than this are discarded
    
    This function returns:
    * dictionary with SHV mutations, truncated_resistance_hits, spurious_resistance_hits, _acquired mutations
    """
    hits_dict = collections.defaultdict(list) 

    alignment_hits = align_query_to_ref(ref_file, assembly, ref_index=minimap2_index, min_identity=min_identity, min_query_coverage=min_spurious_coverage)
    alignment_hits = cull_redundant_hits(alignment_hits)
    # print(alignment_hits)

    # calculate alignment coverage
    for hit in alignment_hits:
        alignment_length = hit.query_end - hit.query_start
        coverage = (alignment_length / hit.query_length) * 100
        
        if coverage >= min_spurious_coverage:
            if hit.percent_identity < 100.0:
                aa_result = check_for_exact_aa_match(ref_file, hit, assembly)
                if aa_result is not None:
                    hit.query_name = aa_result
                    exact_match = True
                else:
                    exact_match = False
            else:
                aa_result = None
                exact_match = True

            if hit.query_name not in gene_info:
                continue

            hit_allele, hit_class, hit_bla_class = gene_info[hit.query_name]
            
            hit_bla_class, shv_muts, class_changing_muts, omega_loop_seq = \
                    check_for_shv_mutations(hit, hit_allele, hit_bla_class, exact_match)
            
            if hit_class == 'Bla' and hit_bla_class:
                hit_class = hit_bla_class

            # SHV mutations
            for mut in shv_muts:
                mut_str, mut_metadata = mut
                mut_metadata['Genetic Variation Type'] = 'Protein variant detected'
                hits_dict['SHV_mutations'].append([mut_str, mut_metadata])

            if omega_loop_seq is not None:
                omega_str, omega_metadata = omega_loop_seq
                omega_metadata['Genetic Variation Type'] = 'Protein variant detected'
                hits_dict['SHV_mutations'].append([f'blaSHV:p.{omega_str}', omega_metadata])

                # hits_dict['SHV_mutations'].append([f'blaSHV:p.{omega_loop_seq}', {'Genetic Variation Type': 'Protein variant detected'}])
            if not hits_dict['SHV_mutations']:
                del hits_dict['SHV_mutations']
    
            #---- aac Logic ---
            # aac_pattern = r"^aac\(6'\)-Ib(?:-cr.*)?$"
            aac_pattern = r"^aac\(6'\)-Ib(?:-cr.*|\d.*)?$"
            is_aac = re.fullmatch(aac_pattern, hit.query_name.split('__')[2]) is not None
            target_classes = []

            if is_aac:
                # Check for Exact Match (100% Identity AND 100% Coverage)
                if hit.percent_identity >= 100.0 and coverage >= 100.0:
                    
                    # if Exact match and has "-cr" suffix 
                    if "-cr" in hit.query_name:
                        target_classes = ['AGly_acquired', 'Flq_acquired']
                    
                    # Exact match and no "-cr" suffix 
                    else:
                        target_classes = ['AGly_acquired']
                        
                # Inexact Match
                else:
                    aac_muts = check_for_aac_mutations(hit)
                    has_cr_muts = any("Trp102Arg" in m for m in aac_muts) and \
                                  any("Asp179Tyr" in m for m in aac_muts)
                    
                    if has_cr_muts:
                        target_classes = ['Flq_acquired', 'AGly_acquired']
                    else:
                        target_classes = ['AGly_acquired']

            else:
                if not (hit_class.endswith('_chr') or hit_class.endswith('_mutations')):
                    hit_class += '_acquired'
                target_classes = [hit_class]

            trunc_cov = 100.0
            if aa_result is not None:
                hit_allele += '^'
            else:
                if hit.percent_identity < 100:
                    hit_allele += '*'
                if alignment_length < hit.query_length:    
                    hit_allele += '?'
                trunc_suffix, trunc_cov, _ = truncation_check(hit)
                hit_allele += trunc_suffix
                
            if class_changing_muts:
                hit_allele += ' +' + ' +'.join(class_changing_muts)

            hit_data = {
                'Input Sequence ID': hit.ref_name,
                'Input Gene Length': hit.ref_length,
                'Input Gene Start': hit.ref_start +1,
                'Input Gene Stop': hit.ref_end,
                'Reference Gene Length': hit.query_length,
                'Reference Gene Start': hit.query_start +1,
                'Reference Gene Stop': hit.query_end,
                'Sequence Identity': f"{hit.percent_identity:.2f}",
                'Coverage': f"{coverage:.2f}",
                'Strand Orientation': hit.strand
            }

            for t_class in target_classes:
                if coverage >= min_coverage and hit.percent_identity >= min_identity and trunc_cov >= 90.0:
                    if t_class.endswith('_acquired') or t_class.endswith('_chr'):
                        hit_data['Genetic Variation Type'] = 'Gene presence detected'
                    hits_dict[t_class].append([hit_allele, hit_data])
                elif coverage >= min_coverage and hit.percent_identity >= min_identity and trunc_cov < 90.0:
                    hit_data['Genetic Variation Type'] = 'Gene presence detected'
                    hits_dict['truncated_resistance_hits'].append([hit_allele, hit_data])
                else:
                    hit_data['Genetic Variation Type'] = 'Gene presence detected'
                    hits_dict['spurious_resistance_hits'].append([hit_allele, hit_data])

    return hits_dict



def check_for_aac_mutations(hit):
    
    loci = [(102, 'W'), (179, 'D')]

    aa_map = {
        'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe', 'G': 'Gly',
        'H': 'His', 'I': 'Ile', 'K': 'Lys', 'L': 'Leu', 'M': 'Met', 'N': 'Asn',
        'P': 'Pro', 'Q': 'Gln', 'R': 'Arg', 'S': 'Ser', 'T': 'Thr', 'V': 'Val',
        'W': 'Trp', 'Y': 'Tyr'
    }

    # AAL93141.1 (aac(6’)-Ib ) Reference sequence
    aac_ref = 'VTNSNDSVTLRLMTEHDLAMLYEWLNRSHIVEWWGGEEARPTL' \
              'ADVQEQYLPSVLAQESVTPYIAMLNGEPIGYAQSYVALGSGDGWWEEETDPGVRGIDQLL' \
              'ANASQLGKGLGTKLVRALVELLFNDPEVTKIQTDPSPSNLRAIRCYEKAGFERQGTVTTP' \
              'DGPAVYMVQTRQAFERTRSVA'

    protein_aligner = Align.PairwiseAligner()
    protein_aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    protein_aligner.open_gap_score = -10
    protein_aligner.extend_gap_score = -0.5

    query_seq = hit.ref_seq
    query_translation = translate_nucl_to_prot(query_seq)

    if not query_translation:
        return []

    alignments = protein_aligner.align(aac_ref, query_translation)
    prt_alignment = alignments[0]


    mapping = get_mapping_by_query_pos(prt_alignment)

    mutations = []
    for pos, wt_base in loci:
        if pos in mapping:
            target_base, query_base = mapping[pos]
            
            if query_base != wt_base and query_base not in ['-', '.', '*']:
                wt_name = aa_map.get(wt_base, wt_base)
                mut_name = aa_map.get(query_base, query_base)
                mutation = f"{wt_name}{pos}{mut_name}"
                mutations.append(mutation)

    return mutations



def get_mapping_by_query_pos(alignment):
    """
    Maps alignment indices to query positions, starting at 16.
    """
    aligned_target = alignment[0]
    aligned_query = alignment[1]
    
    query_pos_map = {}
    query_counter = 15 
    
    for i in range(len(aligned_query)):
        t_base = aligned_target[i]
        q_base = aligned_query[i]
        
        if q_base != '-':
            query_counter += 1
            query_pos_map[query_counter] = (t_base, q_base)
            
    return query_pos_map
 
    

# def get_mapping_by_query_pos(alignment):
#     """
#     """
#     aligned_target = alignment[0]
#     aligned_query = alignment[1]
    
#     query_pos_map = {}
#     query_counter = 0
    
#     for i in range(len(aligned_query)):
#         t_base = aligned_target[i]
#         q_base = aligned_query[i]
        
#         if q_base != '-':
#             query_counter += 1
#             query_pos_map[query_counter] = (t_base, q_base)
            
#     return query_pos_map

