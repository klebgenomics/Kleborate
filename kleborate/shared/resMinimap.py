"""
Search for resistance genes, summarise by class (one class per column)

Copyright 2023 Kat Holt, Ryan Wick (rrwick@gmail.com), Mary Maranga (gathonimaranga@gmail.com)
https://github.com/klebgenomics/KleborateModular/

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
 
from .alignment import align_query_to_ref, cull_redundant_hits, is_exact_aa_match, translate_nucl_to_prot, check_for_exact_aa_match, truncation_check
from .misc import load_fasta, reverse_complement
from kleborate.modules.klebsiella_pneumo_complex__amr.shv_mutations import*
from kleborate.modules.klebsiella_pneumo_complex__amr.qrdr_mutations import*
from kleborate.modules.klebsiella_pneumo_complex__amr.omp_mutations import*
from kleborate.modules.klebsiella_pneumo_complex__amr.col_mutations import*


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
    
    hits_dict = collections.defaultdict(list)  # key = class, value = list
    alignment_hits = align_query_to_ref(ref_file, assembly,ref_index=minimap2_index,  min_identity=min_identity, min_query_coverage=min_spurious_coverage)
    alignment_hits = cull_redundant_hits(alignment_hits)
    
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
                

            hit_allele, hit_class, hit_bla_class = gene_info[hit.query_name]
            

            hit_bla_class, shv_muts, class_changing_muts, omega_loop_seq = \
                    check_for_shv_mutations(hit, hit_allele, hit_bla_class, exact_match)
            
            # checks if the variable hit_class contains the string value 'Bla'
            # If it does, then the value of hit_class is replaced with the value of hit_bla_class.
            if hit_class == 'Bla':
                hit_class = hit_bla_class

            # append the list shv_muts to hits_dict['SHV_mutations']
            hits_dict['SHV_mutations'] += shv_muts
            if omega_loop_seq is not None:
                hits_dict['SHV_mutations'].append(f'omega-loop={omega_loop_seq}')

            # checks if the list associated with the 'SHV_mutations' key in the hits_dict dictionary is empty
            # If it is, the key 'SHV_mutations' and its associated value are deleted from the dictionary.
            if not hits_dict['SHV_mutations']:
                del hits_dict['SHV_mutations']

            # checks if the value in hit_class does not end with the strings '_chr' or '_mutations'.
            # append '_acquired'  to hit_class
            if not (hit_class.endswith('_chr') or hit_class.endswith('_mutations')):
                hit_class += '_acquired'

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
            
            
            # If the hit is decent (above the min coverage and identity thresholds), it goes in the
            # column for the class.
            if coverage >= min_coverage and hit.percent_identity >= min_identity and trunc_cov >= 90.0:
                hits_dict[hit_class].append(hit_allele)
                
            # If the hit is decent but the gene is truncated, it goes in the
            # truncated_resistance_hits column.
            elif coverage >= min_coverage and hit.percent_identity >= min_identity and trunc_cov < 90.0:
                hits_dict['truncated_resistance_hits'].append(hit_allele)
                
            # If the hit is bad (below the min coverage and identity thresholds but above the
            # thresholds for spurious hits) then it goes in the spurious hit column.
            else:
                hits_dict['spurious_resistance_hits'].append(hit_allele)
    
    return hits_dict
