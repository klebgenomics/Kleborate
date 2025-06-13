"""
Copyright 2025 Mary Maranga
https://github.com/katholt/Kleborate/

This file is part of Kleborate. Kleborate is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Kleborate is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Kleborate. If
not, see <http://www.gnu.org/licenses/>.
"""

from Bio.Seq import Seq
from Bio import Align
from Bio.Align import substitution_matrices
from ...shared.alignment import align_query_to_ref, truncation_check, translate_nucl_to_prot, translate_nucl_to_prot, get_bases_per_ref_pos, get_last_base_aa_before_gap


def check_omp_genes(hits_dict, assembly, omp, min_identity, min_coverage):
    """
    Checks for OmpK35 and OmpK36 mutations.
    """

    best_ompk35_cov, best_ompk36_cov = 0.0, 0.0
    ompk35_hit_data, ompk36_hit_data = None, None
    ompk35_frameshift, ompk35_deletion = None, None
    ompk36_frameshift, ompk36_deletion = None, None

    ompk36_loci = {'OmpK36': [(25, 'C')]}

    aa_map = {
        'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe', 'G': 'Gly',
        'H': 'His', 'I': 'Ile', 'K': 'Lys', 'L': 'Leu', 'M': 'Met', 'N': 'Asn',
        'P': 'Pro', 'Q': 'Gln', 'R': 'Arg', 'S': 'Ser', 'T': 'Thr', 'V': 'Val',
        'W': 'Trp', 'Y': 'Tyr'
    }

    # define the aligner
    aligner = Align.PairwiseAligner(mode='global', match_score=5, mismatch_score=-4,
                                    open_gap_score=-10, extend_gap_score=-0.5)

    alignment_hits = align_query_to_ref(omp, assembly, min_query_coverage=None, min_identity=None)

    if 'Omp_mutations' not in hits_dict:
        hits_dict['Omp_mutations'] = []

    for hit in alignment_hits:
        _, coverage, translation = truncation_check(hit)
        dna_hit_cov = hit.query_cov

        hit_data = {
            'Input_sequence_ID': hit.ref_name,
            'Input_gene_length': hit.ref_length,
            'Input_gene_start': hit.ref_start,
            'Input_gene_stop': hit.ref_end,
            'Reference_gene_length': hit.query_length,
            'Reference_gene_start': hit.query_start,
            'Reference_gene_stop': hit.query_end,
            'Sequence_identity': f"{hit.percent_identity:.2f}%",
            'Coverage': f"{coverage:.2f}%",  # uses protein level coverage
            'Strand_orientation': hit.strand
        }

        if hit.query_name == 'OmpK35':

            if coverage > best_ompk35_cov:
                best_ompk35_cov = coverage
                ompk35_hit_data = hit_data
                ompk35_dna_cov = dna_hit_cov
                ompk35_ref_seq = hit.query_seq

                # ---- Frameshift/Deletion Logic for OmpK35 ----
                if best_ompk35_cov < 99.0 and ompk35_dna_cov > 99.0:
                    # checks for Frameshift
                    if translation:
                        ompk35_ref_trans = translate_nucl_to_prot(ompk35_ref_seq)
                        ompk35_prot_align = aligner.align(ompk35_ref_trans, translation)
                        aa_pos, aa_code = get_last_base_aa_before_gap(ompk35_prot_align[0])
                        fs_report = f"OmpK35:p.{aa_map[aa_code]}{aa_pos}fs"
                        ompk35_frameshift = (fs_report, {**hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'})

                elif best_ompk35_cov < 99.0 and ompk35_dna_cov < 99.0:
                    # checks for Deletion
                    ompk35_assembly_seq = hit.ref_seq
                    ompk35_dna_alignments = aligner.align(ompk35_ref_seq, ompk35_assembly_seq)
                    pos, base = get_last_base_aa_before_gap(ompk35_dna_alignments[0])
                    del_report = f"OmpK35:c.{base}{pos}del"
                    ompk35_deletion = (del_report, {**hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'})

        elif hit.query_name == 'OmpK36':
            if coverage > best_ompk36_cov:
                best_ompk36_cov = coverage
                ompk36_hit_data = hit_data
                ompk36_dna_cov = dna_hit_cov

                ompk36_ref_seq = hit.query_seq
                ompk36_assembly_seq = hit.ref_seq
                ompk36_dna_alignments = aligner.align(ompk36_ref_seq, ompk36_assembly_seq)

                # ---- Frameshift/Deletion Logic for OmpK36 ----
                if best_ompk36_cov < 99.0 and ompk36_dna_cov > 99.0:
                    # checks for Frameshift
                    if translation:
                        ompk36_ref_trans = translate_nucl_to_prot(ompk36_ref_seq)
                        ompk36_prot_align = aligner.align(ompk36_ref_trans, translation)
                        aa_pos, aa_code = get_last_base_aa_before_gap(ompk36_prot_align[0])
                        fs_report = f"OmpK36:p.{aa_map[aa_code]}{aa_pos}fs"
                        ompk36_frameshift = (fs_report, {**hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'})

                elif best_ompk36_cov < 99.0 and ompk36_dna_cov < 99.0:
                    # checks for Deletion
                    pos, base = get_last_base_aa_before_gap(ompk36_dna_alignments[0])
                    del_report = f"OmpK36:c.{base}{pos}del"
                    ompk36_deletion = (del_report, {**hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'})

                # ----- NUCLEOTIDE VARIANT -----
                bases_per_ref_pos = get_bases_per_ref_pos(ompk36_dna_alignments[0])
                loci = ompk36_loci[hit.query_name]
                for pos, wt_base in loci:
                    assembly_base = bases_per_ref_pos[pos]
                    ref_base = ompk36_ref_seq[pos - 1]
                    if ref_base == wt_base and assembly_base == 'T':
                        hits_dict['Omp_mutations'].append([
                            f"{hit.query_name}:c.{pos}{wt_base.upper()}>{assembly_base.upper()}",
                            {'Genetic_variation_type': 'Nucleotide variant detected', **hit_data}
                        ])

                # ----- PROTEIN VARIANT -----
                if coverage >= min_coverage:
                    if 'GDGDTY' in translation:
                        l3_insertion, insertion_type = 'GDGDTY', 'GD'
                    elif 'GDTDTY' in translation:
                        l3_insertion, insertion_type = 'GDTDTY', 'TD'
                    else:
                        l3_insertion, insertion_type = None, None

                    if l3_insertion:
                        motif_start_index = translation.index(l3_insertion)
                        insertion_relative_index = l3_insertion.index(insertion_type)

                        insertion_start_aa = motif_start_index + insertion_relative_index + 1
                        insertion_end_aa = insertion_start_aa + 1

                        inserted_residues = f"{aa_map[insertion_type[0]]}{aa_map[insertion_type[1]]}"
                        insertion_annotation = f"OmpK36:p.{insertion_start_aa}_{insertion_end_aa}ins{inserted_residues}"

                        hits_dict['Omp_mutations'].append([
                            insertion_annotation,
                            {'Genetic_variation_type': 'Protein variant detected', **hit_data}
                        ])
        else:
            assert False

    # reporting truncations
    truncations = []
    if ompk35_frameshift:
        truncations.append(ompk35_frameshift)
    elif ompk35_deletion:
        truncations.append(ompk35_deletion)

    if ompk36_frameshift:
        truncations.append(ompk36_frameshift)
    elif ompk36_deletion:
        truncations.append(ompk36_deletion)

    for trunc_name, hit_data in truncations:
        data = dict(hit_data) if hit_data else {}
        hits_dict['Omp_mutations'].append([trunc_name, data])



# def check_omp_genes(hits_dict, assembly, omp, min_identity, min_coverage):

#     best_ompk35_cov, best_ompk36_cov = 0.0, 0.0
#     ompk36_loci = {'OmpK36': [(25, 'C')]}
    
#     # define the aligner
#     aligner = Align.PairwiseAligner(mode='global', match_score=5, mismatch_score=-4, open_gap_score = -10, extend_gap_score = -0.5)
    
#     alignment_hits = align_query_to_ref(omp, assembly, min_query_coverage=None, min_identity=None)
    
#     ompk35_hit = False
#     ompk36_hit = False
    
#     for hit in alignment_hits:
#         _, coverage, translation = truncation_check(hit)
        
#         if hit.query_name == 'OmpK35':
#             ompk35_hit = True
#             if coverage > best_ompk35_cov:
#                 best_ompk35_cov = coverage

#         elif hit.query_name == 'OmpK36':
#             ompk36_hit = True
#             query_seq = hit.query_seq
#             assembly_seq = hit.ref_seq
#             alignments = aligner.align(query_seq , assembly_seq)
#             bases_per_ref_pos = get_bases_per_ref_pos(alignments[0])
#             loci = ompk36_loci[hit.query_name]
#             for pos, wt_base in loci:
#                 assembly_base = bases_per_ref_pos[pos]
#                 query_base = query_seq[pos-1]
#                 if query_base == wt_base and assembly_base == 'T':
#                     hits_dict['Omp_mutations'].append(f"{hit.query_name}_{wt_base.lower()}{pos}{assembly_base.lower()}")

#             if coverage > best_ompk36_cov:
#                 best_ompk36_cov = coverage

#             if coverage >= min_coverage:
#                 if 'GDGDTY' in translation:
#                     hits_dict['Omp_mutations'].append('OmpK36GD')
                    
#                 elif 'GDTDTY' in translation:
#                     hits_dict['Omp_mutations'].append('OmpK36TD')
                
#         else:
#             assert False

#     truncations = []
#     if ompk35_hit and best_ompk35_cov < 90.0:
#         truncations.append('OmpK35-' + ('%.0f' % best_ompk35_cov) + '%')
   
#     if ompk36_hit and best_ompk36_cov < 90.0:
#         truncations.append('OmpK36-' + ('%.0f' % best_ompk36_cov) + '%')
    

#     if truncations:
#         if 'Omp_mutations' not in hits_dict:
#             hits_dict['Omp_mutations'] = []
#         hits_dict['Omp_mutations'] += truncations


