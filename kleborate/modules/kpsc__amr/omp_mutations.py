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
import traceback
from Bio import Align
from Bio.Align import substitution_matrices
from ...shared.alignment import align_query_to_ref, truncation_check, translate_nucl_to_prot, get_bases_per_ref_pos, find_start_deletion_in_alignment, deletion_checks, get_frameshift_info
from ...shared.misc import load_fasta, reverse_complement


def check_omp_genes(hits_dict, assembly, omp, min_identity, min_coverage):
    """
    Checks for OmpK35 and OmpK36 mutations.
    """
    best_ompk35_cov, best_ompk36_cov = 0.0, 0.0
    
    ompk35_hit_data, ompk36_hit_data = None, None
    ompk35_frameshift, ompk35_deletion, ompk35_insertion = None, None, None
    ompk36_frameshift, ompk36_deletion, ompk36_insertion = None, None, None

    ompk36_loci = {'OmpK36': [(25, 'C')]}
    ref_seqs = dict(load_fasta(omp))

    aa_map = {'*': 'Ter'}

    # Protein aligner
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.mode = 'global'
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    # DNA aligner
    dna_aligner = Align.PairwiseAligner(
        mode='global',
        # match_score=5,
        # mismatch_score=-10,
        open_gap_score=-10,
        extend_gap_score=-0.5
    )

    alignment_hits = align_query_to_ref(omp, assembly, min_query_coverage=None, min_identity=None)

    if 'Omp_mutations' not in hits_dict:
        hits_dict['Omp_mutations'] = []

    for hit in alignment_hits:
        _, coverage, translation = truncation_check(hit)
        dna_hit_cov = hit.query_cov

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

        # --- Frameshift and Deletion checks in OmpK35 ---
        if hit.query_name == 'OmpK35':
            if coverage >= best_ompk35_cov:
                best_ompk35_cov = coverage
                ompk35_dna_cov = dna_hit_cov
                ompk35_hit_data = hit_data
                ompk35_ref_seq = ref_seqs['OmpK35']
                ompk35_query_seq = hit.ref_seq

                
                ompk35_frameshift, ompk35_deletion = None, None

                if coverage > 110.0:
                    ompk35_aln = dna_aligner.align(ompk35_ref_seq, ompk35_query_seq)
                    insertion_pos = insertion_checks(ompk35_aln)
                    ompK35_insertion = (
                        f"ompK35:p.{insertion_pos}",
                        {**hit_data, 'Genetic Variation Type': 'Inactivating mutation detected'}
                    )


                elif coverage == 0.0:
                    aln = dna_aligner.align(ompk35_ref_seq, ompk35_query_seq)[0]
                    deleted_base_pos = find_start_deletion_in_alignment(aln)
                    ompk35_deletion = (
                        f"ompK35:{deleted_base_pos}",
                        {**hit_data, 'Genetic Variation Type': 'Inactivating mutation detected'}
                    )
                
                elif best_ompk35_cov < 90.0 and ompk35_dna_cov > 90.0:
                    if translation:
                        ompk35_ref_trans = translate_nucl_to_prot(ompk35_ref_seq)
                        ompk35_query_trans = translate_nucl_to_prot(ompk35_query_seq)
                        ompk35_prot_align = aligner.align(ompk35_ref_trans, ompk35_query_trans)
                        fs_info = get_frameshift_info(ompk35_prot_align[0])
                        if fs_info is not None:
                            aa_pos, ref_aa, alt_aa, fs_len = fs_info
                            alt_str = aa_map.get(alt_aa, alt_aa)
                            fs_report = f"ompK35:p.{ref_aa}{aa_pos}{alt_str}" if fs_len == 0 else f"ompK35:p.{ref_aa}{aa_pos}{alt_str}fsTer{fs_len}"
                            ompk35_frameshift = (fs_report, {**hit_data, 'Genetic Variation Type': 'Inactivating mutation detected'})
                
                elif best_ompk35_cov < 90.0 and ompk35_dna_cov < 90.0:
                    ompk35_dna_alignments = dna_aligner.align(ompk35_ref_seq, ompk35_query_seq)
                    del_info = deletion_checks(ompk35_dna_alignments[0], ompk35_ref_seq)
                    if del_info is not None:
                        pos, base = del_info
                        ompk35_deletion = (f"ompK35:c.{base}{pos}del", {**hit_data, 'Genetic Variation Type': 'Inactivating mutation detected'})

        # --- Frameshift and Deletion checks in OmpK36 ---
        elif hit.query_name == 'OmpK36':
            if coverage >= best_ompk36_cov:
                best_ompk36_cov = coverage
                ompk36_dna_cov = dna_hit_cov
                ompk36_hit_data = hit_data
                ompk36_ref_seq = ref_seqs['OmpK36']
                ompk36_query_seq = hit.ref_seq

                ompk36_frameshift, ompk36_deletion = None, None
                ompk36_dna_alignments = dna_aligner.align(ompk36_ref_seq, ompk36_query_seq)

                if coverage > 110.0:
                    insertion_pos = insertion_checks(ompk36_dna_alignments )
                    ompK36_insertion = (
                        f"ompK35:p.{insertion_pos}",
                        {**hit_data, 'Genetic Variation Type': 'Inactivating mutation detected'}
                    )


                elif coverage == 0.0:
                    aln = ompk36_dna_alignments[0]
                    deleted_base_pos = find_start_deletion_in_alignment(aln)
                    ompk36_deletion = (
                        f"ompK36:{deleted_base_pos}",
                        {**hit_data, 'Genetic Variation Type': 'Inactivating mutation detected'}
                    )
                
                elif best_ompk36_cov < 90.0 and ompk36_dna_cov > 90.0:
                    if translation:
                        ompk36_ref_trans = translate_nucl_to_prot(ompk36_ref_seq)
                        ompk36_query_trans = translate_nucl_to_prot(ompk36_query_seq)
                        ompk36_prot_align = aligner.align(ompk36_ref_trans, ompk36_query_trans)
                        fs_info = get_frameshift_info(ompk36_prot_align[0])
                        if fs_info is not None:
                            aa_pos, ref_aa, alt_aa, fs_len = fs_info
                            alt_str = aa_map.get(alt_aa, alt_aa)
                            fs_report = f"ompK36:p.{ref_aa}{aa_pos}{alt_str}" if fs_len == 0 else f"ompK36:p.{ref_aa}{aa_pos}{alt_str}fsTer{fs_len}"
                            ompk36_frameshift = (fs_report, {**hit_data, 'Genetic Variation Type': 'Inactivating mutation detected'})
                
                elif best_ompk36_cov < 90.0 and ompk36_dna_cov < 90.0:
                    del_info = deletion_checks(ompk36_dna_alignments[0], ompk36_ref_seq)
                    if del_info is not None:
                        pos, base = del_info
                        ompk36_deletion = (f"ompK36:c.{base}{pos}del", {**hit_data, 'Genetic Variation Type': 'Inactivating mutation detected'})

                # ----- OmpK36 Nucleotide variant -----
                bases_per_ref_pos = get_bases_per_ref_pos(ompk36_dna_alignments[0])
                loci = ompk36_loci.get(hit.query_name, [])
                for pos, wt_base in loci:
                    assembly_base = bases_per_ref_pos[pos]
                    ref_base = ompk36_ref_seq[pos - 1]
                    if ref_base == wt_base and assembly_base == 'T':

                        hits_dict['Omp_mutations'].append([
                            f"{hit.query_name[0].lower() + hit.query_name[1:]}:c.{pos}{wt_base.upper()}>{assembly_base.upper()}",
                            {'Genetic Variation Type': 'Nucleotide variant detected', **hit_data}
                        ])
                        # hits_dict['Omp_mutations'].append([
                        #     f"{hit.query_name}:c.{pos}{wt_base.upper()}>{assembly_base.upper()}",
                        #     {'Genetic Variation Type': 'Nucleotide variant detected', **hit_data}
                        # ])
                # ----- OmpK36 L3 insertion checks -----
                if coverage >= min_coverage and translation:
                    l3_insertion, insertion_type = None, None
                    if 'GDGDTY' in translation:
                        l3_insertion, insertion_type = 'GDGDTY', 'GD'
                    elif 'GDTDTY' in translation:
                        l3_insertion, insertion_type = 'GDTDTY', 'TD'
                    elif 'GDDTY' in translation:
                        l3_insertion, insertion_type = 'GDDTY', 'D'

                    if l3_insertion:
                        motif_start_index = translation.index(l3_insertion)
                        insertion_relative_index = l3_insertion.index(insertion_type)
                        insertion_start_aa = motif_start_index + insertion_relative_index + 1
                        insertion_end_aa = insertion_start_aa + 1
                        insertion_annotation = f"ompK36:p.{insertion_start_aa}_{insertion_end_aa}ins{insertion_type}"
                        
                        hits_dict['Omp_mutations'].append([
                            insertion_annotation,
                            {'Genetic Variation Type': 'Inactivating mutation detected', **hit_data}
                        ])

    truncations = []
    if ompk35_insertion:
        truncations.append(ompk35_insertion)
    elif ompk35_frameshift:
        truncations.append(ompk35_frameshift)
    elif ompk35_deletion:
        truncations.append(ompk35_deletion)

    if ompk36_insertion:
        truncations.append(ompk36_insertion)
    elif ompk36_frameshift:
        truncations.append(ompk36_frameshift)
    elif ompk36_deletion:
        truncations.append(ompk36_deletion)

    for trunc_name, t_hit_data in truncations:
        data = dict(t_hit_data) if t_hit_data else {}
        hits_dict['Omp_mutations'].append([trunc_name, data])

    # --- check if OmpK35 and OmpK36 gene is deleted ---
    if ompk35_hit_data is None:
        hits_dict['Omp_mutations'].append([
            "ompK35:del", 
            {"Genetic Variation Type": "Gene deletion detected", "Coverage": "0.00"}
        ])

    if ompk36_hit_data is None:
        hits_dict['Omp_mutations'].append([
            "ompK36:del", 
            {"Genetic Variation Type": "Gene deletion detected", "Coverage": "0.00"}
        ])





def insertion_checks(alignment):
    """
    Identifies insertion in the alignment 
    """
    aligned_ref, aligned_query = alignment[0], alignment[1]

    insertions = [] 
    ref_pos = 0
    insertion_len = 0
    insertion_start = None

    # Scan the alignment for insertions
    for ref_aln, query_aln in zip(aligned_ref, aligned_query):
        if ref_aln != '-':
            ref_pos += 1

        if ref_aln == '-' and query_aln != '-':
            if insertion_len == 0:
                insertion_start = ref_pos
            insertion_len += 1
        else:
            if insertion_len > 0:
                insertions.append((insertion_start, insertion_len))
                insertion_len = 0
                insertion_start = None

    if insertion_len > 0:
        insertions.append((insertion_start, insertion_len))

    if not insertions:
        return None
    
    insertion_start_pos, insertion_length = max(insertions, key=lambda x: x[1])
    
    insertion_end_pos = insertion_start_pos + 1

    return f"{insertion_start_pos}_{insertion_end_pos}insN[{insertion_length}]"




# def check_omp_genes(hits_dict, assembly, omp, min_identity, min_coverage):
#     """
#     Checks for OmpK35 and OmpK36 mutations.
#     """
#     best_ompk35_cov, best_ompk36_cov = 0.0, 0.0
    
#     ompk35_hit_data, ompk36_hit_data = None, None
#     ompk35_frameshift, ompk35_deletion = None, None
#     ompk36_frameshift, ompk36_deletion = None, None

#     ompk36_loci = {'OmpK36': [(25, 'C')]}
#     ref_seqs = dict(load_fasta(omp))

#     aa_map = {'*': 'Ter'}

#     # Protein aligner
#     aligner = Align.PairwiseAligner()
#     aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
#     aligner.mode = 'global'
#     aligner.open_gap_score = -10
#     aligner.extend_gap_score = -0.5

#     # DNA aligner
#     dna_aligner = Align.PairwiseAligner(
#         mode='global',
#         match_score=5,
#         mismatch_score=-10,
#         open_gap_score=-10,
#         extend_gap_score=-0.5
#     )

#     alignment_hits = align_query_to_ref(omp, assembly, min_query_coverage=None, min_identity=None)
#     # print(alignment_hits)

#     if 'Omp_mutations' not in hits_dict:
#         hits_dict['Omp_mutations'] = []

#     for hit in alignment_hits:
#         _, coverage, translation = truncation_check(hit)
#         # print(f"{hit.query_name} coverage: {coverage}")
#         # print(f"{hit.query_name} seq: {hit.ref_seq}")
#         dna_hit_cov = hit.query_cov

#         hit_data = {
#             'Input_sequence_ID': hit.ref_name,
#             'Input_gene_length': hit.ref_length,
#             'Input_gene_start': hit.ref_start,
#             'Input_gene_stop': hit.ref_end,
#             'Reference_gene_length': hit.query_length,
#             'Reference_gene_start': hit.query_start,
#             'Reference_gene_stop': hit.query_end,
#             'Sequence_identity': f"{hit.percent_identity:.2f}%",
#             'Coverage': f"{coverage:.2f}%",
#             'Strand_orientation': hit.strand
#         }

#         # --- Frameshift and Deletion checks in OmpK35 ---
#         if hit.query_name == 'OmpK35':
#             if coverage >= best_ompk35_cov:
#                 best_ompk35_cov = coverage
#                 ompk35_dna_cov = dna_hit_cov
#                 ompk35_hit_data = hit_data
#                 ompk35_ref_seq = ref_seqs['OmpK35']
#                 ompk35_query_seq = hit.ref_seq

                
#                 ompk35_frameshift, ompk35_deletion = None, None

#                 if coverage == 0.0:
#                     aln = dna_aligner.align(ompk35_ref_seq, ompk35_query_seq)[0]
#                     deleted_base_pos = find_start_deletion_in_alignment(aln)
#                     ompk35_deletion = (
#                         f"OmpK35:{deleted_base_pos}",
#                         {**hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'}
#                     )
                
#                 elif best_ompk35_cov < 90.0 and ompk35_dna_cov > 90.0:
#                     if translation:
#                         ompk35_ref_trans = translate_nucl_to_prot(ompk35_ref_seq)
#                         ompk35_query_trans = translate_nucl_to_prot(ompk35_query_seq)
#                         ompk35_prot_align = aligner.align(ompk35_ref_trans, ompk35_query_trans)
#                         fs_info = get_frameshift_info(ompk35_prot_align[0])
#                         if fs_info is not None:
#                             aa_pos, ref_aa, alt_aa, fs_len = fs_info
#                             alt_str = aa_map.get(alt_aa, alt_aa)
#                             fs_report = f"OmpK35:p.{ref_aa}{aa_pos}{alt_str}" if fs_len == 0 else f"OmpK35:p.{ref_aa}{aa_pos}{alt_str}fsTer{fs_len}"
#                             ompk35_frameshift = (fs_report, {**hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'})
                
#                 elif best_ompk35_cov < 90.0 and ompk35_dna_cov < 90.0:
#                     ompk35_dna_alignments = dna_aligner.align(ompk35_ref_seq, ompk35_query_seq)
#                     # print(ompk35_dna_alignments[0])
#                     del_info = deletion_checks(ompk35_dna_alignments[0], ompk35_ref_seq)
#                     if del_info is not None:
#                         pos, base = del_info
#                         ompk35_deletion = (f"OmpK35:c.{base}{pos}del", {**hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'})

#         # --- Frameshift and Deletion checks in OmpK36 ---
#         elif hit.query_name == 'OmpK36':
#             if coverage >= best_ompk36_cov:
#                 best_ompk36_cov = coverage
#                 ompk36_dna_cov = dna_hit_cov
#                 ompk36_hit_data = hit_data
#                 ompk36_ref_seq = ref_seqs['OmpK36']
#                 ompk36_query_seq = hit.ref_seq

#                 ompk36_frameshift, ompk36_deletion = None, None
#                 ompk36_dna_alignments = dna_aligner.align(ompk36_ref_seq, ompk36_query_seq)

#                 if coverage == 0.0:
#                     aln = ompk36_dna_alignments[0]
#                     deleted_base_pos = find_start_deletion_in_alignment(aln)
#                     ompk36_deletion = (
#                         f"OmpK36:{deleted_base_pos}",
#                         {**hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'}
#                     )
                
#                 elif best_ompk36_cov < 90.0 and ompk36_dna_cov > 90.0:
#                     if translation:
#                         ompk36_ref_trans = translate_nucl_to_prot(ompk36_ref_seq)
#                         ompk36_query_trans = translate_nucl_to_prot(ompk36_query_seq)
#                         ompk36_prot_align = aligner.align(ompk36_ref_trans, ompk36_query_trans)
#                         fs_info = get_frameshift_info(ompk36_prot_align[0])
#                         if fs_info is not None:
#                             aa_pos, ref_aa, alt_aa, fs_len = fs_info
#                             alt_str = aa_map.get(alt_aa, alt_aa)
#                             fs_report = f"OmpK36:p.{ref_aa}{aa_pos}{alt_str}" if fs_len == 0 else f"OmpK36:p.{ref_aa}{aa_pos}{alt_str}fsTer{fs_len}"
#                             ompk36_frameshift = (fs_report, {**hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'})
                
#                 elif best_ompk36_cov < 90.0 and ompk36_dna_cov < 90.0:
#                     del_info = deletion_checks(ompk36_dna_alignments[0], ompk36_ref_seq)
#                     if del_info is not None:
#                         pos, base = del_info
#                         ompk36_deletion = (f"OmpK36:c.{base}{pos}del", {**hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'})

#                 # ----- OmpK36 Nucleotide variant -----
#                 bases_per_ref_pos = get_bases_per_ref_pos(ompk36_dna_alignments[0])
#                 loci = ompk36_loci.get(hit.query_name, [])
#                 for pos, wt_base in loci:
#                     assembly_base = bases_per_ref_pos[pos]
#                     ref_base = ompk36_ref_seq[pos - 1]
#                     if ref_base == wt_base and assembly_base == 'T':
#                         hits_dict['Omp_mutations'].append([
#                             f"{hit.query_name}:c.{pos}{wt_base.upper()}>{assembly_base.upper()}",
#                             {'Genetic_variation_type': 'Nucleotide variant detected', **hit_data}
#                         ])

#                 # ----- OmpK36 L3 insertion checks -----
#                 if coverage >= min_coverage and translation:
#                     l3_insertion, insertion_type = None, None
#                     if 'GDGDTY' in translation:
#                         l3_insertion, insertion_type = 'GDGDTY', 'GD'
#                     elif 'GDTDTY' in translation:
#                         l3_insertion, insertion_type = 'GDTDTY', 'TD'
#                     elif 'GDDTY' in translation:
#                         l3_insertion, insertion_type = 'GDDTY', 'D'

#                     if l3_insertion:
#                         motif_start_index = translation.index(l3_insertion)
#                         insertion_relative_index = l3_insertion.index(insertion_type)
#                         insertion_start_aa = motif_start_index + insertion_relative_index + 1
#                         insertion_end_aa = insertion_start_aa + 1
#                         insertion_annotation = f"OmpK36:p.{insertion_start_aa}_{insertion_end_aa}ins{insertion_type}"
                        
#                         hits_dict['Omp_mutations'].append([
#                             insertion_annotation,
#                             {'Genetic_variation_type': 'Protein variant detected', **hit_data}
#                         ])

#     truncations = []
#     if ompk35_frameshift:
#         truncations.append(ompk35_frameshift)
#     elif ompk35_deletion:
#         truncations.append(ompk35_deletion)

#     if ompk36_frameshift:
#         truncations.append(ompk36_frameshift)
#     elif ompk36_deletion:
#         truncations.append(ompk36_deletion)

#     for trunc_name, t_hit_data in truncations:
#         data = dict(t_hit_data) if t_hit_data else {}
#         hits_dict['Omp_mutations'].append([trunc_name, data])

#     # --- check if OmpK35 and OmpK36 gene is deleted ---
#     if ompk35_hit_data is None:
#         hits_dict['Omp_mutations'].append([
#             "OmpK35:del", 
#             {"Genetic_variation_type": "Gene deleted", "Coverage": "0.00%"}
#         ])

#     if ompk36_hit_data is None:
#         hits_dict['Omp_mutations'].append([
#             "OmpK36:del", 
#             {"Genetic_variation_type": "Gene deleted", "Coverage": "0.00%"}
#         ])


