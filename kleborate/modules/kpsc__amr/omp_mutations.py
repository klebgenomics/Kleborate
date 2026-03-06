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
from ...shared.alignment import align_query_to_ref, truncation_check, translate_nucl_to_prot, translate_nucl_to_prot, get_bases_per_ref_pos, find_start_deletion_in_alignment, deletion_checks, get_frameshift_info
from ...shared.misc import load_fasta, reverse_complement


def check_omp_genes(hits_dict, assembly, omp, min_identity, min_coverage):
    """
    Checks for OmpK35 and OmpK36 mutations.
    """
    best_ompk35_cov, best_ompk36_cov = 0.0, 0.0
    ompk35_hit_data, ompk36_hit_data = None, None
    ompk35_frameshift, ompk35_deletion = None, None
    ompk36_frameshift, ompk36_deletion = None, None

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
        match_score=5,
        mismatch_score=-10,
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
            'Input_sequence_ID': hit.ref_name,
            'Input_gene_length': hit.ref_length,
            'Input_gene_start': hit.ref_start,
            'Input_gene_stop': hit.ref_end,
            'Reference_gene_length': hit.query_length,
            'Reference_gene_start': hit.query_start,
            'Reference_gene_stop': hit.query_end,
            'Sequence_identity': f"{hit.percent_identity:.2f}%",
            'Coverage': f"{coverage:.2f}%",
            'Strand_orientation': hit.strand
        }

        # --- Frameshift and Deletion checks in OmpK35---
        if hit.query_name == 'OmpK35':
            ompk35_ref_seq = ref_seqs['OmpK35']
            ompk35_query_seq = hit.ref_seq
            ompk35_hit_data = hit_data

            if coverage == 0.0:
                aln = dna_aligner.align(ompk35_ref_seq, ompk35_query_seq)[0]
                deleted_base_pos = find_start_deletion_in_alignment(aln)
                deletion_report = f"OmpK35_{deleted_base_pos}"

                ompk35_deletion = (
                    deletion_report,
                    {**hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'}
                )
                continue  

            if coverage > best_ompk35_cov:
                best_ompk35_cov = coverage
                ompk35_dna_cov = dna_hit_cov

                if best_ompk35_cov < 90.0 and ompk35_dna_cov > 90.0:
                    if translation:
                        ompk35_ref_trans = translate_nucl_to_prot(ompk35_ref_seq)
                        ompk35_query_trans = translate_nucl_to_prot(ompk35_query_seq)
                        ompk35_prot_align = aligner.align(ompk35_ref_trans, ompk35_query_trans)

                        fs_info = get_frameshift_info(ompk35_prot_align[0])
                        if fs_info is not None:
                            aa_pos, ref_aa, alt_aa, fs_len = fs_info
                            alt_str = aa_map.get(alt_aa, alt_aa)
                            if fs_len == 0:
                                fs_report = f"OmpK35_{ref_aa}{aa_pos}{alt_str}"
                            else:
                                fs_report = f"OmpK35_{ref_aa}{aa_pos}{alt_str}fsTer{fs_len}"

                            ompk35_frameshift = (
                                fs_report,
                                {**hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'}
                            )
                    
                elif best_ompk35_cov < 90.0 and ompk35_dna_cov < 90.0:
                    ompk35_dna_alignments = dna_aligner.align(ompk35_ref_seq, ompk35_query_seq)
                    del_info = deletion_checks(ompk35_dna_alignments[0], ompk35_ref_seq)
                    if del_info is not None:
                        pos, base = del_info
                        del_report = f"OmpK35:c.{base}{pos}del"
                        ompk35_deletion = (
                            del_report,
                            {**hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'}
                        )

        # --- Frameshift and Deletion checks in OmpK36 ---
        elif hit.query_name == 'OmpK36':
            ompk36_ref_seq = ref_seqs['OmpK36']
            ompk36_query_seq = hit.ref_seq
            ompk36_hit_data = hit_data

            if coverage == 0.0:
                aln = dna_aligner.align(ompk36_ref_seq, ompk36_query_seq)[0]
                deleted_base_pos = find_start_deletion_in_alignment(aln)
                deletion_report = f"OmpK36_{deleted_base_pos}"

                ompk36_deletion = (
                    deletion_report,
                    {**hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'}
                )
                continue 

            if coverage > best_ompk36_cov:
                best_ompk36_cov = coverage
                ompk36_dna_cov = dna_hit_cov
                ompk36_dna_alignments = dna_aligner.align(ompk36_ref_seq, ompk36_query_seq)

                if best_ompk36_cov < 90.0 and ompk36_dna_cov > 90.0:
                    if translation:
                        ompk36_ref_trans = translate_nucl_to_prot(ompk36_ref_seq)
                        ompk36_query_trans = translate_nucl_to_prot(ompk36_query_seq)
                        ompk36_prot_align = aligner.align(ompk36_ref_trans, ompk36_query_trans)

                        fs_info = get_frameshift_info(ompk36_prot_align[0])
                        if fs_info is not None:
                            aa_pos, ref_aa, alt_aa, fs_len = fs_info
                            alt_str = aa_map.get(alt_aa, alt_aa)
                            if fs_len == 0:
                                fs_report = f"OmpK36_{ref_aa}{aa_pos}{alt_str}"
                            else:
                                fs_report = f"OmpK36_{ref_aa}{aa_pos}{alt_str}fsTer{fs_len}"

                            ompk36_frameshift = (
                                fs_report,
                                {**hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'}
                            )
                        
                elif best_ompk36_cov < 90.0 and ompk36_dna_cov < 90.0:
                    del_info = deletion_checks(ompk36_dna_alignments[0], ompk36_ref_seq)
                    if del_info is not None:
                        pos, base = del_info
                        del_report = f"OmpK36:c.{base}{pos}del"
                        ompk36_deletion = (
                            del_report,
                            {**hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'}
                        )

                # ----- OmpK36 Nucleotide variant -----
                bases_per_ref_pos = get_bases_per_ref_pos(ompk36_dna_alignments[0])
                loci = ompk36_loci.get(hit.query_name, [])
                for pos, wt_base in loci:
                    assembly_base = bases_per_ref_pos[pos]
                    ref_base = ompk36_ref_seq[pos - 1]
                    if ref_base == wt_base and assembly_base == 'T':
                        hits_dict['Omp_mutations'].append([
                            f"{hit.query_name}:c.{pos}{wt_base.upper()}>{assembly_base.upper()}",
                            {'Genetic_variation_type': 'Nucleotide variant detected', **hit_data}
                        ])

                # L3 insertion checks
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
                        inserted_residues = insertion_type
                        insertion_annotation = f"OmpK36:p.{insertion_start_aa}_{insertion_end_aa}ins{inserted_residues}"
                        
                        hits_dict['Omp_mutations'].append([
                            insertion_annotation,
                            {'Genetic_variation_type': 'Protein variant detected', **hit_data}
                        ])

    truncations = []
    if ompk35_frameshift:
        truncations.append(ompk35_frameshift)
    elif ompk35_deletion:
        truncations.append(ompk35_deletion)

    if ompk36_frameshift:
        truncations.append(ompk36_frameshift)
    elif ompk36_deletion:
        truncations.append(ompk36_deletion)

    for trunc_name, t_hit_data in truncations:
        data = dict(t_hit_data) if t_hit_data else {}
        hits_dict['Omp_mutations'].append([trunc_name, data])

    # --- NEW: Check for completely missing genes ---
    if ompk35_hit_data is None:
        hits_dict['Omp_mutations'].append([
            "OmpK35:del", 
            {"Genetic_variation_type": "Gene absent", "Coverage": "0.00%"}
        ])

    if ompk36_hit_data is None:
        hits_dict['Omp_mutations'].append([
            "OmpK36:del", 
            {"Genetic_variation_type": "Gene absent", "Coverage": "0.00%"}
        ])

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

#     aa_map = { '*': 'Ter'}

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

#         # --- Frameshift and Deletion checks in OmpK35---
#         if hit.query_name == 'OmpK35':
#             ompk35_ref_seq = ref_seqs['OmpK35']
#             ompk35_query_seq = hit.ref_seq
#             ompk35_hit_data = hit_data

#             # sometimes the hit does not start at the first base of the gene( gene is considered 0% cov)
#             # check for deletion at the start of the alignment 
#             if coverage == 0.0:
#                 aln = dna_aligner.align(ompk35_ref_seq, ompk35_query_seq)[0]
#                 # print(aln)
#                 deleted_base_pos = find_start_deletion_in_alignment(aln)
#                 deletion_report = f"OmpK35_{deleted_base_pos}"

#                 ompk35_deletion = (
#                     deletion_report,
#                     {**hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'}
#                 )
#                 continue  

#             if coverage > best_ompk35_cov:
#                 best_ompk35_cov = coverage
#                 ompk35_dna_cov = dna_hit_cov

#                 # check for frameshift mutation
#                 if best_ompk35_cov < 90.0 and ompk35_dna_cov > 90.0:
#                     if translation:
#                         ompk35_ref_trans = translate_nucl_to_prot(ompk35_ref_seq)
#                         ompk35_query_trans = translate_nucl_to_prot(ompk35_query_seq)
#                         ompk35_prot_align = aligner.align(ompk35_ref_trans, ompk35_query_trans)

#                         fs_info = get_frameshift_info(ompk35_prot_align[0])
#                         if fs_info is not None:
#                             aa_pos, ref_aa, alt_aa, fs_len = fs_info
#                             alt_str = aa_map.get(alt_aa, alt_aa)
#                             if fs_len == 0:
#                                 fs_report = f"OmpK35_{ref_aa}{aa_pos}{alt_str}"
#                             else:
#                                 fs_report = f"OmpK35_{ref_aa}{aa_pos}{alt_str}fsTer{fs_len}"

#                             ompk35_frameshift = (
#                                 fs_report,
#                                 {**hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'}
#                             )
                    
#                 # check for deletion
#                 elif best_ompk35_cov < 90.0 and ompk35_dna_cov < 90.0:
#                     ompk35_dna_alignments = dna_aligner.align(ompk35_ref_seq, ompk35_query_seq)
#                     del_info = deletion_checks(ompk35_dna_alignments[0], ompk35_ref_seq)
#                     if del_info is not None:
#                         pos, base = del_info
#                         del_report = f"OmpK35:c.{base}{pos}del"
#                         ompk35_deletion = (
#                             del_report,
#                             {**hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'}
#                         )

#         # --- Frameshift and Deletion checks in OmpK35 ---
#         elif hit.query_name == 'OmpK36':
#             ompk36_ref_seq = ref_seqs['OmpK36']
#             ompk36_query_seq = hit.ref_seq
#             ompk36_hit_data = hit_data
#             # check for deletion at the start of the alignment
#             if coverage == 0.0:
#                 aln = dna_aligner.align(ompk36_ref_seq, ompk36_query_seq)[0]
#                 deleted_base_pos = find_start_deletion_in_alignment(aln)
#                 deletion_report = f"OmpK36_{deleted_base_pos}"

#                 ompk36_deletion = (
#                     deletion_report,
#                     {**hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'}
#                 )
#                 continue 

#             if coverage > best_ompk36_cov:
#                 best_ompk36_cov = coverage
#                 ompk36_dna_cov = dna_hit_cov
#                 ompk36_dna_alignments = dna_aligner.align(ompk36_ref_seq, ompk36_query_seq)
#                 # check for frameshift mutation
#                 if best_ompk36_cov < 90.0 and ompk36_dna_cov > 90.0:
#                     if translation:
#                         ompk36_ref_trans = translate_nucl_to_prot(ompk36_ref_seq)
#                         ompk36_query_trans = translate_nucl_to_prot(ompk36_query_seq)
#                         ompk36_prot_align = aligner.align(ompk36_ref_trans, ompk36_query_trans)

#                         fs_info = get_frameshift_info(ompk36_prot_align[0])
#                         if fs_info is not None:
#                             aa_pos, ref_aa, alt_aa, fs_len = fs_info
#                             alt_str = aa_map.get(alt_aa, alt_aa)
#                             if fs_len == 0:
#                                 fs_report = f"OmpK36_{ref_aa}{aa_pos}{alt_str}"
#                             else:
#                                 fs_report = f"OmpK36_{ref_aa}{aa_pos}{alt_str}fsTer{fs_len}"

#                             ompk36_frameshift = (
#                                 fs_report,
#                                 {**hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'}
#                             )
                        
#                 # check for deletion
#                 elif best_ompk36_cov < 90.0 and ompk36_dna_cov < 90.0:
#                     del_info = deletion_checks(ompk36_dna_alignments[0], ompk36_ref_seq)
#                     if del_info is not None:
#                         pos, base = del_info
#                         del_report = f"OmpK36:c.{base}{pos}del"
#                         ompk36_deletion = (
#                             del_report,
#                             {**hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'}
#                         )

#                 # ----- OmpK 36 Nucleotide variant -----
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

#                 # L3 insertion checks (GD / TD / D)
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

#                         inserted_residues = insertion_type
#                         insertion_annotation = (
#                             f"OmpK36:p.{insertion_start_aa}_{insertion_end_aa}ins{inserted_residues}"
#                         )
#                         hits_dict['Omp_mutations'].append([
#                             insertion_annotation,
#                             {'Genetic_variation_type': 'Protein variant detected', **hit_data}
#                         ])

#         else:
#             continue

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



# def get_frameshift_info(alignment):
#     """
#     Detects frameshift mutations in a protein alignment and reports
#     premature stop codons.
#     """
#     reference_protein = alignment[0]
#     query_protein = alignment[1]

#     ref_aa_position = 0
#     seq_len = len(reference_protein)

#     SCAN_WINDOW = 20
#     ANCHOR_LENGTH = 3

#     for i in range(seq_len):
#         ref_aa = reference_protein[i]
#         query_aa = query_protein[i]

#         if ref_aa != '-':
#             ref_aa_position += 1

#         # Detect (mismatch or indel)
#         if ref_aa != query_aa:

#             is_frameshift = True

#             # stop codon in query
#             if query_aa == '*':
#                 is_frameshift = True
#             else:
#                 search_end = min(seq_len, i + SCAN_WINDOW)
#                 for k in range(i + 1, search_end):
#                     if k + ANCHOR_LENGTH > seq_len:
#                         break

#                     ref_fragment = reference_protein[k:k + ANCHOR_LENGTH]
#                     query_fragment = query_protein[k:k + ANCHOR_LENGTH]

#                     if ref_fragment == query_fragment and '-' not in ref_fragment:
#                         is_frameshift = False
#                         break

#             if not is_frameshift:
#                 continue

#             # Determine true frameshift position and amino acid
#             frameshift_ref_pos = ref_aa_position
#             frameshift_ref_aa = ref_aa

#             ref_remaining = reference_protein[i:].replace('-', '').replace('.', '')
#             query_remaining = query_protein[i:].replace('-', '').replace('.', '')

#             if ref_aa == '-':
#                 frameshift_ref_aa = ref_remaining[0] if ref_remaining else '-'
#                 frameshift_ref_pos = ref_aa_position + 1

#             frameshift_query_aa = query_remaining[0] if query_remaining else '*'

#             ref_idx = 0
#             query_idx = 0

#             while (
#                 frameshift_ref_aa == frameshift_query_aa
#                 and frameshift_ref_aa not in ('-', '*')
#             ):
#                 ref_idx += 1
#                 if ref_idx < len(ref_remaining):
#                     frameshift_ref_aa = ref_remaining[ref_idx]
#                     frameshift_ref_pos += 1
#                 else:
#                     frameshift_ref_aa = '-'

#                 query_idx += 1
#                 if query_idx < len(query_remaining):
#                     frameshift_query_aa = query_remaining[query_idx]
#                 else:
#                     frameshift_query_aa = '*'

#             # Determine length of frameshift until stop codon
#             query_after_frameshift = query_remaining[query_idx:]

#             if '*' in query_after_frameshift:
#                 frameshift_length = query_after_frameshift.index('*')
#             else:
#                 frameshift_length = len(query_after_frameshift)

#             return (
#                 frameshift_ref_pos,
#                 frameshift_ref_aa,
#                 frameshift_query_aa,
#                 frameshift_length
#             )

#     return None





# def deletion_checks(alignment, ref_seq):
#     """
#     Identifies deletion.
#     Returns: (position, deleted_base).
#     """
#     aligned_ref, aligned_query = alignment[0], alignment[1]

#     deletions = [] 
#     ref_pos = 0
#     deletion_len = 0
#     deletion_start = None

#     # Scan alignment for deletions
#     for ref_aln, query_aln in zip(aligned_ref, aligned_query):
#         if ref_aln != '-':
#             ref_pos += 1

#         if query_aln == '-' and ref_aln != '-':
#             if deletion_len == 0:
#                 deletion_start = ref_pos
#             deletion_len += 1
#         else:
#             if deletion_len > 0:
#                 deletions.append((deletion_start, deletion_len))
#                 deletion_len = 0
#                 deletion_start = None

    
#     if deletion_len > 0:
#         deletions.append((deletion_start, deletion_len))

#     if ref_pos < len(ref_seq):
#         trunc_start = ref_pos + 1
#         trunc_len = len(ref_seq) - ref_pos
#         deletions.append((trunc_start, trunc_len))

#     if not deletions:
#         return None

#     best_start, best_len = max(deletions, key=lambda x: x[1])

#     deleted_base = ref_seq[best_start - 1]

#     return best_start, deleted_base




# def find_start_deletion_in_alignment(alignment):
#     """
#     Finds deletion if the hit does not start at the first base of the gene
#     """
#     # retrieves the alignment coordinates for the Reference sequence
#     ref_alignment= alignment.aligned[0]

#     # Get the start index of the first aligned block (0-based index)
#     match_index = ref_alignment[0][0]

#     # Convert 0-based Python index to 1-based
#     start_del = 1
#     end_del = match_index

#     if start_del == end_del:
#         return f"c.({start_del})del"
#     else:
#         return f"c.({start_del}_{end_del})del"

