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
from ...shared.alignment import align_query_to_ref, truncation_check, translate_nucl_to_prot, find_start_deletion_in_alignment, deletion_checks, get_frameshift_info
from ...shared.misc import load_fasta, reverse_complement


def check_for_mgrb_pmrb_gene_truncations(hits_dict, assembly, trunc, min_ident):

    best_mgrb_cov, best_pmrb_cov = 0.0, 0.0
    mgrb_hit_data, pmrb_hit_data = None, None
    mgrb_frameshift, mgrb_deletion = None, None
    pmrb_frameshift, pmrb_deletion = None, None
    mgrb_query_seq, pmrb_query_seq = None, None  

    start_codons = {'TTG', 'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTG'}
    ref_seqs = dict(load_fasta(trunc))

    aa_map = { '*': 'Ter'}

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

    alignment_hits = align_query_to_ref(trunc, assembly, None, min_identity=None)

    if 'Col_mutations' not in hits_dict:
        hits_dict['Col_mutations'] = []

    for hit in alignment_hits:
        assert hit.query_name == 'pmrB' or hit.query_name == 'mgrB'
        _, coverage, translation = truncation_check(hit)
        dna_hit_cov = hit.query_cov

        hit_data = {
            'Input Sequence ID': hit.ref_name,
            'Input Gene Length': hit.ref_length,
            'Input Gene Start': hit.ref_start +1,
            'Input Gene Stop': hit.ref_end,
            'Reference Gene Length': hit.query_length,
            'Reference Gene Start': hit.query_start + 1,
            'Reference Gene Stop': hit.query_end,
            'Sequence Identity': f"{hit.percent_identity:.2f}",
            'Coverage': f"{coverage:.2f}",
            'Strand Orientation': hit.strand
        }

        # --- mgrB checks ---
        if hit.query_name == 'mgrB':
            if coverage >= best_mgrb_cov:
                best_mgrb_cov = coverage
                mgrb_dna_cov = dna_hit_cov
                mgrb_hit_data = hit_data
                mgrb_query_seq = hit.ref_seq
                mgrb_ref_seq = ref_seqs['mgrB']

                mgrb_frameshift, mgrb_deletion = None, None

                if coverage == 0.0:
                    aln = dna_aligner.align(mgrb_ref_seq, mgrb_query_seq)[0]
                    deleted_base_pos = find_start_deletion_in_alignment(aln)
                    mgrb_deletion = (f"mgrB:{deleted_base_pos}", {**hit_data, 'Genetic Variation Type': 'Inactivating mutation detected'})
                
                elif best_mgrb_cov < 90.0 and mgrb_dna_cov > 90.0:
                    if translation:
                        mgrb_ref_trans = translate_nucl_to_prot(mgrb_ref_seq)
                        mgrb_query_trans = translate_nucl_to_prot(mgrb_query_seq)
                        mgrb_prot_align = aligner.align(mgrb_ref_trans, mgrb_query_trans)
                        fs_info = get_frameshift_info(mgrb_prot_align[0])
                        if fs_info is not None:
                            aa_pos, ref_aa, alt_aa, fs_len = fs_info
                            alt_str = aa_map.get(alt_aa, alt_aa)
                            fs_report = f"mgrB:p.{ref_aa}{aa_pos}{alt_str}" if fs_len == 0 else f"mgrB:p.{ref_aa}{aa_pos}{alt_str}fsTer{fs_len}"
                            mgrb_frameshift = (fs_report, {**hit_data, 'Genetic Variation Type': 'Inactivating mutation detected'})
                
                elif best_mgrb_cov < 90.0 and mgrb_dna_cov < 90.0:
                    mgrb_dna_alignments = dna_aligner.align(mgrb_ref_seq, mgrb_query_seq)
                    del_info = deletion_checks(mgrb_dna_alignments[0], mgrb_ref_seq)
                    if del_info is not None:
                        pos, base = del_info
                        mgrb_deletion = (f"mgrB:c.{base}{pos}del", {**hit_data, 'Genetic Variation Type': 'Inactivating mutation detected'})

        # --- pmrB checks ---
        elif hit.query_name == 'pmrB':
            if coverage >= best_pmrb_cov:
                best_pmrb_cov = coverage
                pmrb_dna_cov = dna_hit_cov 
                pmrb_hit_data = hit_data
                pmrb_query_seq = hit.ref_seq
                pmrb_ref_seq = ref_seqs['pmrB']

                pmrb_frameshift, pmrb_deletion = None, None

                if coverage == 0.0:
                    aln = dna_aligner.align(pmrb_ref_seq, pmrb_query_seq)[0]
                    deleted_base_pos = find_start_deletion_in_alignment(aln)
                    pmrb_deletion = (f"pmrB:{deleted_base_pos}", {**hit_data, 'Genetic Variation Type': 'Inactivating mutation detected'})
                
                elif best_pmrb_cov < 90.0 and pmrb_dna_cov > 90.0:
                    if translation:
                        pmrb_ref_trans = translate_nucl_to_prot(pmrb_ref_seq)
                        pmrb_query_trans = translate_nucl_to_prot(pmrb_query_seq)
                        pmrb_prot_align = aligner.align(pmrb_ref_trans, pmrb_query_trans)
                        fs_info = get_frameshift_info(pmrb_prot_align[0])
                        if fs_info is not None:
                            aa_pos, ref_aa, alt_aa, fs_len = fs_info
                            alt_str = aa_map.get(alt_aa, alt_aa)
                            fs_report = f"pmrB:p.{ref_aa}{aa_pos}{alt_str}" if fs_len == 0 else f"pmrB:p.{ref_aa}{aa_pos}{alt_str}fsTer{fs_len}"
                            pmrb_frameshift = (fs_report, {**hit_data, 'Genetic Variation Type': 'Inactivating mutation detected'})   

                elif best_pmrb_cov < 90.0 and pmrb_dna_cov < 90.0:
                    pmrb_dna_alignments = dna_aligner.align(pmrb_ref_seq, pmrb_query_seq)
                    del_info = deletion_checks(pmrb_dna_alignments[0], pmrb_ref_seq)
                    if del_info is not None:
                        pos, base = del_info
                        pmrb_deletion = (f"pmrB:c.{base}{pos}del", {**hit_data, 'Genetic Variation Type': 'Inactivating mutation detected'})

    # --- Truncation reporting logic ---
    truncations = []
    
    if mgrb_frameshift:
        truncations.append(mgrb_frameshift)
    elif mgrb_deletion:
        truncations.append(mgrb_deletion)
    elif mgrb_hit_data is None:
        truncations.append(("mgrB:del", {"Genetic Variation Type": "Gene deletion detected", "Coverage": "0.00"}))
    else:
        if mgrb_query_seq:
            mgrb_start_codon = mgrb_query_seq[:3].upper()
            if mgrb_start_codon not in start_codons:
                truncations.append(('mgrB:p.(Met1?)', {**mgrb_hit_data, 'Genetic Variation Type': 'Inactivating mutation detected'}))


    if pmrb_frameshift:
        truncations.append(pmrb_frameshift)
    elif pmrb_deletion:
        truncations.append(pmrb_deletion)
    elif pmrb_hit_data is None:
        truncations.append(("pmrB:del", {"Genetic Variation Type": "Gene deletion detected", "Coverage": "0.00"}))
    else:
        if pmrb_query_seq:
            pmrb_start_codon = pmrb_query_seq[:3].upper()
            if pmrb_start_codon not in start_codons:
                truncations.append(('pmrB:p.(Met1?)', {**pmrb_hit_data, 'Genetic Variation Type': 'Inactivating mutation detected'}))

    for trunc_name, hit_metadata in truncations:
        data = dict(hit_metadata) if hit_metadata else {}
        hits_dict['Col_mutations'].append([trunc_name, data])





# def check_for_mgrb_pmrb_gene_truncations(hits_dict, assembly, trunc, min_ident):

#     best_mgrb_cov, best_pmrb_cov = 0.0, 0.0
#     mgrb_hit_data, pmrb_hit_data = None, None
#     mgrb_frameshift, mgrb_deletion = None, None
#     pmrb_frameshift, pmrb_deletion = None, None
#     mgrb_query_seq, pmrb_query_seq = None, None  

#     start_codons = {'TTG', 'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTG'}
#     ref_seqs = dict(load_fasta(trunc))

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

#     alignment_hits = align_query_to_ref(trunc, assembly, None, min_identity=None)
#     # print(alignment_hits)

#     if 'Col_mutations' not in hits_dict:
#         hits_dict['Col_mutations'] = []

#     for hit in alignment_hits:
#         # print(hit.query_name, hit.ref_seq)
#         assert hit.query_name == 'pmrB' or hit.query_name == 'mgrB'
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

#         # --- mgrB checks ---
#         if hit.query_name == 'mgrB':
#             if coverage >= best_mgrb_cov:
#                 best_mgrb_cov = coverage
#                 mgrb_dna_cov = dna_hit_cov
#                 mgrb_hit_data = hit_data
#                 mgrb_query_seq = hit.ref_seq
#                 mgrb_ref_seq = ref_seqs['mgrB']

#                 mgrb_frameshift, mgrb_deletion = None, None

#                 if coverage == 0.0:
#                     aln = dna_aligner.align(mgrb_ref_seq, mgrb_query_seq)[0]
#                     deleted_base_pos = find_start_deletion_in_alignment(aln)
#                     mgrb_deletion = (f"mgrB:{deleted_base_pos}", {**hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'})
                
#                 elif best_mgrb_cov < 90.0 and mgrb_dna_cov > 90.0:
#                     if translation:
#                         mgrb_ref_trans = translate_nucl_to_prot(mgrb_ref_seq)
#                         mgrb_query_trans = translate_nucl_to_prot(mgrb_query_seq)
#                         mgrb_prot_align = aligner.align(mgrb_ref_trans, mgrb_query_trans)
#                         fs_info = get_frameshift_info(mgrb_prot_align[0])
#                         if fs_info is not None:
#                             aa_pos, ref_aa, alt_aa, fs_len = fs_info
#                             alt_str = aa_map.get(alt_aa, alt_aa)
#                             fs_report = f"mgrB:p.{ref_aa}{aa_pos}{alt_str}" if fs_len == 0 else f"mgrB:p.{ref_aa}{aa_pos}{alt_str}fsTer{fs_len}"
#                             mgrb_frameshift = (fs_report, {**hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'})
                
#                 elif best_mgrb_cov < 90.0 and mgrb_dna_cov < 90.0:
#                     mgrb_dna_alignments = dna_aligner.align(mgrb_ref_seq, mgrb_query_seq)
#                     del_info = deletion_checks(mgrb_dna_alignments[0], mgrb_ref_seq)
#                     if del_info is not None:
#                         pos, base = del_info
#                         mgrb_deletion = (f"mgrB:c.{base}{pos}del", {**hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'})

#         # --- pmrB checks ---
#         elif hit.query_name == 'pmrB':
#             if coverage >= best_pmrb_cov:
#                 best_pmrb_cov = coverage
#                 pmrb_dna_cov = dna_hit_cov 
#                 pmrb_hit_data = hit_data
#                 pmrb_query_seq = hit.ref_seq
#                 pmrb_ref_seq = ref_seqs['pmrB']

#                 pmrb_frameshift, pmrb_deletion = None, None

#                 if coverage == 0.0:
#                     aln = dna_aligner.align(pmrb_ref_seq, pmrb_query_seq)[0]
#                     deleted_base_pos = find_start_deletion_in_alignment(aln)
#                     pmrb_deletion = (f"pmrB:{deleted_base_pos}", {**hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'})
                
#                 elif best_pmrb_cov < 90.0 and pmrb_dna_cov > 90.0:
#                     if translation:
#                         pmrb_ref_trans = translate_nucl_to_prot(pmrb_ref_seq)
#                         pmrb_query_trans = translate_nucl_to_prot(pmrb_query_seq)
#                         pmrb_prot_align = aligner.align(pmrb_ref_trans, pmrb_query_trans)
#                         fs_info = get_frameshift_info(pmrb_prot_align[0])
#                         if fs_info is not None:
#                             aa_pos, ref_aa, alt_aa, fs_len = fs_info
#                             alt_str = aa_map.get(alt_aa, alt_aa)
#                             fs_report = f"pmrB:p.{ref_aa}{aa_pos}{alt_str}" if fs_len == 0 else f"pmrB:p.{ref_aa}{aa_pos}{alt_str}fsTer{fs_len}"
#                             pmrb_frameshift = (fs_report, {**hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'})   

#                 elif best_pmrb_cov < 90.0 and pmrb_dna_cov < 90.0:
#                     pmrb_dna_alignments = dna_aligner.align(pmrb_ref_seq, pmrb_query_seq)
#                     del_info = deletion_checks(pmrb_dna_alignments[0], pmrb_ref_seq)
#                     if del_info is not None:
#                         pos, base = del_info
#                         pmrb_deletion = (f"pmrB:c.{base}{pos}del", {**hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'})

#     # --- Truncation reporting logic ---
#     truncations = []
    
#     if mgrb_frameshift:
#         truncations.append(mgrb_frameshift)
#     elif mgrb_deletion:
#         truncations.append(mgrb_deletion)
#     elif mgrb_hit_data is None:
#         truncations.append(("mgrB:del", {"Genetic_variation_type": "Gene deleted", "Coverage": "0.00%"}))
#     else:
#         # Check start codon only if no truncation found and gene exists
#         if mgrb_query_seq:
#             mgrb_start_codon = mgrb_query_seq[:3].upper()
#             if mgrb_start_codon not in start_codons:
#                 truncations.append(('mgrB$', {**mgrb_hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'}))


#     if pmrb_frameshift:
#         truncations.append(pmrb_frameshift)
#     elif pmrb_deletion:
#         truncations.append(pmrb_deletion)
#     elif pmrb_hit_data is None:
#         truncations.append(("pmrB:del", {"Genetic_variation_type": "Gene deleted", "Coverage": "0.00%"}))
#     else:
#         if pmrb_query_seq:
#             pmrb_start_codon = pmrb_query_seq[:3].upper()
#             if pmrb_start_codon not in start_codons:
#                 truncations.append(('pmrB$', {**pmrb_hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'}))

#     for trunc_name, hit_metadata in truncations:
#         data = dict(hit_metadata) if hit_metadata else {}
#         hits_dict['Col_mutations'].append([trunc_name, data])