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
from ...shared.alignment import align_query_to_ref, truncation_check, translate_nucl_to_prot, get_last_base_aa_before_gap


def check_for_mgrb_pmrb_gene_truncations(hits_dict, assembly, trunc, min_ident):

    best_mgrb_cov, best_pmrb_cov = 0.0, 0.0
    mgrb_hit_data, pmrb_hit_data = None, None
    mgrb_frameshift, mgrb_deletion = None, None
    pmrb_frameshift, pmrb_deletion = None, None
    mgrb_hit, pmrb_hit = None, None  

    start_codons = {'TTG', 'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTG'}

    aa_map = {
        'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe', 'G': 'Gly',
        'H': 'His', 'I': 'Ile', 'K': 'Lys', 'L': 'Leu', 'M': 'Met', 'N': 'Asn',
        'P': 'Pro', 'Q': 'Gln', 'R': 'Arg', 'S': 'Ser', 'T': 'Thr', 'V': 'Val',
        'W': 'Trp', 'Y': 'Tyr'
    }

    # define the aligner
    aligner = Align.PairwiseAligner(mode='global', match_score=5, mismatch_score=-4,
                                    open_gap_score=-10, extend_gap_score=-0.5)

    alignment_hits = align_query_to_ref(trunc, assembly, None, min_identity=None)

    if 'Col_mutations' not in hits_dict:
        hits_dict['Col_mutations'] = []

    for hit in alignment_hits:
        assert hit.query_name == 'pmrB' or hit.query_name == 'mgrB'
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

        if hit.query_name == 'mgrB':
            if coverage > best_mgrb_cov:
                best_mgrb_cov = coverage
                mgrb_hit_data = hit_data
                mgrb_dna_cov = dna_hit_cov
                mgrb_hit = hit.ref_seq
                mgrb_ref_seq = hit.query_seq

                # ---- FRAMESHIFT/DELETION LOGIC for mgrB ----
                if best_mgrb_cov < 99.0 and mgrb_dna_cov > 99.0:
                    # checks for Frameshift
                    if translation:
                        mgrb_ref_trans = translate_nucl_to_prot(mgrb_ref_seq)
                        mgrb_prot_align= aligner.align(mgrb_ref_trans, translation)
                        aa_pos, aa_code = get_last_base_aa_before_gap(mgrb_prot_align[0])
                        fs_report = f"mgrb:p.{aa_map[aa_code]}{aa_pos}fs"
                        mgrb_frameshift = (fs_report, {**hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'})

                elif best_mgrb_cov < 99.0 and mgrb_dna_cov < 99.0:
                    # checks for Deletion
                    mgrb_dna_alignments = aligner.align(mgrb_ref_seq, mgrb_hit)
                    pos, base = get_last_base_aa_before_gap(mgrb_dna_alignments[0])
                    del_report = f"mgrb:c.{base}{pos}del"
                    mgrb_deletion = (del_report, {**hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'})

        elif hit.query_name == 'pmrB':
            if coverage > best_pmrb_cov:
                best_pmrb_cov = coverage
                pmrb_hit_data = hit_data
                pmrb_dna_cov = dna_hit_cov 
                pmrb_hit = hit.ref_seq 
                pmrb_ref_seq = hit.query_seq

                # ---- FRAMESHIFT/DELETION LOGIC for pmrB ----
                if best_pmrb_cov < 99.0 and pmrb_dna_cov > 99.0:
                    if translation:
                        pmrb_ref_trans = translate_nucl_to_prot(pmrb_ref_seq)
                        pmrb_prot_align= aligner.align(pmrb_ref_trans, translation)
                        aa_pos, aa_code = get_last_base_aa_before_gap(pmrb_prot_align[0])

                        fs_report = f"pmrb:p.{aa_map[aa_code]}{aa_pos}fs"
                        pmrb_frameshift = (fs_report, {**hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'})

                elif best_pmrb_cov < 99.0 and pmrb_dna_cov < 99.0:
                    # DELETION
                    pmrb_dna_alignments = aligner.align(pmrb_ref_seq, pmrb_hit)
                    pos, base = get_last_base_aa_before_gap(pmrb_dna_alignments[0])
                    del_report = f"pmrb:c.{base}{pos}del"
                    pmrb_deletion = (del_report, {**hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'})
        else:
            assert False

    # Summarise detected truncations and add to results
    truncations = []
    if mgrb_frameshift:
        truncations.append(mgrb_frameshift)
    elif mgrb_deletion:
        truncations.append(mgrb_deletion)
    else:
        # Check if the first 3 bp of the hit is a known start codon for mgrB
        if mgrb_hit:
            mgrb_start_codon = mgrb_hit[:3].upper()
            if mgrb_start_codon not in start_codons:
                truncations.append(('mgrB$', {**mgrb_hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'}))

    if pmrb_frameshift:
        truncations.append(pmrb_frameshift)
    elif pmrb_deletion:
        truncations.append(pmrb_deletion)
    else:
        # Check if the first 3 bp of the hit is a known start codon for pmrB
        if pmrb_hit:
            pmrb_start_codon = pmrb_hit[:3].upper()
            if pmrb_start_codon not in start_codons:
                truncations.append(('pmrB$', {**pmrb_hit_data, 'Genetic_variation_type': 'Inactivating mutation detected'}))

    for trunc_name, hit_data in truncations:
        data = dict(hit_data) if hit_data else {}
        hits_dict['Col_mutations'].append([trunc_name, data])



# def check_for_mgrb_pmrb_gene_truncations(hits_dict, assembly, trunc, min_ident):
#     best_mgrb_cov, best_pmrb_cov = 0.0, 0.0
#     mgrb_hit, pmrb_hit = None, None
#     start_codons = {'TTG', 'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTG'}

#     alignment_hits = align_query_to_ref(trunc, assembly, None, min_identity=None)
#     for hit in alignment_hits:
#         assert hit.query_name == 'pmrB' or hit.query_name == 'mgrB'
#         _, coverage, _ = truncation_check(hit)
        
#         if hit.query_name == 'mgrB' and coverage > best_mgrb_cov:
#             best_mgrb_cov = coverage
#             mgrb_hit = hit.ref_seq
#         elif hit.query_name == 'pmrB' and coverage > best_pmrb_cov:
#             best_pmrb_cov = coverage
#             pmrb_hit = hit.ref_seq
            

#     truncations = []
#     if best_mgrb_cov < 90.0:
#         truncations.append('MgrB-' + ('%.0f' % best_mgrb_cov) + '%')
#     else:
#         # Check if the first 3 bp of the hit is a known start codon for mgrB
#         if mgrb_hit:
#             hit_start_codon = mgrb_hit[:3].upper()
#             if hit_start_codon not in start_codons:
#                 truncations.append('mgrB$')

#     if best_pmrb_cov < 90.0:
#         truncations.append('PmrB-' + ('%.0f' % best_pmrb_cov) + '%')
#     else:
#         # Check if the first 3 bp of the hit is a known start codon for pmrB
#         if pmrb_hit:
#             hit_start_codon = pmrb_hit[:3].upper()
#             if hit_start_codon not in start_codons:
#                 truncations.append('pmrB$')

#     if truncations:
#         hits_dict['Col_mutations'] += truncations
