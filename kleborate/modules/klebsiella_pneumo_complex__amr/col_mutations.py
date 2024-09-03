"""
Copyright 2024 Kat Holt
Copyright 2024 Ryan Wick (rrwick@gmail.com)
Copyright 2024 (gathonimaranga@gmail.com)
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
from Bio.Align import substitution_matrices
from ...shared.alignment import align_query_to_ref, truncation_check


def check_for_mgrb_pmrb_gene_truncations(hits_dict, assembly, trunc, min_ident):
    best_mgrb_cov, best_pmrb_cov = 0.0, 0.0
    mgrb_hit, pmrb_hit = None, None
    start_codons = {'TTG', 'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTG'}

    alignment_hits = align_query_to_ref(trunc, assembly, None, min_identity=None)
    for hit in alignment_hits:
        assert hit.query_name == 'pmrB' or hit.query_name == 'mgrB'
        _, coverage, _ = truncation_check(hit)
        
        if hit.query_name == 'mgrB' and coverage > best_mgrb_cov:
            best_mgrb_cov = coverage
            mgrb_hit = hit.ref_seq
        elif hit.query_name == 'pmrB' and coverage > best_pmrb_cov:
            best_pmrb_cov = coverage
            pmrb_hit = hit.ref_seq
            

    truncations = []
    if best_mgrb_cov < 90.0:
        truncations.append('MgrB-' + ('%.0f' % best_mgrb_cov) + '%')
    else:
        # Check if the first 3 bp of the hit is a known start codon for mgrB
        if mgrb_hit:
            hit_start_codon = mgrb_hit[:3].upper()
            if hit_start_codon not in start_codons:
                truncations.append('mgrB$')

    if best_pmrb_cov < 90.0:
        truncations.append('PmrB-' + ('%.0f' % best_pmrb_cov) + '%')
    else:
        # Check if the first 3 bp of the hit is a known start codon for pmrB
        if pmrb_hit:
            hit_start_codon = pmrb_hit[:3].upper()
            if hit_start_codon not in start_codons:
                truncations.append('pmrB$')

    if truncations:
        hits_dict['Col_mutations'] += truncations



# def check_for_mgrb_pmrb_gene_truncations(hits_dict, assembly, trunc,min_ident):
#     best_mgrb_cov, best_pmrb_cov = 0.0, 0.0

#     alignment_hits = align_query_to_ref(trunc, assembly, None, min_identity=None)
#     for hit in alignment_hits:
#         assert hit.query_name == 'pmrB' or hit.query_name == 'mgrB'
#         _, coverage, _ = truncation_check(hit)
        
#         if hit.query_name == 'mgrB' and coverage > best_mgrb_cov:
#             best_mgrb_cov = coverage
#         elif hit.query_name == 'pmrB' and coverage > best_pmrb_cov:
#             best_pmrb_cov = coverage
            

#     truncations = []
#     if best_mgrb_cov < 90.0:
#         truncations.append('MgrB-' + ('%.0f' % best_mgrb_cov) + '%')
#     if best_pmrb_cov < 90.0:
#         truncations.append('PmrB-' + ('%.0f' % best_pmrb_cov) + '%')

#     if truncations:
#         hits_dict['Col_mutations'] += truncations
