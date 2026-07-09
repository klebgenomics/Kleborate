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
from ...shared.alignment import align_query_to_ref, truncation_check, get_bases_per_ref_pos

def check_for_qrdr_mutations(hits_dict, assembly, qrdr, min_identity, min_coverage):
    """
    This function checks for qrdr mutations

    This function returns:
    * a hits dictionary with Fluoroquinolone(Qrdr) mutations
    """

    qrdr_loci = {'GyrA': [(83, 'S'), (87, 'D')],
                 'ParC': [(80, 'S'), (84, 'E')]}

    gyra_ref = 'MSDLAREITPVNIEEELKNSYLDYAMSVIVGRALPDVRDGLKPVHRRVLYAMNVLGNDWN' \
               'KAYKKSARVVGDVIGKYHPHGDSAVYDTIVRMAQPFSLRYMLVDGQGNFGSIDGDSAAAM'

    parc_ref = 'MSDMAERLALHEFTENAYLNYSMYVIMDRALPFIGDGLKPVQRRIVYAMSELGLNASAKF' \
               'KKSARTVGDVLGKYHPHGDSACYEAMVLMAQPFSYRYPLVDGQGNWGAPDDPKSFAAMRY'

    # Define PairwiseAligner
    protein_aligner = Align.PairwiseAligner()
    protein_aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    protein_aligner.open_gap_score = -10
    protein_aligner.extend_gap_score = -0.5

    snps = []

    alignment_hits = align_query_to_ref(qrdr, assembly, min_query_coverage=None, min_identity=min_identity)
    
    for hit in alignment_hits:
        _, coverage, translation = truncation_check(hit)

        hit_data = {
                'Input Sequence ID': hit.ref_name,
                'Input Gene Length': hit.ref_length,
                'Input Gene Start': hit.ref_start,
                'Input Gene Stop': hit.ref_end,
                'Reference Gene Length': hit.query_length,
                'Reference Gene Start': hit.query_start +1,
                'Reference Gene Stop': hit.query_end,
                'Sequence Identity': f"{hit.percent_identity:.2f}",
                'Coverage': f"{hit.query_cov:.2f}", 
                'Strand Orientation':hit.strand
            }


        if coverage > min_coverage:
            if hit.query_name == 'GyrA':
                alignments = protein_aligner.align(gyra_ref, translation)
            elif hit.query_name == 'ParC':
                alignments = protein_aligner.align(parc_ref, translation)
            else:
                assert False

            coords = alignments[0].coordinates
            ref_protein_start  = int(coords[0, 0]) + 1
            ref_protein_stop   = int(coords[0, -1])
            ref_protein_length = ref_protein_stop - ref_protein_start + 1

            input_protein_start  = int(coords[1, 0]) + 1
            input_protein_stop   = int(coords[1, -1])
            input_protein_length = int(coords[1, -1]) - input_protein_start + 1

            hit_data.update({
                'Input Protein Length':     input_protein_length,
                'Input Protein Start':      input_protein_start,
                'Input Protein Stop':       input_protein_stop,
                'Reference Protein Length': ref_protein_length,
                'Reference Protein Start':  ref_protein_start,
                'Reference Protein Stop':   ref_protein_stop,
            })

            bases_per_ref_pos = get_bases_per_ref_pos(alignments[0])
            loci = qrdr_loci[hit.query_name]

            for pos, wt_base in loci:
                assembly_base = bases_per_ref_pos[pos]

                if pos in bases_per_ref_pos and assembly_base != wt_base \
                        and assembly_base != '-' and assembly_base != '.':
                    
                    # mutation = f"{hit.query_name}:p.{wt_base}{pos}{assembly_base}"
                    mutation = f"{hit.query_name[0].lower() + hit.query_name[1:]}:p.{wt_base}{pos}{assembly_base}"
                    snps.append([mutation, {'Genetic Variation Type': 'Protein variant detected'},hit_data])

    if snps:
        hits_dict['Flq_mutations'].extend(snps)



# def check_for_qrdr_mutations(hits_dict, assembly, qrdr, min_identity, min_coverage):
    
#     """
#     This function checks for qrdr mutations
    
#     This function returns:
#     * a hits dictionary with Fluoroquinolone(Qrdr) mutations
#     """

#     qrdr_loci = {'GyrA': [(83, 'S'), (87, 'D')],
#                      'ParC': [(80, 'S'), (84, 'E')]}

#     gyra_ref = 'MSDLAREITPVNIEEELKNSYLDYAMSVIVGRALPDVRDGLKPVHRRVLYAMNVLGNDWN' \
#                'KAYKKSARVVGDVIGKYHPHGDSAVYDTIVRMAQPFSLRYMLVDGQGNFGSIDGDSAAAM'
#     parc_ref = 'MSDMAERLALHEFTENAYLNYSMYVIMDRALPFIGDGLKPVQRRIVYAMSELGLNASAKF' \
#                'KKSARTVGDVLGKYHPHGDSACYEAMVLMAQPFSYRYPLVDGQGNWGAPDDPKSFAAMRY'

#     # define PairwiseAligner
#     protein_aligner = Align.PairwiseAligner()
#     protein_aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
#     protein_aligner.open_gap_score = -10
#     protein_aligner.extend_gap_score = -0.5

#     snps = []

#     alignment_hits = align_query_to_ref(qrdr, assembly, min_query_coverage=None, min_identity=min_identity) 
#     for hit in alignment_hits:
#         _, coverage, translation = truncation_check(hit)
        
#         if coverage > min_coverage:
#             if hit.query_name == 'GyrA':
#                 alignments = protein_aligner.align(gyra_ref, translation)
#             elif hit.query_name == 'ParC':
#                 alignments = protein_aligner.align(parc_ref, translation)
#             else:
#                 assert False
#             bases_per_ref_pos = get_bases_per_ref_pos(alignments[0])
#             loci = qrdr_loci[hit.query_name]

#             for pos, wt_base in loci:
#                 assembly_base = bases_per_ref_pos[pos]
#                 if pos in bases_per_ref_pos and assembly_base != wt_base \
#                         and assembly_base != '-' and assembly_base != '.':
#                     snps.append(hit.query_name + '-' + str(pos) + assembly_base)
        
#     if snps:
#         hits_dict['Flq_mutations'] += snps
