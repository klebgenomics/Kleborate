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

    # Amino acid single-letter to three-letter code mapping
    aa_map = {
        'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe', 'G': 'Gly',
        'H': 'His', 'I': 'Ile', 'K': 'Lys', 'L': 'Leu', 'M': 'Met', 'N': 'Asn',
        'P': 'Pro', 'Q': 'Gln', 'R': 'Arg', 'S': 'Ser', 'T': 'Thr', 'V': 'Val',
        'W': 'Trp', 'Y': 'Tyr'
    }

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
                'Input_sequence_ID': hit.ref_name,
                'Input_gene_length': hit.ref_length,
                'Input_gene_start': hit.ref_start,
                'Input_gene_stop': hit.ref_end,
                'Reference_gene_length': hit.query_length,
                'Reference_gene_start': hit.query_start,
                'Reference_gene_stop': hit.query_end,
                'Sequence_identity': f"{hit.percent_identity:.2f}%",
                'Coverage': f"{hit.query_cov:.2f}%", # Dna hit coverage
                'Strand_orientation':hit.strand
            }


        if coverage > min_coverage:
            if hit.query_name == 'GyrA':
                alignments = protein_aligner.align(gyra_ref, translation)
            elif hit.query_name == 'ParC':
                alignments = protein_aligner.align(parc_ref, translation)
            else:
                assert False

            bases_per_ref_pos = get_bases_per_ref_pos(alignments[0])
            loci = qrdr_loci[hit.query_name]

            for pos, wt_base in loci:
                assembly_base = bases_per_ref_pos[pos]

                if pos in bases_per_ref_pos and assembly_base != wt_base \
                        and assembly_base != '-' and assembly_base != '.':
                    # Directly map the amino acids to their three-letter codes
                    mutation = f"{hit.query_name}:p.{aa_map.get(wt_base)}{pos}{aa_map.get(assembly_base)}"
                    snps.append([mutation, {'Genetic_variation_type': 'Protein variant detected'},hit_data])

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
