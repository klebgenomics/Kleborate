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
from Bio import Align
from Bio.Align import substitution_matrices
from ...shared.misc import reverse_complement

def check_for_shv_mutations(hit, hit_allele, bla_class, exact_match):
    
    # Don't do anything on non-SHV genes.
    if 'SHV' not in hit_allele:
        return bla_class, [], [], None

    nucl_seq = hit.ref_seq

    # If there are any ambiguous bases in the sequence, then we can't do the translation.
    ambiguous_bases = set(b for b in nucl_seq) - {'A', 'C', 'G', 'T'}
    if ambiguous_bases:
        return bla_class, [], [], None

    # BioPython doesn't like it if the sequence isn't a multiple of 3.
    nucl_seq = nucl_seq[:len(nucl_seq) // 3 * 3]

    coding_dna = Seq(nucl_seq)
    translation = str(coding_dna.translate(table='Bacterial', to_stop=True))

    shv_1_ref = 'MRYIRLCIISLLATLPLAVHASPQPLEQIKLSESQLSGRVGMIEMDLASGRTLTAWRADERFPMMSTFKVVLCGAVLAR' \
                'VDAGDEQLERKIHYRQQDLVDYSPVSEKHLADGMTVGELCAAAITMSDNSAANLLLATVGGPAGLTAFLRQIGDNVTRL' \
                'DRWETELNEALPGDARDTTTPASMAATLRKLLTSQRLSARSQRQLLQWMVDDRVAGPLIRSVLPAGWFIADKTGAGERG' \
                'ARGIVALLGPNNKAERIVVIYLRDTPASMAERNQQIAGIGAALIEHWQR'

    # define the aligner
    protein_aligner = Align.PairwiseAligner()
    protein_aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    protein_aligner.open_gap_score = -10
    protein_aligner.extend_gap_score = -0.5


    alignments = protein_aligner.align(shv_1_ref, translation)

    # If we didn't get any global amino acid alignments, then it's not appropriate to look for SHV
    # mutations in this hit.
    if not alignments:
        return bla_class, [], [], None

    alignment = alignments[0]
    ref_aligned, hit_aligned = alignment
    score = alignment.score
    

    # If the identity of the alignment is too low, then it's not appropriate to look for SHV
    # mutations in this hit.
    identity = get_percent_identity(ref_aligned, hit_aligned)
    
    if identity < 0.9:
        return bla_class, [], [], None
    
    
    # Mutations at these sites will lead to an ESBL class:
    pos_169_mut, pos_169_aa = get_mut(ref_aligned, hit_aligned, 164, 169, 'L')
    pos_179_mut, pos_179_aa = get_mut(ref_aligned, hit_aligned, 174, 179, 'D')
    pos_238_mut, pos_238_aa = get_mut(ref_aligned, hit_aligned, 233, 238, 'G')
    pos_148_mut, pos_148_aa = get_mut(ref_aligned, hit_aligned, 143, 148, 'L')
    
    # Mutations at site Ambler-240 will lead to an ESBL class, but only if site Ambler-35 is also mutated.
    pos_035_mut, pos_035_aa = get_mut(ref_aligned, hit_aligned,  30,  35, 'L')
    pos_240_mut, pos_240_aa = get_mut(ref_aligned, hit_aligned, 234, 240, 'E')

    has_esbl = (pos_169_mut or pos_179_mut or pos_238_mut or pos_148_mut or 
                                            (pos_240_mut and pos_035_mut))

    esbl_mutations = [pos_169_mut, pos_179_mut, pos_238_mut, pos_148_mut, pos_240_mut]

    # Mutations at these sites will lead to inhibition:
    pos_069_mut, pos_069_aa = get_mut(ref_aligned, hit_aligned,  64,  69, 'M')
    pos_130_mut, pos_130_aa = get_mut(ref_aligned, hit_aligned, 125, 130, 'S')
    pos_234_mut, pos_234_aa = get_mut(ref_aligned, hit_aligned, 229, 234, 'K')
    pos_235_mut, pos_235_aa = get_mut(ref_aligned, hit_aligned, 230, 235, 'T')

    has_inhr = (pos_069_mut or pos_130_mut or pos_234_mut or pos_235_mut)
    inhr_mutations = [pos_069_mut, pos_130_mut, pos_234_mut, pos_235_mut]

    # Mutations at these sites don't change the class, but will still be reported:
    pos_025_mut, pos_025_aa = get_mut(ref_aligned, hit_aligned,  20,  25, 'A')
    pos_035_mut, pos_035_aa = get_mut(ref_aligned, hit_aligned,  30,  35, 'L')
    pos_146_mut, pos_146_aa = get_mut(ref_aligned, hit_aligned, 141, 146, 'A')
    pos_156_mut, pos_156_aa = get_mut(ref_aligned, hit_aligned, 151, 156, 'G')

    # Mutations in the omega loop, for tracking only (not class modification):
    pos_164_mut, pos_164_aa = get_mut(ref_aligned, hit_aligned, 159, 164, 'R')
    pos_165_mut, pos_165_aa = get_mut(ref_aligned, hit_aligned, 160, 165, 'W')
    pos_166_mut, pos_166_aa = get_mut(ref_aligned, hit_aligned, 161, 166, 'E')
    pos_167_mut, pos_167_aa = get_mut(ref_aligned, hit_aligned, 162, 167, 'T')
    pos_168_mut, pos_168_aa = get_mut(ref_aligned, hit_aligned, 163, 168, 'E')
    pos_170_mut, pos_170_aa = get_mut(ref_aligned, hit_aligned, 165, 170, 'N')
    pos_171_mut, pos_171_aa = get_mut(ref_aligned, hit_aligned, 166, 171, 'E')
    pos_172_mut, pos_172_aa = get_mut(ref_aligned, hit_aligned, 167, 172, 'A')
    pos_173_mut, pos_173_aa = get_mut(ref_aligned, hit_aligned, 168, 173, 'L')
    pos_174_mut, pos_174_aa = get_mut(ref_aligned, hit_aligned, 169, 174, 'P')
    pos_175_mut, pos_175_aa = get_mut(ref_aligned, hit_aligned, 170, 175, 'G')
    pos_176_mut, pos_176_aa = get_mut(ref_aligned, hit_aligned, 171, 176, 'D')
    pos_177_mut, pos_177_aa = get_mut(ref_aligned, hit_aligned, 172, 177, 'A')
    pos_178_mut, pos_178_aa = get_mut(ref_aligned, hit_aligned, 173, 178, 'R')

    omega_loop_seq = ''.join([pos_164_aa, pos_165_aa, pos_166_aa, pos_167_aa, pos_168_aa,
                              pos_169_aa, pos_170_aa, pos_171_aa, pos_172_aa, pos_173_aa,
                              pos_174_aa, pos_175_aa, pos_176_aa, pos_177_aa, pos_178_aa,
                              pos_179_aa])
                              
    if omega_loop_seq == 'RWETELNEALPGDARD':  # if it's the same as SHV-1
        omega_loop_seq = None

    shv_mutations = [pos_025_mut, pos_035_mut, pos_069_mut, pos_130_mut, pos_146_mut, pos_148_mut,
                     pos_156_mut, pos_164_mut, pos_165_mut, pos_166_mut, pos_167_mut, pos_168_mut,
                     pos_169_mut, pos_170_mut, pos_171_mut, pos_172_mut, pos_173_mut, pos_174_mut,
                     pos_175_mut, pos_176_mut, pos_177_mut, pos_178_mut, pos_179_mut, pos_234_mut,
                     pos_235_mut, pos_238_mut, pos_240_mut]
                     
    shv_mutations = [m for m in shv_mutations if m]

    if exact_match:
        new_bla_class = bla_class
    else:
        new_bla_class = get_new_bla_class(has_esbl, has_inhr)
    class_changing_mutations = \
        get_class_changing_mutations(bla_class, new_bla_class, esbl_mutations, inhr_mutations)

    return new_bla_class, shv_mutations, class_changing_mutations, omega_loop_seq


def get_percent_identity(ref_aligned, hit_aligned):
    matches = 0
    for i, a in enumerate(ref_aligned):
        b = hit_aligned[i]
        if a == b:
            matches += 1
    return matches / len(ref_aligned)



def get_mut(ref_aligned, hit_aligned, ref_pos, ambler_pos, ref_aa):
    """
    ref_pos:    the index of the AA in SHV-1 (0-based because we're looking it up in Python)
    ambler_pos: the index of the AA in the Ambler alignment (1-based because it's just for the name)
    ref_aa:     the AA in the SHV-1 sequence
    """
    ref_no_gaps, hit_no_gaps = [], []
    for i, a in enumerate(ref_aligned):
        b = hit_aligned[i]
        if a != '-':
            ref_no_gaps.append(a)
            hit_no_gaps.append(b)
    assert len(ref_no_gaps) == 286
    assert ref_no_gaps[ref_pos] == ref_aa
    hit_aa = hit_no_gaps[ref_pos]
    if ref_aa != hit_aa and hit_aa != '-':
        mutation_notation = f'{ambler_pos}{hit_aa}'
    else:
        mutation_notation = ''
    return mutation_notation, hit_aa


def get_new_bla_class(esbl, inhr):
    if not esbl and not inhr:
        return 'Bla_chr'
    if esbl and not inhr:
        return 'Bla_ESBL'
    if inhr and not esbl:
        return 'Bla_inhR'
    if esbl and inhr:
        return 'Bla_ESBL_inhR'


def get_class_changing_mutations(bla_class, new_bla_class, esbl_mutations, inhr_mutations):
    if bla_class == new_bla_class:
        return []
    class_changing_mutations = []
    if ('ESBL' in bla_class and 'ESBL' not in new_bla_class) or \
            ('ESBL' not in bla_class and 'ESBL' in new_bla_class):
        class_changing_mutations += esbl_mutations
    if ('inhR' in bla_class and 'inhR' not in new_bla_class) or \
            ('inhR' not in bla_class and 'inhR' in new_bla_class):
        class_changing_mutations += inhr_mutations
    return [m for m in class_changing_mutations if m]
