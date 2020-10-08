"""
Copyright 2020 Kat Holt
Copyright 2020 Ryan Wick (rrwick@gmail.com)
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


def truncation_check(hit, cov_threshold=90.0):
    """
    Checks to see if the gene is truncated at the amino acid level.
    """
    nucl_seq, ref_start, _ = hit.get_seq_start_end_pos_strand()

    # The hit must start at the first base of the gene. If not, the gene is considered 0%.
    if ref_start != 1:
        return '-0%', 0.0, ''

    # If there are any ambiguous bases in the sequence, then they will break translation, probably
    # resulting in truncation call.
    ambiguous_bases = set(b for b in nucl_seq) - {'A', 'C', 'G', 'T'}
    for b in ambiguous_bases:
        nucl_seq = nucl_seq.split(b)[0]

    # BioPython doesn't like it if the sequence isn't a multiple of 3.
    nucl_seq = nucl_seq[:len(nucl_seq) // 3 * 3]

    # The assumption is that the reference allele is a full CDS with a stop codon at the end. This
    # isn't always true (the reference sequence is sometimes broken) but will serve to make our
    # denominator for coverage.
    ref_aa_length = (hit.ref_length - 3) // 3

    coding_dna = Seq(nucl_seq)
    translation = str(coding_dna.translate(table='Bacterial', to_stop=True))

    coverage = 100.0 * len(translation) / ref_aa_length
    if coverage >= cov_threshold:
        return '', coverage, translation
    else:
        return '-{:.0f}%'.format(coverage), coverage, translation
