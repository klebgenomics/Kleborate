"""
Copyright 2018 Kat Holt
Copyright 2018 Ryan Wick (rrwick@gmail.com)
https://github.com/katholt/Kleborate/

This file is part of Kleborate. Kleborate is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Kleborate is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Kleborate. If
not, see <http://www.gnu.org/licenses/>.
"""

from .misc import load_fasta


def get_contig_stat_results(contigs):
    contig_count, n50, longest_contig, ambiguous_bases = get_contig_stats(contigs)
    return {'contig_count': str(contig_count),
            'N50': str(n50),
            'largest_contig': str(longest_contig),
            'ambiguous_bases': ambiguous_bases}


def get_contig_stats(assembly):
    """
    Returns various contig length metrics.
    """
    fasta = load_fasta(assembly)

    characters = set()
    for _, seq in fasta:
        characters |= set(b for b in seq)
    characters -= {'A', 'C', 'G', 'T'}
    if characters:
        ambiguous_bases = 'yes'
    else:
        ambiguous_bases = 'no'

    contig_lengths = sorted([len(x[1]) for x in fasta])
    if not contig_lengths:
        return 0, 0, 0
    longest = contig_lengths[-1]

    half_total_length = sum(contig_lengths) / 2
    total_so_far = 0
    segment_lengths = contig_lengths[::-1]
    for length in segment_lengths:
        total_so_far += length
        if total_so_far >= half_total_length:
            n50 = length
            break
    else:
        n50 = 0

    return len(contig_lengths), n50, longest, ambiguous_bases
