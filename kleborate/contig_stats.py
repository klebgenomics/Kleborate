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

import collections

from .misc import load_fasta


def get_contig_stat_results(contigs, kp_complex):
    contig_count, n50, longest_contig, total_size, ambiguous_bases = get_contig_stats(contigs)
    qc_warnings = get_qc_warnings(total_size, n50, ambiguous_bases, kp_complex)
    return {'contig_count': str(contig_count),
            'N50': str(n50),
            'largest_contig': str(longest_contig),
            'total_size': str(total_size),
            'ambiguous_bases': ambiguous_bases,
            'QC_warnings': qc_warnings}


def get_contig_stats(assembly):
    fasta = load_fasta(assembly)

    base_counts = collections.defaultdict(int)
    for _, seq in fasta:
        for b in seq:
            base_counts[b] += 1
    base_counts.pop('A', None)
    base_counts.pop('C', None)
    base_counts.pop('G', None)
    base_counts.pop('T', None)
    ambiguous_base_count = sum(base_counts.values())
    if ambiguous_base_count:
        ambiguous_bases = 'yes (' + str(ambiguous_base_count) + ')'
    else:
        ambiguous_bases = 'no'

    contig_lengths = sorted([len(x[1]) for x in fasta])
    if not contig_lengths:
        return 0, 0, 0
    longest = contig_lengths[-1]

    total_size = sum(contig_lengths)
    half_total_length = total_size / 2
    total_so_far = 0
    segment_lengths = contig_lengths[::-1]
    for length in segment_lengths:
        total_so_far += length
        if total_so_far >= half_total_length:
            n50 = length
            break
    else:
        n50 = 0

    return len(contig_lengths), n50, longest, total_size, ambiguous_bases


def get_qc_warnings(total_size, n50, ambiguous_bases, kp_complex):
    warnings = []

    # Genome size should be 5 to 6.5 Mbp for KpSC genomes. For non-KpSC genomes, we don't check
    # the genome size.
    if kp_complex:
        if total_size < 5000000:
            warnings.append('total_size')
        elif total_size > 6500000:
            warnings.append('total_size')

    if n50 < 10000:
        warnings.append('N50')

    if 'yes' in ambiguous_bases:
        warnings.append('ambiguous_bases')

    if warnings:
        return ','.join(warnings)
    else:
        return '-'
