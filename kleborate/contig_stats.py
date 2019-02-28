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

import sys
import gzip


def get_contig_stat_results(contigs):
    contig_count, n50, longest_contig = get_contig_stats(contigs)
    return {'contig_count': str(contig_count),
            'N50': str(n50),
            'largest_contig': str(longest_contig)}


def get_contig_stats(assembly):
    """
    Returns various contig length metrics.
    """
    fasta = load_fasta(assembly)
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
    return len(contig_lengths), n50, longest


def load_fasta(filename):
    """Returns the names and sequences for the given fasta file."""
    fasta_seqs = []
    if get_compression_type(filename) == 'gz':
        open_func = gzip.open
    else:  # plain text
        open_func = open
    with open_func(filename, 'rt') as fasta_file:
        name = ''
        sequence = ''
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    fasta_seqs.append((name.split()[0], sequence))
                    sequence = ''
                name = line[1:]
            else:
                sequence += line
        if name:
            fasta_seqs.append((name.split()[0], sequence))
    return fasta_seqs


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)

    unknown_file = open(filename, 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == 'bz2':
        sys.exit('cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('cannot use zip format - use gzip instead')
    return compression_type
