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

import unittest

from kleborate.contig_stats import get_contig_stats


class TestContigStats(unittest.TestCase):

    def test_count_1(self):
        contig_count, _, _, _, _ = get_contig_stats('test/test_contig_stats/contig_stats_1.fasta')
        self.assertEqual(contig_count, 4)

    def test_count_2(self):
        contig_count, _, _, _, _ = get_contig_stats('test/test_contig_stats/contig_stats_2.fasta')
        self.assertEqual(contig_count, 3)

    def test_n50_1(self):
        _, n50, _, _, _ = get_contig_stats('test/test_contig_stats/contig_stats_1.fasta')
        self.assertEqual(n50, 40)

    def test_n50_2(self):
        _, n50, _, _, _ = get_contig_stats('test/test_contig_stats/contig_stats_2.fasta')
        self.assertEqual(n50, 200)

    def test_longest_1(self):
        _, _, longest_contig, _, _ = get_contig_stats('test/test_contig_stats/contig_stats_1.fasta')
        self.assertEqual(longest_contig, 45)

    def test_longest_2(self):
        _, _, longest_contig, _, _ = get_contig_stats('test/test_contig_stats/contig_stats_2.fasta')
        self.assertEqual(longest_contig, 200)

    def test_ambiguous_bases_1(self):
        _, _, _, _, ambiguous = get_contig_stats('test/test_contig_stats/contig_stats_1.fasta')
        self.assertEqual(ambiguous, 'no')

    def test_ambiguous_bases_2(self):
        _, _, _, _, ambiguous = get_contig_stats('test/test_contig_stats/contig_stats_2.fasta')
        self.assertEqual(ambiguous, 'yes (1)')

    def test_ambiguous_bases_3(self):
        _, _, _, _, ambiguous = get_contig_stats('test/test_contig_stats/contig_stats_3.fasta')
        self.assertEqual(ambiguous, 'no')

    def test_ambiguous_bases_4(self):
        _, _, _, _, ambiguous = get_contig_stats('test/test_contig_stats/contig_stats_4.fasta')
        self.assertEqual(ambiguous, 'yes (4)')

    def test_total_size_1(self):
        _, _, _, total_size, _ = get_contig_stats('test/test_contig_stats/contig_stats_1.fasta')
        self.assertEqual(total_size, 115)

    def test_total_size_2(self):
        _, _, _, total_size, _ = get_contig_stats('test/test_contig_stats/contig_stats_2.fasta')
        self.assertEqual(total_size, 260)

    def test_total_size_3(self):
        _, _, _, total_size, _ = get_contig_stats('test/test_contig_stats/contig_stats_3.fasta')
        self.assertEqual(total_size, 260)

    def test_total_size_4(self):
        _, _, _, total_size, _ = get_contig_stats('test/test_contig_stats/contig_stats_4.fasta')
        self.assertEqual(total_size, 260)
