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
import unittest

from kleborate.__main__ import get_rmp_mlst_results, get_rmpa2_results


class TestRmp(unittest.TestCase):
    """
    Tests RmpADC ST calls.
    """

    def setUp(self):
        self.data_dir = 'test/test_rmp/data'
        Args = collections.namedtuple('Args', ['min_coverage', 'min_identity'])
        self.args = Args(min_coverage=80.0, min_identity=90.0)

    def test_rmp_random(self):
        """
        This test has just random sequence and should give no rmp call.
        """
        results = get_rmp_mlst_results(self.data_dir, 'test/test_rmp/test_random.fasta', self.args)
        self.assertEqual(results['rmpA'], '-')
        self.assertEqual(results['rmpD'], '-')
        self.assertEqual(results['rmpC'], '-')
        self.assertEqual(results['RmpADC'], '-')
        self.assertEqual(results['RmST'], '0')

    def test_rmp_19(self):
        """
        This test is an exact match for RmST19, an ST without any truncated alleles.
        """
        results = get_rmp_mlst_results(self.data_dir, 'test/test_rmp/test_rmp_1.fasta', self.args)
        self.assertEqual(results['rmpA'], '7')
        self.assertEqual(results['rmpD'], '8')
        self.assertEqual(results['rmpC'], '3')
        self.assertEqual(results['RmpADC'], '2A')
        self.assertEqual(results['RmST'], '19')

    def test_rmp_21(self):
        """
        This test is an exact match for RmST21, an ST with a truncated allele.
        """
        results = get_rmp_mlst_results(self.data_dir, 'test/test_rmp/test_rmp_2.fasta', self.args)
        self.assertEqual(results['rmpA'], '7')
        self.assertEqual(results['rmpD'], '36')
        self.assertEqual(results['rmpC'], '4-17%')
        self.assertEqual(results['RmpADC'], '2A')
        self.assertEqual(results['RmST'], '21')

    def test_rmp_novel(self):
        """
        This test is an exact match for alleles, but an unknown combination.
        """
        results = get_rmp_mlst_results(self.data_dir, 'test/test_rmp/test_rmp_3.fasta', self.args)
        self.assertEqual(results['rmpA'], '8')
        self.assertEqual(results['rmpD'], '1')
        self.assertEqual(results['rmpC'], '9-11%')
        self.assertEqual(results['RmpADC'], 'rmp unknown')
        self.assertEqual(results['RmST'], '0')

    def test_rmp_missing_1(self):
        """
        This test is a match for RmST1, except that it's missing rmpC.
        """
        results = get_rmp_mlst_results(self.data_dir, 'test/test_rmp/test_rmp_4.fasta', self.args)
        self.assertEqual(results['rmpA'], '41-54%')
        self.assertEqual(results['rmpD'], '20')
        self.assertEqual(results['rmpC'], '-')
        self.assertEqual(results['RmpADC'], '2A (incomplete)')
        self.assertEqual(results['RmST'], '1-1LV')

    def test_rmpa2_random(self):
        """
        This test has just random sequence and should give no rmpA2 call.
        """
        results = get_rmpa2_results(self.data_dir, 'test/test_rmp/test_random.fasta', self.args)
        self.assertEqual(results['rmpA2'], '-')

    def test_rmpa2_exact(self):
        """
        This test has an exact match for rmpA2_1.
        """
        results = get_rmpa2_results(self.data_dir, 'test/test_rmp/test_rmpa2_1.fasta', self.args)
        self.assertEqual(results['rmpA2'], 'rmpA2_1')

    def test_rmpa2_exact_truncated(self):
        """
        This test has an exact match for rmpA2_7.
        """
        results = get_rmpa2_results(self.data_dir, 'test/test_rmp/test_rmpa2_2.fasta', self.args)
        self.assertEqual(results['rmpA2'], 'rmpA2_7-47%')

    def test_rmpa2_inexact(self):
        """
        This test has an inexact match for rmpA2_1.
        """
        results = get_rmpa2_results(self.data_dir, 'test/test_rmp/test_rmpa2_3.fasta', self.args)
        self.assertEqual(results['rmpA2'], 'rmpA2_1*')
