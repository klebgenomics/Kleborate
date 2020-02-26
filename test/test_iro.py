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

import collections
import unittest

from kleborate.kleborate import get_iro_mlst_results


class TestIro(unittest.TestCase):
    """
    Tests salmochelin ST calls.
    """

    def setUp(self):
        self.data_dir = 'test/test_iro/data'
        Args = collections.namedtuple('Args', ['min_coverage', 'min_identity'])
        self.args = Args(min_coverage=80.0, min_identity=90.0)

    def test_iro_random(self):
        """
        This test has just random sequence and should give no iro call.
        """
        results = get_iro_mlst_results(self.data_dir, 'test/test_iro/test_random.fasta', self.args)
        self.assertEqual(results['iroB'], '-')
        self.assertEqual(results['iroC'], '-')
        self.assertEqual(results['iroD'], '-')
        self.assertEqual(results['iroN'], '-')
        self.assertEqual(results['Salmochelin'], '-')
        self.assertEqual(results['SmST'], '0')

    def test_iro_exact_33(self):
        """
        This test is an exact match for SmST33.
        """
        results = get_iro_mlst_results(self.data_dir, 'test/test_iro/test_iro_1.fasta', self.args)
        self.assertEqual(results['iroB'], '22-29%')
        self.assertEqual(results['iroC'], '33-4%')
        self.assertEqual(results['iroD'], '12')
        self.assertEqual(results['iroN'], '15')
        self.assertEqual(results['Salmochelin'], 'iro 4')
        self.assertEqual(results['SmST'], '33')

    def test_iro_inexact(self):
        """
        This test is an inexact match for SmST33. There are single base changes in iroC and in
        iroN.
        """
        results = get_iro_mlst_results(self.data_dir, 'test/test_iro/test_iro_2.fasta', self.args)
        self.assertEqual(results['iroB'], '22-29%')
        self.assertEqual(results['iroC'], '33*-4%')
        self.assertEqual(results['iroD'], '12')
        self.assertEqual(results['iroN'], '15*')
        self.assertEqual(results['Salmochelin'], 'iro 4')
        self.assertEqual(results['SmST'], '33-2LV')

    def test_iro_novel(self):
        """
        This test is an exact match for alleles, but an unknown combination.
        """
        results = get_iro_mlst_results(self.data_dir, 'test/test_iro/test_iro_3.fasta', self.args)
        self.assertEqual(results['iroB'], '11')
        self.assertEqual(results['iroC'], '1')
        self.assertEqual(results['iroD'], '18')
        self.assertEqual(results['iroN'], '5')
        self.assertEqual(results['Salmochelin'], 'iro unknown')
        self.assertEqual(results['SmST'], '0')

    def test_iro_incomplete(self):
        """
        This test is an exact match for only one allele.
        """
        results = get_iro_mlst_results(self.data_dir, 'test/test_iro/test_iro_4.fasta', self.args)
        self.assertEqual(results['iroB'], '-')
        self.assertEqual(results['iroC'], '1')
        self.assertEqual(results['iroD'], '-')
        self.assertEqual(results['iroN'], '-')
        self.assertEqual(results['Salmochelin'], '-')
        self.assertEqual(results['SmST'], '0')

    def test_iro_exact_1(self):
        """
        This test is an exact match for SmST1.
        """
        results = get_iro_mlst_results(self.data_dir, 'test/test_iro/test_iro_5.fasta', self.args)
        self.assertEqual(results['iroB'], '1')
        self.assertEqual(results['iroC'], '1')
        self.assertEqual(results['iroD'], '1')
        self.assertEqual(results['iroN'], '1')
        self.assertEqual(results['Salmochelin'], 'iro 1')
        self.assertEqual(results['SmST'], '1')

    def test_iro_missing_1(self):
        """
        This test is a match for SmST1, except that it's missing iroN.
        """
        results = get_iro_mlst_results(self.data_dir, 'test/test_iro/test_iro_6.fasta', self.args)
        self.assertEqual(results['iroB'], '1')
        self.assertEqual(results['iroC'], '1')
        self.assertEqual(results['iroD'], '1')
        self.assertEqual(results['iroN'], '-')
        self.assertEqual(results['Salmochelin'], 'iro 1 (incomplete)')
        self.assertEqual(results['SmST'], '1-1LV')
