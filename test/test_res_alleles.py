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

from kleborate.__main__ import get_output_headers, get_resistance_results


class TestResAlleles(unittest.TestCase):
    """
    Tests calling of resistance via alleles.
    """

    def setUp(self):
        self.data_dir = 'test/test_res_alleles/data'
        Args = collections.namedtuple('Args', ['resistance', 'kaptive_k', 'kaptive_o',
                                               'min_coverage', 'min_identity',
                                               'min_spurious_coverage', 'min_spurious_identity'])
        self.args = Args(resistance=True, kaptive_k=False, kaptive_o=False,
                         min_coverage=80.0, min_identity=90.0,
                         min_spurious_coverage=40.0, min_spurious_identity=80.0)
        _, _, self.res_headers = get_output_headers(self.args, self.data_dir)

    def test_res_01(self):
        results = get_resistance_results(self.data_dir, 'test/test_res_alleles/01.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Tet'], '-')
        self.assertEqual(results['Bla'], 'ABC-1')
        self.assertEqual(results['Bla_ESBL'], '-')
        self.assertEqual(results['spurious_resistance_hits'], '-')

    def test_res_02(self):
        results = get_resistance_results(self.data_dir, 'test/test_res_alleles/02.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Tet'], '-')
        self.assertEqual(results['Bla'], 'ABC-2')
        self.assertEqual(results['Bla_ESBL'], '-')
        self.assertEqual(results['spurious_resistance_hits'], '-')

    def test_res_03(self):
        results = get_resistance_results(self.data_dir, 'test/test_res_alleles/03.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Tet'], '-')
        self.assertEqual(results['Bla'], '-')
        self.assertEqual(results['Bla_ESBL'], 'ABC-3')
        self.assertEqual(results['spurious_resistance_hits'], '-')

    def test_res_04(self):
        results = get_resistance_results(self.data_dir, 'test/test_res_alleles/04.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Tet'], '-')
        self.assertEqual(results['Bla'], '-')
        self.assertEqual(results['Bla_ESBL'], 'ABC-4')
        self.assertEqual(results['spurious_resistance_hits'], '-')

    def test_res_05(self):
        results = get_resistance_results(self.data_dir, 'test/test_res_alleles/05.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Tet'], '-')
        self.assertTrue(results['Bla'] == 'ABC-1^' or results['Bla'] == 'ABC-2^')
        self.assertEqual(results['Bla_ESBL'], '-')
        self.assertEqual(results['spurious_resistance_hits'], '-')

    def test_res_06(self):
        results = get_resistance_results(self.data_dir, 'test/test_res_alleles/06.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Tet'], '-')
        self.assertEqual(results['Bla'], '-')
        self.assertTrue(results['Bla_ESBL'] == 'ABC-3^' or results['Bla_ESBL'] == 'ABC-4^')
        self.assertEqual(results['spurious_resistance_hits'], '-')

    def test_res_07(self):
        results = get_resistance_results(self.data_dir, 'test/test_res_alleles/07.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Tet'], '-')
        self.assertEqual(results['Bla'], 'ABC-2*')
        self.assertEqual(results['Bla_ESBL'], '-')
        self.assertEqual(results['spurious_resistance_hits'], '-')

    def test_res_08(self):
        results = get_resistance_results(self.data_dir, 'test/test_res_alleles/08.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Tet'], '-')
        self.assertEqual(results['Bla'], 'ABC-1;ABC-2')
        self.assertEqual(results['Bla_ESBL'], '-')
        self.assertEqual(results['spurious_resistance_hits'], '-')

    def test_res_09(self):
        results = get_resistance_results(self.data_dir, 'test/test_res_alleles/09.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Tet'], '-')
        self.assertTrue(results['Bla'] == 'ABC-1^' or results['Bla'] == 'ABC-2^')
        self.assertTrue(results['Bla_ESBL'] == 'ABC-3^' or results['Bla_ESBL'] == 'ABC-4^')
        self.assertEqual(results['spurious_resistance_hits'], '-')

    def test_res_10(self):
        results = get_resistance_results(self.data_dir, 'test/test_res_alleles/10.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Tet'], '-')
        self.assertEqual(results['Bla'], 'ABC-1;ABC-2')
        self.assertEqual(results['Bla_ESBL'], 'ABC-3*;ABC-4')
        self.assertEqual(results['spurious_resistance_hits'], '-')

    def test_res_11(self):
        """
        This test has the ABC-1 gene but with a stop codon in the middle.
        """
        results = get_resistance_results(self.data_dir, 'test/test_res_alleles/11.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Tet'], '-')
        self.assertEqual(results['Bla'], 'ABC-1*-50%')
        self.assertEqual(results['Bla_ESBL'], '-')
        self.assertEqual(results['spurious_resistance_hits'], '-')

    def test_res_12(self):
        """
        This test has the ABC-1 gene but missing the first (start) codon.
        """
        results = get_resistance_results(self.data_dir, 'test/test_res_alleles/12.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Tet'], '-')
        self.assertEqual(results['Bla'], 'ABC-1?-0%')
        self.assertEqual(results['Bla_ESBL'], '-')
        self.assertEqual(results['spurious_resistance_hits'], '-')

    def test_res_13(self):
        """
        This test has the ABC-1 gene but at only 50% coverage, so it goes to spurious hits.
        """
        results = get_resistance_results(self.data_dir, 'test/test_res_alleles/13.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Tet'], '-')
        self.assertEqual(results['Bla'], '-')
        self.assertEqual(results['Bla_ESBL'], '-')
        self.assertEqual(results['spurious_resistance_hits'], 'ABC-1?-50%')

