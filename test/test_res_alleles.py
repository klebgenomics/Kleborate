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
import tempfile
import unittest

from kleborate.kleborate import get_data_path, gunzip_contigs_if_necessary, get_output_headers, \
    get_resistance_results


class TestResAlleles(unittest.TestCase):
    """
    Tests calling of resistance via alleles.
    """

    def setUp(self):
        self.data_dir = 'test/res_test/data'
        Args = collections.namedtuple('Args', ['species', 'resistance', 'kaptive_k', 'kaptive_o'])
        self.args = Args(species=False, resistance=True, kaptive_k=False, kaptive_o=False)
        _, _, self.res_headers = get_output_headers(self.args, self.data_dir)

    def test_res_01(self):
        results = get_resistance_results(self.data_dir, 'test/res_test/01.fasta', self.args,
                                         self.res_headers)
        self.assertEqual(results['Flq'], 'ABC-1')
        self.assertEqual(results['Tet'], '-')
        self.assertEqual(results['Bla'], '-')
        self.assertEqual(results['Bla_Carb'], '-')
        self.assertEqual(results['Bla_ESBL'], '-')
        self.assertEqual(results['Bla_broad'], '-')

    def test_res_02(self):
        results = get_resistance_results(self.data_dir, 'test/res_test/02.fasta', self.args,
                                         self.res_headers)
        self.assertEqual(results['Flq'], 'ABC-2')
        self.assertEqual(results['Tet'], '-')
        self.assertEqual(results['Bla'], '-')
        self.assertEqual(results['Bla_Carb'], '-')
        self.assertEqual(results['Bla_ESBL'], '-')
        self.assertEqual(results['Bla_broad'], '-')

    def test_res_03(self):
        """
        This query allele is equally distant from both ABC-1 and ABC-2 (two substitutions), but an
        exact amino acid match for each. While both alleles are an acceptable answer, we should
        return ABC-1 because it's the lower numbered allele.
        """
        results = get_resistance_results(self.data_dir, 'test/res_test/03.fasta', self.args,
                                         self.res_headers)
        self.assertEqual(results['Flq'], 'ABC-1*')
        self.assertEqual(results['Tet'], '-')
        self.assertEqual(results['Bla'], '-')
        self.assertEqual(results['Bla_Carb'], '-')
        self.assertEqual(results['Bla_ESBL'], '-')
        self.assertEqual(results['Bla_broad'], '-')

    def test_res_04(self):
        """
        This query allele has no perfect nucleotide matches but it matches ABC-3's amino acids.
        """
        results = get_resistance_results(self.data_dir, 'test/res_test/04.fasta', self.args,
                                         self.res_headers)
        self.assertEqual(results['Flq'], 'ABC-3^')
        self.assertEqual(results['Tet'], '-')
        self.assertEqual(results['Bla'], '-')
        self.assertEqual(results['Bla_Carb'], '-')
        self.assertEqual(results['Bla_ESBL'], '-')
        self.assertEqual(results['Bla_broad'], '-')
