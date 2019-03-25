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
from kleborate.kleborate import get_data_path, get_output_headers, get_resistance_results


class Args(object):
    def __init__(self):
        self.resistance = True
        self.species = False
        self.kaptive_k = False
        self.kaptive_o = False


class TestResOmp(unittest.TestCase):
    """
    Tests calling of carbapenem resistance via the OmpK35/OmpK36 genes.
    """

    def setUp(self):
        self.args = Args()
        self.data_dir = get_data_path()
        _, _, self.res_headers = get_output_headers(self.args, self.data_dir)

    def test_both_genes_intact(self):
        results = get_resistance_results(self.data_dir, 'test/sequences/test_res_omp_1.fasta',
                                         self.args, self.res_headers)
        self.assertEqual(results['Bla_Carb'], '-')

    def test_ompk35_frameshift(self):
        """
        A frameshift in OmpK35 should cause an early stop and lead to a carbapenem resistance call.
        """
        results = get_resistance_results(self.data_dir, 'test/sequences/test_res_omp_2.fasta',
                                         self.args, self.res_headers)
        self.assertTrue('OmpK35-' in results['Bla_Carb'])

    def test_ompk35_early_stop(self):
        """
        This tests an early stop mutation (without a frameshift) in OmpK35.
        """
        results = get_resistance_results(self.data_dir, 'test/sequences/test_res_omp_3.fasta',
                                         self.args, self.res_headers)
        self.assertTrue('OmpK35-' in results['Bla_Carb'])

    def test_ompk36_missing(self):
        results = get_resistance_results(self.data_dir, 'test/sequences/test_res_omp_4.fasta',
                                         self.args, self.res_headers)
        self.assertTrue('OmpK36-' in results['Bla_Carb'])

    def test_ompk36gd(self):
        results = get_resistance_results(self.data_dir, 'test/sequences/test_res_omp_5.fasta',
                                         self.args, self.res_headers)
        self.assertTrue('OmpK36GD' in results['Bla_Carb'])

    def test_ompk36td(self):
        results = get_resistance_results(self.data_dir, 'test/sequences/test_res_omp_6.fasta',
                                         self.args, self.res_headers)
        self.assertTrue('OmpK36TD' in results['Bla_Carb'])
