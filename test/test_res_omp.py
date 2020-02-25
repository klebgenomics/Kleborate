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

from kleborate.kleborate import get_output_headers, get_resistance_results


class TestResOmp(unittest.TestCase):
    """
    Tests calling of carbapenem resistance via the OmpK35/OmpK36 genes.
    """

    def setUp(self):
        self.data_dir = 'test/test_res_omp/data'
        Args = collections.namedtuple('Args', ['resistance', 'kaptive_k', 'kaptive_o',
                                               'min_coverage', 'min_identity',
                                               'min_spurious_coverage', 'min_spurious_identity'])
        self.args = Args(resistance=True, kaptive_k=False, kaptive_o=False,
                         min_coverage=80.0, min_identity=90.0,
                         min_spurious_coverage=40.0, min_spurious_identity=80.0)
        _, _, self.res_headers = get_output_headers(self.args, self.data_dir)

    def test_both_genes_intact(self):
        results = get_resistance_results(self.data_dir, 'test/test_res_omp/test_res_omp_1.fasta',
                                         self.args, self.res_headers, True)
        self.assertEqual(results['Omp'], '-')

    def test_ompk35_frameshift(self):
        """
        A frameshift in OmpK35 should cause an early stop and lead to a carbapenem resistance call.
        """
        results = get_resistance_results(self.data_dir, 'test/test_res_omp/test_res_omp_2.fasta',
                                         self.args, self.res_headers, True)
        self.assertTrue('OmpK35-' in results['Omp'])

    def test_ompk35_early_stop(self):
        """
        This tests an early stop mutation (without a frameshift) in OmpK35.
        """
        results = get_resistance_results(self.data_dir, 'test/test_res_omp/test_res_omp_3.fasta',
                                         self.args, self.res_headers, True)
        self.assertTrue('OmpK35-' in results['Omp'])

    def test_ompk36_missing(self):
        results = get_resistance_results(self.data_dir, 'test/test_res_omp/test_res_omp_4.fasta',
                                         self.args, self.res_headers, True)
        self.assertTrue('OmpK36-' in results['Omp'])

    def test_ompk36gd(self):
        results = get_resistance_results(self.data_dir, 'test/test_res_omp/test_res_omp_5.fasta',
                                         self.args, self.res_headers, True)
        self.assertTrue('OmpK36GD' in results['Omp'])

    def test_ompk36td(self):
        results = get_resistance_results(self.data_dir, 'test/test_res_omp/test_res_omp_6.fasta',
                                         self.args, self.res_headers, True)
        self.assertTrue('OmpK36TD' in results['Omp'])

    def test_ompk36td_non_kp_complex(self):
        """
        Setting the Kp complex variable to False should turn off the OmpK tests.
        """
        results = get_resistance_results(self.data_dir, 'test/test_res_omp/test_res_omp_6.fasta',
                                         self.args, self.res_headers, False)
        self.assertTrue('OmpK36TD' not in results['Omp'])
