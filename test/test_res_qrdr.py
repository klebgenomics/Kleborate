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
from kleborate.kleborate import get_output_headers, get_resistance_results


class Args(object):
    def __init__(self):
        self.resistance = True
        self.kaptive_k = False
        self.kaptive_o = False


class TestResGyrAParC(unittest.TestCase):
    """
    Tests calling of fluoroquinolone resistance via snps in the GyrA and ParC genes.
    """

    def setUp(self):
        self.args = Args()
        self.data_dir = 'test/test_res_qrdr/data'
        _, _, self.res_headers = get_output_headers(self.args, self.data_dir)

    def test_no_mutations(self):
        results = get_resistance_results(self.data_dir, 'test/test_res_qrdr/test_res_qrdr_1.fasta',
                                         self.args, self.res_headers, True)
        self.assertEqual(results['Flq'], '-')

    def test_gyra(self):
        results = get_resistance_results(self.data_dir, 'test/test_res_qrdr/test_res_qrdr_2.fasta',
                                         self.args, self.res_headers, True)
        self.assertTrue('GyrA-83C' in results['Flq'])

    def test_parc(self):
        results = get_resistance_results(self.data_dir, 'test/test_res_qrdr/test_res_qrdr_3.fasta',
                                         self.args, self.res_headers, True)
        self.assertTrue('ParC-84D' in results['Flq'])

    def test_gyra_and_parc(self):
        results = get_resistance_results(self.data_dir, 'test/test_res_qrdr/test_res_qrdr_4.fasta',
                                         self.args, self.res_headers, True)
        self.assertTrue('GyrA-83C' in results['Flq'])
        self.assertTrue('ParC-84D' in results['Flq'])

    def test_gyra_and_parc_non_kp_complex(self):
        """
        Setting the Kp complex variable to False should turn off the QRDR tests.
        """
        results = get_resistance_results(self.data_dir, 'test/test_res_qrdr/test_res_qrdr_4.fasta',
                                         self.args, self.res_headers, False)
        self.assertTrue('GyrA-83C' not in results['Flq'])
        self.assertTrue('ParC-84D' not in results['Flq'])
