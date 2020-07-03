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
from kleborate.__main__ import get_data_path, get_output_headers, get_summary_results


class TestResScore(unittest.TestCase):

    def setUp(self):
        self.data_dir = get_data_path()
        Args = collections.namedtuple('Args', ['resistance', 'kaptive_k', 'kaptive_o'])
        self.args = Args(resistance=True, kaptive_k=False, kaptive_o=False)
        _, _, self.res_headers = get_output_headers(self.args, self.data_dir)
        self.results = {res_class: '-' for res_class in self.res_headers}
        self.results['Yersiniabactin'] = '-'
        self.results['Colibactin'] = '-'
        self.results['Aerobactin'] = '-'

    def test_nothing(self):
        summary_results = get_summary_results(self.results, self.res_headers)
        self.assertEqual(summary_results['virulence_score'], '0')
        self.assertEqual(summary_results['resistance_score'], '0')
        self.assertEqual(summary_results['num_resistance_classes'], '0')
        self.assertEqual(summary_results['num_resistance_genes'], '0')

    def test_res_counts_1(self):
        self.results['AGly'] = 'a'
        self.results['Flq'] = 'b'
        self.results['Tet'] = 'c'
        summary_results = get_summary_results(self.results, self.res_headers)
        self.assertEqual(summary_results['num_resistance_classes'], '3')
        self.assertEqual(summary_results['num_resistance_genes'], '3')

    def test_res_counts_2(self):
        self.results['AGly'] = 'a'
        self.results['Flq'] = 'b;c'
        self.results['Tet'] = 'd;e;f'
        summary_results = get_summary_results(self.results, self.res_headers)
        self.assertEqual(summary_results['num_resistance_classes'], '3')
        self.assertEqual(summary_results['num_resistance_genes'], '6')

    def test_res_counts_3(self):
        """
        Intrinsic Bla genes should not add to the counts.
        """
        self.results['AGly'] = 'a'
        self.results['Flq'] = 'b;c'
        self.results['Tet'] = 'd;e;f'
        self.results['Bla'] = 'g;SHV-3'
        summary_results = get_summary_results(self.results, self.res_headers)
        self.assertEqual(summary_results['num_resistance_classes'], '3')
        self.assertEqual(summary_results['num_resistance_genes'], '7')

    def test_res_counts_4(self):
        """
        Bla genes in columns other than 'Bla' (e.g. 'Bla_ESBL') should add to the counts.
        """
        self.results['Tet'] = 'd;e;f'
        self.results['Bla'] = 'g;SHV-3'
        self.results['Bla_ESBL'] = 'SHV-4'
        summary_results = get_summary_results(self.results, self.res_headers)
        self.assertEqual(summary_results['num_resistance_classes'], '2')
        self.assertEqual(summary_results['num_resistance_genes'], '5')

    def test_res_counts_5(self):
        """
        Bla genes in columns other than 'Bla' (e.g. 'Bla_broad') should add to the counts.
        """
        self.results['Tet'] = 'd;e;f'
        self.results['Bla_broad'] = 'SHV-5'
        summary_results = get_summary_results(self.results, self.res_headers)
        self.assertEqual(summary_results['num_resistance_classes'], '2')
        self.assertEqual(summary_results['num_resistance_genes'], '4')

    def test_res_counts_6(self):
        """
        Omp genes should not add to the counts.
        """
        self.results['AGly'] = 'a'
        self.results['Flq'] = 'b;c'
        self.results['Tet'] = 'd;e;f'
        self.results['Omp_truncations'] = 'g;h'
        summary_results = get_summary_results(self.results, self.res_headers)
        self.assertEqual(summary_results['num_resistance_classes'], '3')
        self.assertEqual(summary_results['num_resistance_genes'], '6')

    def test_res_counts_7(self):
        """
        Mutations should add to the class count but not to the gene count (because that is for
        acquired resistance genes).
        """
        self.results['AGly'] = 'a'
        self.results['Flq'] = 'GyrA-83I'
        self.results['Tet'] = 'd;e;f'
        summary_results = get_summary_results(self.results, self.res_headers)
        self.assertEqual(summary_results['num_resistance_classes'], '3')
        self.assertEqual(summary_results['num_resistance_genes'], '4')

    def test_res_counts_8(self):
        """
        Mutations should add to the class count but not to the gene count (because that is for
        acquired resistance genes).
        """
        self.results['AGly'] = 'a'
        self.results['Flq'] = 'b;ParC-80I'
        self.results['Tet'] = 'd;e;f'
        summary_results = get_summary_results(self.results, self.res_headers)
        self.assertEqual(summary_results['num_resistance_classes'], '3')
        self.assertEqual(summary_results['num_resistance_genes'], '5')

    def test_res_counts_9(self):
        """
        Truncations should add to the class count but not to the gene count (because that is for
        acquired resistance genes).
        """
        self.results['AGly'] = 'a'
        self.results['Tet'] = 'b;c'
        self.results['Col'] = 'MgrB-20%'
        summary_results = get_summary_results(self.results, self.res_headers)
        self.assertEqual(summary_results['num_resistance_classes'], '3')
        self.assertEqual(summary_results['num_resistance_genes'], '3')

    def test_res_counts_10(self):
        """
        Truncations should add to the class count but not to the gene count (because that is for
        acquired resistance genes).
        """
        self.results['AGly'] = 'a'
        self.results['Tet'] = 'b;c'
        self.results['Col'] = 'd;PmrB-20%'
        summary_results = get_summary_results(self.results, self.res_headers)
        self.assertEqual(summary_results['num_resistance_classes'], '3')
        self.assertEqual(summary_results['num_resistance_genes'], '4')

    def test_res_counts_11(self):
        """
        Spurious hits should not add to either the gene or class counts.
        """
        self.results['AGly'] = 'a'
        self.results['Flq'] = 'b'
        self.results['Tet'] = 'c'
        self.results['spurious_resistance_hits'] = 'd;e;f'
        summary_results = get_summary_results(self.results, self.res_headers)
        self.assertEqual(summary_results['num_resistance_classes'], '3')
        self.assertEqual(summary_results['num_resistance_genes'], '3')

    def test_res_score_1(self):
        self.results['Bla'] = 'a'
        summary_results = get_summary_results(self.results, self.res_headers)
        self.assertEqual(summary_results['resistance_score'], '0')

    def test_res_score_2(self):
        self.results['Bla'] = 'a'
        self.results['Bla_ESBL'] = 'b'
        summary_results = get_summary_results(self.results, self.res_headers)
        self.assertEqual(summary_results['resistance_score'], '1')

    def test_res_score_3(self):
        self.results['Bla'] = 'a'
        self.results['Bla_ESBL_inhR'] = 'b'
        summary_results = get_summary_results(self.results, self.res_headers)
        self.assertEqual(summary_results['resistance_score'], '1')

    def test_res_score_4(self):
        self.results['Bla'] = 'a'
        self.results['Bla_Carb'] = 'b'
        summary_results = get_summary_results(self.results, self.res_headers)
        self.assertEqual(summary_results['resistance_score'], '2')

    def test_res_score_5(self):
        self.results['Bla'] = 'a'
        self.results['Bla_ESBL'] = 'b'
        self.results['Bla_Carb'] = 'c'
        summary_results = get_summary_results(self.results, self.res_headers)
        self.assertEqual(summary_results['resistance_score'], '2')

    def test_res_score_6(self):
        self.results['Bla'] = 'a'
        self.results['Bla_ESBL'] = 'b'
        self.results['Bla_Carb'] = 'c'
        self.results['Col'] = 'd'
        summary_results = get_summary_results(self.results, self.res_headers)
        self.assertEqual(summary_results['resistance_score'], '3')

    def test_res_score_7(self):
        self.results['Col'] = 'a'
        summary_results = get_summary_results(self.results, self.res_headers)
        self.assertEqual(summary_results['resistance_score'], '0')

    def test_vir_score_1(self):
        self.results['Yersiniabactin'] = 'a'
        summary_results = get_summary_results(self.results, self.res_headers)
        self.assertEqual(summary_results['virulence_score'], '1')

    def test_vir_score_2(self):
        self.results['Colibactin'] = 'a'
        summary_results = get_summary_results(self.results, self.res_headers)
        self.assertEqual(summary_results['virulence_score'], '2')

    def test_vir_score_3(self):
        self.results['Colibactin'] = 'a'
        self.results['Yersiniabactin'] = 'b'
        summary_results = get_summary_results(self.results, self.res_headers)
        self.assertEqual(summary_results['virulence_score'], '2')

    def test_vir_score_4(self):
        self.results['Aerobactin'] = 'a'
        summary_results = get_summary_results(self.results, self.res_headers)
        self.assertEqual(summary_results['virulence_score'], '3')

    def test_vir_score_5(self):
        self.results['Aerobactin'] = 'a'
        self.results['Yersiniabactin'] = 'b'
        summary_results = get_summary_results(self.results, self.res_headers)
        self.assertEqual(summary_results['virulence_score'], '4')

    def test_vir_score_6(self):
        self.results['Aerobactin'] = 'a'
        self.results['Colibactin'] = 'b'
        summary_results = get_summary_results(self.results, self.res_headers)
        self.assertEqual(summary_results['virulence_score'], '5')

    def test_vir_score_7(self):
        self.results['Aerobactin'] = 'a'
        self.results['Colibactin'] = 'b'
        self.results['Yersiniabactin'] = 'c'
        summary_results = get_summary_results(self.results, self.res_headers)
        self.assertEqual(summary_results['virulence_score'], '5')
