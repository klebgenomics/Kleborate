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
from kleborate.kleborate import get_resource_paths, get_chromosome_mlst_results


class TestMlst(unittest.TestCase):
    """
    Tests chromosomal MLST calls.
    """

    def setUp(self):
        self.data_folder, self.mlstblast, self.resblast, self.clusterblast = get_resource_paths()

    def test_chromosome_random(self):
        """
        This test has just random sequence and should give no MLST call.
        """
        results = get_chromosome_mlst_results(self.mlstblast, self.data_folder,
                                              'test/test_random.fasta')
        self.assertEqual(results['gapA'], '-')
        self.assertEqual(results['infB'], '-')
        self.assertEqual(results['mdh'], '-')
        self.assertEqual(results['pgi'], '-')
        self.assertEqual(results['phoE'], '-')
        self.assertEqual(results['rpoB'], '-')
        self.assertEqual(results['tonB'], '-')
        self.assertEqual(results['ST'], '0')

    def test_chromosome_exact(self):
        """
        This test is an exact match for ST23.
        """
        results = get_chromosome_mlst_results(self.mlstblast, self.data_folder,
                                              'test/test_mlst_1.fasta')
        self.assertEqual(results['gapA'], '2')
        self.assertEqual(results['infB'], '1')
        self.assertEqual(results['mdh'], '1')
        self.assertEqual(results['pgi'], '1')
        self.assertEqual(results['phoE'], '9')
        self.assertEqual(results['rpoB'], '4')
        self.assertEqual(results['tonB'], '12')
        self.assertEqual(results['ST'], 'ST23')

    def test_chromosome_inexact(self):
        """
        This test is is one base off from ST23 (single substitution in mdh_1 of G to A).
        """
        results = get_chromosome_mlst_results(self.mlstblast, self.data_folder,
                                              'test/test_mlst_2.fasta')
        self.assertEqual(results['gapA'], '2')
        self.assertEqual(results['infB'], '1')
        self.assertEqual(results['mdh'], '1*')
        self.assertEqual(results['pgi'], '1')
        self.assertEqual(results['phoE'], '9')
        self.assertEqual(results['rpoB'], '4')
        self.assertEqual(results['tonB'], '12')
        self.assertEqual(results['ST'], 'ST23-1LV')

    def test_unknown_ST(self):
        """
        This test has exact matches for alleles but in an unknown ST combination
        """
        results = get_chromosome_mlst_results(self.mlstblast, self.data_folder,
                                              'test/test_mlst_3.fasta')
        self.assertEqual(results['gapA'], '44')
        self.assertEqual(results['infB'], '56')
        self.assertEqual(results['mdh'], '34')
        self.assertEqual(results['pgi'], '19')
        self.assertEqual(results['phoE'], '30')
        self.assertEqual(results['rpoB'], '53')
        self.assertEqual(results['tonB'], '55')
        self.assertEqual(results['ST'], '0')
