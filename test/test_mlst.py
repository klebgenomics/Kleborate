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
import tempfile
import unittest

from kleborate.__main__ import get_chromosome_mlst_results, gunzip_contigs_if_necessary, \
    get_kp_subspecies_based_on_st


class TestMlst(unittest.TestCase):
    """
    Tests chromosomal MLST calls.
    """

    def setUp(self):
        self.data_dir = 'test/test_mlst/data'
        Args = collections.namedtuple('Args', ['min_coverage', 'min_identity'])
        self.args = Args(min_coverage=80.0, min_identity=90.0)

    def test_chromosome_random(self):
        """
        This test has just random sequence and should give no MLST call.
        """
        results = get_chromosome_mlst_results(self.data_dir,
                                              'test/test_mlst/test_random.fasta', True, self.args)
        self.assertEqual(results['gapA'], '-')
        self.assertEqual(results['infB'], '-')
        self.assertEqual(results['mdh'], '-')
        self.assertEqual(results['pgi'], '-')
        self.assertEqual(results['phoE'], '-')
        self.assertEqual(results['rpoB'], '-')
        self.assertEqual(results['tonB'], '-')
        self.assertEqual(results['ST'], '0')
        self.assertEqual(results['Chr_ST'], '0')

    def test_chromosome_exact(self):
        """
        This test is an exact match for ST23.
        """
        results = get_chromosome_mlst_results(self.data_dir,
                                              'test/test_mlst/test_mlst_1.fasta', True, self.args)
        self.assertEqual(results['gapA'], '2')
        self.assertEqual(results['infB'], '1')
        self.assertEqual(results['mdh'], '1')
        self.assertEqual(results['pgi'], '1')
        self.assertEqual(results['phoE'], '9')
        self.assertEqual(results['rpoB'], '4')
        self.assertEqual(results['tonB'], '12')
        self.assertEqual(results['ST'], 'ST23')
        self.assertEqual(results['Chr_ST'], 'ST23')

    def test_chromosome_inexact(self):
        """
        This test is is one base off from ST23 (single substitution in mdh_1 of G to A).
        """
        results = get_chromosome_mlst_results(self.data_dir,
                                              'test/test_mlst/test_mlst_2.fasta', True, self.args)
        self.assertEqual(results['gapA'], '2')
        self.assertEqual(results['infB'], '1')
        self.assertEqual(results['mdh'], '1*')
        self.assertEqual(results['pgi'], '1')
        self.assertEqual(results['phoE'], '9')
        self.assertEqual(results['rpoB'], '4')
        self.assertEqual(results['tonB'], '12')
        self.assertEqual(results['ST'], 'ST23-1LV')
        self.assertEqual(results['Chr_ST'], 'ST23-1LV')

    def test_unknown_ST(self):
        """
        This test has exact matches for alleles but in an unknown ST combination
        """
        results = get_chromosome_mlst_results(self.data_dir,
                                              'test/test_mlst/test_mlst_3.fasta', True, self.args)
        self.assertEqual(results['gapA'], '44')
        self.assertEqual(results['infB'], '56')
        self.assertEqual(results['mdh'], '34')
        self.assertEqual(results['pgi'], '19')
        self.assertEqual(results['phoE'], '30')
        self.assertEqual(results['rpoB'], '53')
        self.assertEqual(results['tonB'], '55')
        self.assertEqual(results['ST'], '0')
        self.assertEqual(results['Chr_ST'], '0')

    def test_83(self):
        contigs = 'test/test_mlst/83.fasta.gz'
        with tempfile.TemporaryDirectory() as tmp_dir:
            contigs = gunzip_contigs_if_necessary(contigs, tmp_dir)
            results = get_chromosome_mlst_results(self.data_dir, contigs, True, self.args)
            self.assertEqual(results['gapA'], '10')
            self.assertEqual(results['infB'], '7')
            self.assertEqual(results['mdh'], '1')
            self.assertEqual(results['pgi'], '1')
            self.assertEqual(results['phoE'], '1')
            self.assertEqual(results['rpoB'], '1')
            self.assertEqual(results['tonB'], '35')
            self.assertEqual(results['ST'], 'ST160')
            self.assertEqual(results['Chr_ST'], 'ST160')

    def test_134(self):
        contigs = 'test/test_mlst/134.fasta.gz'
        with tempfile.TemporaryDirectory() as tmp_dir:
            contigs = gunzip_contigs_if_necessary(contigs, tmp_dir)
            results = get_chromosome_mlst_results(self.data_dir, contigs, True, self.args)
            self.assertEqual(results['gapA'], '2')
            self.assertEqual(results['infB'], '1')
            self.assertEqual(results['mdh'], '2')
            self.assertEqual(results['pgi'], '1')
            self.assertEqual(results['phoE'], '4')
            self.assertEqual(results['rpoB'], '4')
            self.assertEqual(results['tonB'], '4')
            self.assertEqual(results['ST'], 'ST16')
            self.assertEqual(results['Chr_ST'], 'ST16')

    def test_ba779(self):
        contigs = 'test/test_mlst/BA779.fasta.gz'
        with tempfile.TemporaryDirectory() as tmp_dir:
            contigs = gunzip_contigs_if_necessary(contigs, tmp_dir)
            results = get_chromosome_mlst_results(self.data_dir, contigs, True, self.args)
            self.assertEqual(results['gapA'], '2')
            self.assertEqual(results['infB'], '1')
            self.assertEqual(results['mdh'], '1')
            self.assertEqual(results['pgi'], '1')
            self.assertEqual(results['phoE'], '9')
            self.assertEqual(results['rpoB'], '4')
            self.assertEqual(results['tonB'], '12')
            self.assertEqual(results['ST'], 'ST23')
            self.assertEqual(results['Chr_ST'], 'ST23')

    def test_ozaenae(self):
        contigs = 'test/test_mlst/GCF_900451425.1.fna.gz'
        with tempfile.TemporaryDirectory() as tmp_dir:
            contigs = gunzip_contigs_if_necessary(contigs, tmp_dir)
            results = get_chromosome_mlst_results(self.data_dir, contigs, True, self.args)
            self.assertEqual(results['ST'], 'ST90 (subsp. ozaenae)')
            self.assertEqual(results['Chr_ST'], 'ST90')

    def test_rhinoscleromatis(self):
        contigs = 'test/test_mlst/GCF_000163455.1.fna.gz'
        with tempfile.TemporaryDirectory() as tmp_dir:
            contigs = gunzip_contigs_if_necessary(contigs, tmp_dir)
            results = get_chromosome_mlst_results(self.data_dir, contigs, True, self.args)
            self.assertEqual(results['ST'], 'ST67 (subsp. rhinoscleromatis)')
            self.assertEqual(results['Chr_ST'], 'ST67')

    def test_kp_subspecies_1(self):
        # Just test some non-ozaenae non-rhinoscleromatis STs.
        self.assertEqual(get_kp_subspecies_based_on_st('ST1'), 'ST1')
        self.assertEqual(get_kp_subspecies_based_on_st('ST12'), 'ST12')
        self.assertEqual(get_kp_subspecies_based_on_st('ST123'), 'ST123')
        self.assertEqual(get_kp_subspecies_based_on_st('ST1234'), 'ST1234')
        self.assertEqual(get_kp_subspecies_based_on_st('ST1-1LV'), 'ST1-1LV')
        self.assertEqual(get_kp_subspecies_based_on_st('ST12-1LV'), 'ST12-1LV')
        self.assertEqual(get_kp_subspecies_based_on_st('ST123-1LV'), 'ST123-1LV')
        self.assertEqual(get_kp_subspecies_based_on_st('ST1234-1LV'), 'ST1234-1LV')
        self.assertEqual(get_kp_subspecies_based_on_st('ST1-2LV'), 'ST1-2LV')
        self.assertEqual(get_kp_subspecies_based_on_st('ST12-2LV'), 'ST12-2LV')
        self.assertEqual(get_kp_subspecies_based_on_st('ST123-2LV'), 'ST123-2LV')
        self.assertEqual(get_kp_subspecies_based_on_st('ST1234-2LV'), 'ST1234-2LV')

    def test_kp_subspecies_2(self):
        # Test the ozaenae STs.
        self.assertEqual(get_kp_subspecies_based_on_st('ST90'), 'ST90 (subsp. ozaenae)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST91'), 'ST91 (subsp. ozaenae)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST92'), 'ST92 (subsp. ozaenae)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST93'), 'ST93 (subsp. ozaenae)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST95'), 'ST95 (subsp. ozaenae)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST96'), 'ST96 (subsp. ozaenae)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST97'), 'ST97 (subsp. ozaenae)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST381'), 'ST381 (subsp. ozaenae)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST777'), 'ST777 (subsp. ozaenae)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST3193'), 'ST3193 (subsp. ozaenae)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST3766'), 'ST3766 (subsp. ozaenae)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST3768'), 'ST3768 (subsp. ozaenae)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST3771'), 'ST3771 (subsp. ozaenae)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST3781'), 'ST3781 (subsp. ozaenae)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST3782'), 'ST3782 (subsp. ozaenae)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST3784'), 'ST3784 (subsp. ozaenae)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST3802'), 'ST3802 (subsp. ozaenae)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST3803'), 'ST3803 (subsp. ozaenae)')

    def test_kp_subspecies_3(self):
        # Test the ozaenae STs with 1LV - should still result in ozaenae.
        self.assertEqual(get_kp_subspecies_based_on_st('ST90-1LV'), 'ST90-1LV (subsp. ozaenae)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST91-1LV'), 'ST91-1LV (subsp. ozaenae)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST92-1LV'), 'ST92-1LV (subsp. ozaenae)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST93-1LV'), 'ST93-1LV (subsp. ozaenae)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST95-1LV'), 'ST95-1LV (subsp. ozaenae)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST96-1LV'), 'ST96-1LV (subsp. ozaenae)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST97-1LV'), 'ST97-1LV (subsp. ozaenae)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST381-1LV'), 'ST381-1LV (subsp. ozaenae)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST777-1LV'), 'ST777-1LV (subsp. ozaenae)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST3193-1LV'), 'ST3193-1LV (subsp. ozaenae)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST3766-1LV'), 'ST3766-1LV (subsp. ozaenae)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST3768-1LV'), 'ST3768-1LV (subsp. ozaenae)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST3771-1LV'), 'ST3771-1LV (subsp. ozaenae)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST3781-1LV'), 'ST3781-1LV (subsp. ozaenae)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST3782-1LV'), 'ST3782-1LV (subsp. ozaenae)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST3784-1LV'), 'ST3784-1LV (subsp. ozaenae)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST3802-1LV'), 'ST3802-1LV (subsp. ozaenae)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST3803-1LV'), 'ST3803-1LV (subsp. ozaenae)')

    def test_kp_subspecies_4(self):
        # Test the ozaenae STs with 2LV - should no longer result in ozaenae.
        self.assertEqual(get_kp_subspecies_based_on_st('ST90-2LV'), 'ST90-2LV')
        self.assertEqual(get_kp_subspecies_based_on_st('ST91-2LV'), 'ST91-2LV')
        self.assertEqual(get_kp_subspecies_based_on_st('ST92-2LV'), 'ST92-2LV')
        self.assertEqual(get_kp_subspecies_based_on_st('ST93-2LV'), 'ST93-2LV')
        self.assertEqual(get_kp_subspecies_based_on_st('ST95-2LV'), 'ST95-2LV')
        self.assertEqual(get_kp_subspecies_based_on_st('ST96-2LV'), 'ST96-2LV')
        self.assertEqual(get_kp_subspecies_based_on_st('ST97-2LV'), 'ST97-2LV')
        self.assertEqual(get_kp_subspecies_based_on_st('ST381-2LV'), 'ST381-2LV')
        self.assertEqual(get_kp_subspecies_based_on_st('ST777-2LV'), 'ST777-2LV')
        self.assertEqual(get_kp_subspecies_based_on_st('ST3193-2LV'), 'ST3193-2LV')
        self.assertEqual(get_kp_subspecies_based_on_st('ST3766-2LV'), 'ST3766-2LV')
        self.assertEqual(get_kp_subspecies_based_on_st('ST3768-2LV'), 'ST3768-2LV')
        self.assertEqual(get_kp_subspecies_based_on_st('ST3771-2LV'), 'ST3771-2LV')
        self.assertEqual(get_kp_subspecies_based_on_st('ST3781-2LV'), 'ST3781-2LV')
        self.assertEqual(get_kp_subspecies_based_on_st('ST3782-2LV'), 'ST3782-2LV')
        self.assertEqual(get_kp_subspecies_based_on_st('ST3784-2LV'), 'ST3784-2LV')
        self.assertEqual(get_kp_subspecies_based_on_st('ST3802-2LV'), 'ST3802-2LV')
        self.assertEqual(get_kp_subspecies_based_on_st('ST3803-2LV'), 'ST3803-2LV')

    def test_kp_subspecies_5(self):
        # Test the rhinoscleromatis STs.
        self.assertEqual(get_kp_subspecies_based_on_st('ST67'), 'ST67 (subsp. rhinoscleromatis)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST68'), 'ST68 (subsp. rhinoscleromatis)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST69'), 'ST69 (subsp. rhinoscleromatis)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST3772'),
                         'ST3772 (subsp. rhinoscleromatis)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST3819'),
                         'ST3819 (subsp. rhinoscleromatis)')

    def test_kp_subspecies_6(self):
        # Test the rhinoscleromatis STs with 1LV - should still result in rhinoscleromatis.
        self.assertEqual(get_kp_subspecies_based_on_st('ST67-1LV'),
                         'ST67-1LV (subsp. rhinoscleromatis)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST68-1LV'),
                         'ST68-1LV (subsp. rhinoscleromatis)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST69-1LV'),
                         'ST69-1LV (subsp. rhinoscleromatis)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST3772-1LV'),
                         'ST3772-1LV (subsp. rhinoscleromatis)')
        self.assertEqual(get_kp_subspecies_based_on_st('ST3819-1LV'),
                         'ST3819-1LV (subsp. rhinoscleromatis)')

    def test_kp_subspecies_7(self):
        # Test the rhinoscleromatis STs with 2LV - should no longer result in rhinoscleromatis.
        self.assertEqual(get_kp_subspecies_based_on_st('ST67-2LV'), 'ST67-2LV')
        self.assertEqual(get_kp_subspecies_based_on_st('ST68-2LV'), 'ST68-2LV')
        self.assertEqual(get_kp_subspecies_based_on_st('ST69-2LV'), 'ST69-2LV')
        self.assertEqual(get_kp_subspecies_based_on_st('ST3772-2LV'), 'ST3772-2LV')
        self.assertEqual(get_kp_subspecies_based_on_st('ST3819-2LV'), 'ST3819-2LV')
