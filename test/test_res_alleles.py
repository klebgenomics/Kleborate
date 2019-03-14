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
from kleborate import settings


class TestResAlleles(unittest.TestCase):
    """
    Tests calling of resistance via alleles.
    """

    def setUp(self):
        self.data_dir = get_data_path()

        Args = collections.namedtuple('Args', ['species', 'resistance', 'kaptive_k', 'kaptive_o'])
        self.args = Args(species=False, resistance=True, kaptive_k=False, kaptive_o=False)
        _, _, self.res_headers = get_output_headers(self.args, self.data_dir)

        self.inex = settings.inexact_nucleotide_match
        self.inex_aa = settings.inexact_nucleotide_exact_amino_acid

    def test_83(self):
        contigs = 'test/sequences/83.fasta.gz'
        with tempfile.TemporaryDirectory() as tmp_dir:
            contigs = gunzip_contigs_if_necessary(contigs, tmp_dir)
            results = get_resistance_results(self.data_dir, contigs, self.args, self.res_headers)
            for res_class in ['AGly', 'Col', 'Fcyn', 'Flq', 'Gly', 'MLS', 'Ntmdz', 'Phe', 'Rif',
                              'Sul', 'Tet', 'Tmt', 'Bla_Carb', 'Bla_ESBL', 'Bla_ESBL_inhR',
                              'Bla_broad', 'Bla_broad_inhR']:
                self.assertEqual(results[res_class], '-')
            bla_result = 'AmpH{};SHV-119{}'.format(self.inex, self.inex_aa)
            self.assertEqual(results['Bla'], bla_result)

    def test_134(self):
        contigs = 'test/sequences/134.fasta.gz'
        with tempfile.TemporaryDirectory() as tmp_dir:
            contigs = gunzip_contigs_if_necessary(contigs, tmp_dir)
            results = get_resistance_results(self.data_dir, contigs, self.args, self.res_headers)
            for res_class in ['Col', 'Fcyn', 'Gly', 'Ntmdz', 'Bla_Carb', 'Bla_broad_inhR']:
                self.assertEqual(results[res_class], '-')
            agly_result = 'Aac3-IId{};AadA16{};Aph3-Ia;StrA{};StrB'.format(self.inex, self.inex,
                                                                           self.inex, self.inex)
            self.assertEqual(results['AGly'], agly_result)
            self.assertEqual(results['MLS'], 'MphA')
            self.assertEqual(results['Phe'], 'CatA2{}'.format(self.inex))
            self.assertEqual(results['Rif'], 'Arr3')
            self.assertEqual(results['Tet'], 'TetA')
            self.assertEqual(results['Tmt'], 'DfrA27')
            self.assertEqual(results['Bla'], 'AmpH{}'.format(self.inex))
            self.assertEqual(results['Bla_ESBL'], 'CTX-M-15')

            # SulI is in the genome twice and should therefore appear in the results twice.
            self.assertEqual(results['Sul'], 'SulI;SulI;SulII{}'.format(self.inex))

            flq_result = 'GyrA-83F;GyrA-87N;ParC-84K;QnrB17{}'.format(self.inex)
            self.assertEqual(results['Flq'], flq_result)

            bla_broad_result = 'SHV-1{};TEM-135{}'.format(self.inex_aa, self.inex)
            self.assertEqual(results['Bla_broad'], bla_broad_result)
