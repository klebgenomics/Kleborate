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

from kleborate.__main__ import get_output_headers, get_resistance_results
from kleborate.resBLAST import is_exact_aa_match


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
        self.assertEqual(results['Tet_acquired'], '-')
        self.assertEqual(results['Bla_acquired'], 'ABC-1')
        self.assertEqual(results['Bla_ESBL_acquired'], '-')
        self.assertEqual(results['spurious_resistance_hits'], '-')

    def test_res_02(self):
        results = get_resistance_results(self.data_dir, 'test/test_res_alleles/02.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Tet_acquired'], '-')
        self.assertEqual(results['Bla_acquired'], 'ABC-2')
        self.assertEqual(results['Bla_ESBL_acquired'], '-')
        self.assertEqual(results['spurious_resistance_hits'], '-')

    def test_res_03(self):
        results = get_resistance_results(self.data_dir, 'test/test_res_alleles/03.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Tet_acquired'], '-')
        self.assertEqual(results['Bla_acquired'], '-')
        self.assertEqual(results['Bla_ESBL_acquired'], 'ABC-3')
        self.assertEqual(results['spurious_resistance_hits'], '-')

    def test_res_04(self):
        results = get_resistance_results(self.data_dir, 'test/test_res_alleles/04.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Tet_acquired'], '-')
        self.assertEqual(results['Bla_acquired'], '-')
        self.assertEqual(results['Bla_ESBL_acquired'], 'ABC-4')
        self.assertEqual(results['spurious_resistance_hits'], '-')

    def test_res_05(self):
        results = get_resistance_results(self.data_dir, 'test/test_res_alleles/05.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Tet_acquired'], '-')
        self.assertTrue(results['Bla_acquired'] == 'ABC-1^' or results['Bla_acquired'] == 'ABC-2^')
        self.assertEqual(results['Bla_ESBL_acquired'], '-')
        self.assertEqual(results['spurious_resistance_hits'], '-')

    def test_res_06(self):
        results = get_resistance_results(self.data_dir, 'test/test_res_alleles/06.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Tet_acquired'], '-')
        self.assertEqual(results['Bla_acquired'], '-')
        self.assertTrue(results['Bla_ESBL_acquired'] == 'ABC-3^' or
                        results['Bla_ESBL_acquired'] == 'ABC-4^')
        self.assertEqual(results['spurious_resistance_hits'], '-')

    def test_res_07(self):
        results = get_resistance_results(self.data_dir, 'test/test_res_alleles/07.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Tet_acquired'], '-')
        self.assertEqual(results['Bla_acquired'], 'ABC-2*')
        self.assertEqual(results['Bla_ESBL_acquired'], '-')
        self.assertEqual(results['spurious_resistance_hits'], '-')

    def test_res_08(self):
        results = get_resistance_results(self.data_dir, 'test/test_res_alleles/08.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Tet_acquired'], '-')
        self.assertEqual(results['Bla_acquired'], 'ABC-1;ABC-2')
        self.assertEqual(results['Bla_ESBL_acquired'], '-')
        self.assertEqual(results['spurious_resistance_hits'], '-')

    def test_res_09(self):
        results = get_resistance_results(self.data_dir, 'test/test_res_alleles/09.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Tet_acquired'], '-')
        self.assertTrue(results['Bla_acquired'] == 'ABC-1^' or results['Bla_acquired'] == 'ABC-2^')
        self.assertTrue(results['Bla_ESBL_acquired'] == 'ABC-3^' or
                        results['Bla_ESBL_acquired'] == 'ABC-4^')
        self.assertEqual(results['spurious_resistance_hits'], '-')

    def test_res_10(self):
        results = get_resistance_results(self.data_dir, 'test/test_res_alleles/10.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Tet_acquired'], '-')
        self.assertEqual(results['Bla_acquired'], 'ABC-1;ABC-2')
        self.assertEqual(results['Bla_ESBL_acquired'], 'ABC-3*;ABC-4')
        self.assertEqual(results['spurious_resistance_hits'], '-')

    def test_res_11(self):
        """
        This test has the ABC-1 gene but with a stop codon in the middle.
        """
        results = get_resistance_results(self.data_dir, 'test/test_res_alleles/11.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Tet_acquired'], '-')
        self.assertEqual(results['Bla_acquired'], '-')
        self.assertEqual(results['Bla_ESBL_acquired'], '-')
        self.assertEqual(results['truncated_resistance_hits'], 'ABC-1*-50%')

    def test_res_12(self):
        """
        This test has the ABC-1 gene but missing the first (start) codon.
        """
        results = get_resistance_results(self.data_dir, 'test/test_res_alleles/12.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Tet_acquired'], '-')
        self.assertEqual(results['Bla_acquired'], '-')
        self.assertEqual(results['Bla_ESBL_acquired'], '-')
        self.assertEqual(results['truncated_resistance_hits'], 'ABC-1?-0%')

    def test_res_13(self):
        """
        This test has the ABC-1 gene but at only 50% coverage, so it goes to spurious hits.
        """
        results = get_resistance_results(self.data_dir, 'test/test_res_alleles/13.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Tet_acquired'], '-')
        self.assertEqual(results['Bla_acquired'], '-')
        self.assertEqual(results['Bla_ESBL_acquired'], '-')
        self.assertEqual(results['spurious_resistance_hits'], 'ABC-1?-50%')

    def test_res_14(self):
        """
        I added this test to catch a weird BLAST-related bug where an exact amino acid match for
        SHV-51 was present but not being found.
        """
        results = get_resistance_results(self.data_dir, 'test/test_res_alleles/14.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Bla_chr'], 'SHV-51^')
        self.assertEqual(results['Bla_acquired'], '-')
        self.assertEqual(results['Bla_ESBL_acquired'], '-')
        self.assertEqual(results['spurious_resistance_hits'], '-')

    def test_res_15(self):
        """
        Same as the previous test, but with the nucleotide sequence reverse-complemented.
        """
        results = get_resistance_results(self.data_dir, 'test/test_res_alleles/15.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Bla_chr'], 'SHV-51^')
        self.assertEqual(results['Bla_acquired'], '-')
        self.assertEqual(results['Bla_ESBL_acquired'], '-')
        self.assertEqual(results['spurious_resistance_hits'], '-')

    def test_res_16(self):
        """
        I added this test to catch a weird BLAST-related bug where an exact amino acid match for
        sul2 was present but not being found.
        """
        results = get_resistance_results(self.data_dir, 'test/test_res_alleles/16.fasta', self.args,
                                         self.res_headers, True)
        print(results)
        self.assertEqual(results['Sul_acquired'], 'sul2^')

    def test_res_17(self):
        """
        Same as the previous test, but with the nucleotide sequence reverse-complemented.
        """
        results = get_resistance_results(self.data_dir, 'test/test_res_alleles/17.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Sul_acquired'], 'sul2^')

    def test_res_18(self):
        """
        I added this test to catch a bug where there is an exact amino acid match but inexact
        nucleotide match, and it's the very last base that differs (alternate stop codon).
        """
        results = get_resistance_results(self.data_dir, 'test/test_res_alleles/18.fasta', self.args,
                                         self.res_headers, True)
        print(results)
        self.assertEqual(results['Bla_chr'], 'OKP-B-6^')

    def test_res_19(self):
        """
        Same as the previous test, but with the nucleotide sequence reverse-complemented.
        """
        results = get_resistance_results(self.data_dir, 'test/test_res_alleles/19.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Bla_chr'], 'OKP-B-6^')

    def test_exact_aa_match_1(self):
        """
        Simple case: full protein sequence ending in a stop codon.
          seq_1 -> MNKSLIIFGI*
          seq_2 -> MNKSLIIFGI*
        """
        seq_1 = 'ATGAATAAATCGCTAATCATTTTCGGCATCTAA'
        seq_2 = 'ATGAATAAATCGCTAATCATTTTCGGCATCTAA'
        self.assertTrue(is_exact_aa_match(seq_1, seq_2))

    def test_exact_aa_match_2(self):
        """
        These two seqs have a stop codon in the middle but still match exactly.
          seq_1 -> MNKSLIIF*GI
          seq_2 -> MNKSLIIF*GI
        """
        seq_1 = 'ATGAATAAATCGCTAATCATTTTCTAAGGCATC'
        seq_2 = 'ATGAACAAGTCGCTCATCATTTTCTAAGGCATC'
        self.assertTrue(is_exact_aa_match(seq_1, seq_2))

    def test_exact_aa_match_3(self):
        """
        These two seqs have a stop codon in the middle and mismatch after the stop.
          seq_1 -> MNKSLIIF*GI
          seq_2 -> MNKSLIIF*TL
        """
        seq_1 = 'ATGAATAAATCGCTAATCATTTTCTAAGGCATC'
        seq_2 = 'ATGAACAAGTCGCTCATCATTTTCTAAACCTTA'
        self.assertFalse(is_exact_aa_match(seq_1, seq_2))

    def test_exact_aa_match_4(self):
        """
        These two seqs have a stop codon in the middle and mismatch after the stop.
          seq_1 -> MNK*NYRTDYD
          seq_2 -> MNK*VD*LAGP
        """
        seq_1 = 'ATGAATAAATAAAACTATCGTACGGACTACGAC'
        seq_2 = 'ATGAACAAGTAAGTCGACTGACTAGCTGGACCC'
        self.assertFalse(is_exact_aa_match(seq_1, seq_2))

    def test_exact_aa_match_5(self):
        """
        These two seqs are non-multiple-of-three in length and match in their remainder.
          seq_1 -> MNKSLIIFG (+AT)
          seq_2 -> MNKSLIIFG (+AT)
        """
        seq_1 = 'ATGAATAAATCGCTAATCATTTTCGGCAT'
        seq_2 = 'ATGAACAAGTCGCTCATCATTTTCGGGAT'
        self.assertTrue(is_exact_aa_match(seq_1, seq_2))

    def test_exact_aa_match_6(self):
        """
        These two seqs are non-multiple-of-three in length and match in their remainder.
          seq_1 -> MNKSLIIFG (+AT)
          seq_2 -> MNKSLIIFG (+A)
        """
        seq_1 = 'ATGAATAAATCGCTAATCATTTTCGGCAT'
        seq_2 = 'ATGAACAAGTCGCTCATCATTTTCGGGA'
        self.assertTrue(is_exact_aa_match(seq_1, seq_2))

    def test_exact_aa_match_7(self):
        """
        These two seqs are non-multiple-of-three in length and don't match in their remainder.
          seq_1 -> MNKSLIIFG (+AT)
          seq_2 -> MNKSLIIFG (+CC)
        """
        seq_1 = 'ATGAATAAATCGCTAATCATTTTCGGCAT'
        seq_2 = 'ATGAACAAGTCGCTCATCATTTTCGGGCC'
        self.assertTrue(is_exact_aa_match(seq_1, seq_2))

    def test_exact_aa_match_8(self):
        """
        The first seq's protein sequence contains the second seq's protein sequence.
          seq_1 -> DTMNKSLGI*RS
          seq_2 -> MNKSLGI*
        """
        seq_1 = 'GACACGATGAATAAATCGCTAGGCATCTAACGATCA'
        seq_2 = 'ATGAACAAGTCGCTCGGGATCTAA'
        self.assertTrue(is_exact_aa_match(seq_1, seq_2))
