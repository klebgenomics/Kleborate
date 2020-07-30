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


class TestResAlleles(unittest.TestCase):
    """
    Tests calling of resistance via alleles.
    """

    def setUp(self):
        self.data_dir = 'test/test_shv/data'
        Args = collections.namedtuple('Args', ['resistance', 'kaptive_k', 'kaptive_o',
                                               'min_coverage', 'min_identity',
                                               'min_spurious_coverage', 'min_spurious_identity'])
        self.args = Args(resistance=True, kaptive_k=False, kaptive_o=False,
                         min_coverage=80.0, min_identity=90.0,
                         min_spurious_coverage=40.0, min_spurious_identity=80.0)
        _, _, self.res_headers = get_output_headers(self.args, self.data_dir)

    def test_shv_01(self):
        """
        This test has an exact match for SHV-1.
        """
        results = get_resistance_results(self.data_dir, 'test/test_shv/01.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Bla_chr'], 'SHV-1')
        self.assertEqual(results['SHV_mutations'], '')

    def test_shv_02(self):
        """
        This test has a match for SHV-1 with a mutation at site 238 (G -> Y). This changes the
        class to ESBL, so the mutation is included in the
        """
        results = get_resistance_results(self.data_dir, 'test/test_shv/02.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Bla_ESBL_acquired'], 'SHV-1*-238Y')
        self.assertEqual(results['SHV_mutations'], '238Y')

    def test_shv_03(self):
        """
        Same as test 2, but the gene is on the reverse strand.
        """
        results = get_resistance_results(self.data_dir, 'test/test_shv/03.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Bla_ESBL_acquired'], 'SHV-1*-238Y')
        self.assertEqual(results['SHV_mutations'], '238Y')

    def test_shv_04(self):
        """
        This test has a match for SHV-1 with a mutation at site 50 (G -> Y). This doesn't change
        resistance and so won't be reported.
        """
        results = get_resistance_results(self.data_dir, 'test/test_shv/04.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Bla_chr'], 'SHV-1*')
        self.assertEqual(results['SHV_mutations'], '')

    def test_shv_05(self):
        """
        This test has an exact match for SHV-29.
        """
        results = get_resistance_results(self.data_dir, 'test/test_shv/05.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Bla_ESBL_acquired'], 'SHV-29')
        self.assertEqual(results['SHV_mutations'], '238A;35Q')

    def test_shv_06(self):
        """
        This test has SHV-29 plus an inhibition mutation.
        """
        results = get_resistance_results(self.data_dir, 'test/test_shv/06.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Bla_ESBL_inhR_acquired'], 'SHV-29*-234Y')
        self.assertEqual(results['SHV_mutations'], '234Y;238A;35Q')

    def test_shv_07(self):
        """
        This test has SHV-1 with position 238 deleted. Since it's not in the omega loop, this isn't
        reported and doesn't have an effect.
        """
        results = get_resistance_results(self.data_dir, 'test/test_shv/07.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Bla_chr'], 'SHV-1*')
        self.assertEqual(results['SHV_mutations'], '')

    def test_shv_08(self):
        """
        This test has SHV-1 with a synonymous mutation in the omega loop (so not reported).
        """
        results = get_resistance_results(self.data_dir, 'test/test_shv/08.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Bla_chr'], 'SHV-1^')
        self.assertEqual(results['SHV_mutations'], '')

    def test_shv_09(self):
        """
        This test has SHV-1 with a nonsynonymous mutation in the omega loop (so it is reported).
        """
        results = get_resistance_results(self.data_dir, 'test/test_shv/09.fasta', self.args,
                                         self.res_headers, True)
        self.assertEqual(results['Bla_chr'], 'SHV-1*')
        self.assertEqual(results['SHV_mutations'], 'omega-loop=RWETELNEALRGDARD')

