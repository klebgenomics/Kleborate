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


class TestResTem(unittest.TestCase):
    """
    This test is for a particular bug we found, where BLAST can find an exact amino acid match but
    on the wrong strand. This test sequence was coming up as "TEM-15^" from a wrong strand match
    until we fixed the bug (only checking the forward strand).
    """
    def setUp(self):
        self.data_dir = 'test/test_res_tem/data'
        Args = collections.namedtuple('Args', ['resistance', 'kaptive_k', 'kaptive_o',
                                               'min_coverage', 'min_identity',
                                               'min_spurious_coverage', 'min_spurious_identity'])
        self.args = Args(resistance=True, kaptive_k=False, kaptive_o=False,
                         min_coverage=80.0, min_identity=90.0,
                         min_spurious_coverage=40.0, min_spurious_identity=80.0)
        _, _, self.res_headers = get_output_headers(self.args, self.data_dir)

    def test_tem(self):
        results = get_resistance_results(self.data_dir, 'test/test_res_tem/tem.fasta',
                                         self.args, self.res_headers, True)
        self.assertEqual(results['Bla_ESBL'], '-')
        self.assertEqual(results['Bla'], 'TEM-1D^')
