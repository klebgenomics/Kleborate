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
        self.kaptive_k = False
        self.kaptive_o = False


class TestResTem(unittest.TestCase):
    """
    This test is for a particular bug we found, where BLAST can find an exact amino acid match but
    on the wrong strand. This test sequence was coming up as "TEM-15^" from a wrong strand match
    until we fixed the bug (only checking the forward strand).
    """
    def setUp(self):
        self.args = Args()
        self.data_dir = get_data_path()
        _, _, self.res_headers = get_output_headers(self.args, self.data_dir)

    def test_tem(self):
        results = get_resistance_results(self.data_dir, 'test/sequences/tem.fasta',
                                         self.args, self.res_headers, True)
        self.assertEqual(results['Bla_ESBL'], '-')
        self.assertEqual(results['Bla'], 'TEM-1D^')
