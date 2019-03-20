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


class TestResAlleles(unittest.TestCase):
    """
    Tests calling of resistance via alleles.
    """

    def setUp(self):
        self.data_dir = get_data_path()

        Args = collections.namedtuple('Args', ['species', 'resistance', 'kaptive_k', 'kaptive_o'])
        self.args = Args(species=False, resistance=True, kaptive_k=False, kaptive_o=False)
        _, _, self.res_headers = get_output_headers(self.args, self.data_dir)
