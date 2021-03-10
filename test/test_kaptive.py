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

import unittest

from kleborate.kaptive import get_o_type, confidence_meets_threshold


class TestKaptive(unittest.TestCase):

    def test_o_type_01(self):
        self.assertEqual(get_o_type('O1v1'), 'O1')

    def test_o_type_02(self):
        self.assertEqual(get_o_type('O1v2'), 'O1')

    def test_o_type_03(self):
        self.assertEqual(get_o_type('O2v1'), 'O2')

    def test_o_type_04(self):
        self.assertEqual(get_o_type('O2v2'), 'O2')

    def test_o_type_05(self):
        self.assertEqual(get_o_type('O1/O2v1'), 'unknown')

    def test_o_type_06(self):
        self.assertEqual(get_o_type('O1/O2v2'), 'unknown')

    def test_o_type_07(self):
        self.assertEqual(get_o_type('O3/O3a'), 'O3/O3a')

    def test_o_type_08(self):
        self.assertEqual(get_o_type('O3b'), 'O3b')

    def test_o_type_09(self):
        self.assertEqual(get_o_type('O4'), 'O4')

    def test_o_type_10(self):
        self.assertEqual(get_o_type('O5'), 'O5')

    def test_o_type_11(self):
        self.assertEqual(get_o_type('O8'), 'O8')

    def test_o_type_12(self):
        self.assertEqual(get_o_type('O12'), 'O12')

    def test_o_type_13(self):
        self.assertEqual(get_o_type('OL101'), 'OL101')

    def test_o_type_14(self):
        self.assertEqual(get_o_type('OL102'), 'OL102')

    def test_o_type_15(self):
        self.assertEqual(get_o_type('OL103'), 'OL103')

    def test_o_type_16(self):
        self.assertEqual(get_o_type('OL104'), 'OL104')

    def test_confidence_meets_threshold_01(self):
        self.assertTrue(confidence_meets_threshold('None', 'None'))

    def test_confidence_meets_threshold_02(self):
        self.assertTrue(confidence_meets_threshold('Good', 'Good'))

    def test_confidence_meets_threshold_03(self):
        self.assertTrue(confidence_meets_threshold('Very high', 'Very_high'))

    def test_confidence_meets_threshold_04(self):
        self.assertTrue(confidence_meets_threshold('Perfect', 'Perfect'))

    def test_confidence_meets_threshold_05(self):
        self.assertTrue(confidence_meets_threshold('Perfect', 'Good'))

    def test_confidence_meets_threshold_06(self):
        self.assertTrue(confidence_meets_threshold('Perfect', 'Very_high'))

    def test_confidence_meets_threshold_07(self):
        self.assertFalse(confidence_meets_threshold('Very high', 'Perfect'))

    def test_confidence_meets_threshold_08(self):
        self.assertFalse(confidence_meets_threshold('High', 'Very_high'))

    def test_confidence_meets_threshold_09(self):
        self.assertFalse(confidence_meets_threshold('Good', 'High'))

    def test_confidence_meets_threshold_10(self):
        self.assertFalse(confidence_meets_threshold('Low', 'Good'))
