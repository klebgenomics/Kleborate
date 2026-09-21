"""
This file contains tests for Kleborate. To run all tests, go the repo's root directory and run:
  python3 -m pytest

To get code coverage stats:
  coverage run --source . -m pytest && coverage report -m

Copyright 2026 Kat Holt, Ryan Wick (rrwick@gmail.com)
https://github.com/katholt/Kleborate/

This file is part of Kleborate. Kleborate is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Kleborate is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Kleborate. If
not, see <https://www.gnu.org/licenses/>.
"""

from .kpsc__virulence_score import *


def test_prerequisite_modules():
    assert sorted(prerequisite_modules()) == ['klebsiella__abst', 'klebsiella__cbst', 'klebsiella__ybst']


def test_empty_functions():
    # Tests the functions that aren't used in this module.
    assert add_cli_options(None) is None
    assert check_cli_options(None) is None
    assert check_external_programs() == []


def test_get_results_1():
    previous_results = {
        'klebsiella__ybst__spurious_YbST': '-',
        'klebsiella__abst__spurious_AbST': '-',
        'klebsiella__cbst__spurious_CbST': '-',
        'klebsiella__rmst__spurious_RmST': '-',
        'klebsiella__smst__spurious_SmST': '-',
        'klebsiella__ybst__YbST': 0,
        'klebsiella__abst__AbST': 0,
        'klebsiella__cbst__CbST': 0,
        'klebsiella__smst__SmST': 0,
        'klebsiella__rmst__RmST': 0,
        'klebsiella__ybst__Yersiniabactin': '-',
        'klebsiella__abst__Aerobactin': '-',
        'klebsiella__cbst__Colibactin': '-',
        'klebsiella__smst__Salmochelin': '-',
        'klebsiella__rmst__RmpADC': '-',
        'klebsiella__rmst__RmpADC_status': '-'
    }
    result = get_results(None, None, None, previous_results)
    assert result['virulence_score'] == '0'


def test_get_results_2():
    previous_results = {
        'klebsiella__ybst__spurious_YbST': '-',
        'klebsiella__abst__spurious_AbST': '-',
        'klebsiella__cbst__spurious_CbST': '-',
        'klebsiella__rmst__spurious_RmST': '-',
        'klebsiella__smst__spurious_SmST': '-',
        'klebsiella__ybst__YbST': 1,
        'klebsiella__abst__AbST': 0,
        'klebsiella__cbst__CbST': 0,
        'klebsiella__smst__SmST': 0,
        'klebsiella__rmst__RmST': 0,
        'klebsiella__ybst__Yersiniabactin': 'ybt 1',
        'klebsiella__abst__Aerobactin': '-',
        'klebsiella__cbst__Colibactin': '-',
        'klebsiella__smst__Salmochelin': '-',
        'klebsiella__rmst__RmpADC': '-',
        'klebsiella__rmst__RmpADC_status': '-'
    }
    result = get_results(None, None, None, previous_results)
    assert result['virulence_score'] == '1'


def test_get_results_3():
    previous_results = {
        'klebsiella__ybst__spurious_YbST': '-',
        'klebsiella__abst__spurious_AbST': '-',
        'klebsiella__cbst__spurious_CbST': '-',
        'klebsiella__rmst__spurious_RmST': '-',
        'klebsiella__smst__spurious_SmST': '-',
        'klebsiella__ybst__YbST': 0,
        'klebsiella__abst__AbST': 0,
        'klebsiella__cbst__CbST': 1,
        'klebsiella__smst__SmST': 0,
        'klebsiella__rmst__RmST': 0,
        'klebsiella__ybst__Yersiniabactin': '-',
        'klebsiella__abst__Aerobactin': '-',
        'klebsiella__cbst__Colibactin': 'clb 1',
        'klebsiella__smst__Salmochelin': '-',
        'klebsiella__rmst__RmpADC': '-',
        'klebsiella__rmst__RmpADC_status': '-'
    }
    result = get_results(None, None, None, previous_results)
    assert result['virulence_score'] == '2'


def test_get_results_4():
    previous_results = {
        'klebsiella__ybst__spurious_YbST': '-',
        'klebsiella__abst__spurious_AbST': '-',
        'klebsiella__cbst__spurious_CbST': '-',
        'klebsiella__rmst__spurious_RmST': '-',
        'klebsiella__smst__spurious_SmST': '-',
        'klebsiella__ybst__YbST': 1,
        'klebsiella__abst__AbST': 0,
        'klebsiella__cbst__CbST': 1,
        'klebsiella__smst__SmST': 0,
        'klebsiella__rmst__RmST': 0,
        'klebsiella__ybst__Yersiniabactin': 'ybt 1',
        'klebsiella__abst__Aerobactin': '-',
        'klebsiella__cbst__Colibactin': 'clb 1',
        'klebsiella__smst__Salmochelin': '-',
        'klebsiella__rmst__RmpADC': '-',
        'klebsiella__rmst__RmpADC_status': '-'
    }
    result = get_results(None, None, None, previous_results)
    assert result['virulence_score'] == '2'


def test_get_results_5():
    previous_results = {
        'klebsiella__ybst__spurious_YbST': '-',
        'klebsiella__abst__spurious_AbST': '-',
        'klebsiella__cbst__spurious_CbST': '-',
        'klebsiella__rmst__spurious_RmST': '-',
        'klebsiella__smst__spurious_SmST': '-',
        'klebsiella__ybst__YbST': 0,
        'klebsiella__abst__AbST': 1,
        'klebsiella__cbst__CbST': 0,
        'klebsiella__smst__SmST': 0,
        'klebsiella__rmst__RmST': 0,
        'klebsiella__ybst__Yersiniabactin': '-',
        'klebsiella__abst__Aerobactin': 'iuc 1',
        'klebsiella__cbst__Colibactin': '-',
        'klebsiella__smst__Salmochelin': '-',
        'klebsiella__rmst__RmpADC': '-',
        'klebsiella__rmst__RmpADC_status': '-'
    }
    result = get_results(None, None, None, previous_results)
    assert result['virulence_score'] == '3'


def test_get_results_6():
    previous_results = {
        'klebsiella__ybst__spurious_YbST': '-',
        'klebsiella__abst__spurious_AbST': '-',
        'klebsiella__cbst__spurious_CbST': '-',
        'klebsiella__rmst__spurious_RmST': '-',
        'klebsiella__smst__spurious_SmST': '-',
        'klebsiella__ybst__YbST': 1,
        'klebsiella__abst__AbST': 1,
        'klebsiella__cbst__CbST': 0,
        'klebsiella__smst__SmST': 0,
        'klebsiella__rmst__RmST': 0,
        'klebsiella__ybst__Yersiniabactin': 'ybt 1',
        'klebsiella__abst__Aerobactin': 'iuc 1',
        'klebsiella__cbst__Colibactin': '-',
        'klebsiella__smst__Salmochelin': '-',
        'klebsiella__rmst__RmpADC': '-',
        'klebsiella__rmst__RmpADC_status': '-'
    }
    result = get_results(None, None, None, previous_results)
    assert result['virulence_score'] == '4'


def test_get_results_7():
    previous_results = {
        'klebsiella__ybst__spurious_YbST': '-',
        'klebsiella__abst__spurious_AbST': '-',
        'klebsiella__cbst__spurious_CbST': '-',
        'klebsiella__rmst__spurious_RmST': '-',
        'klebsiella__smst__spurious_SmST': '-',
        'klebsiella__ybst__YbST': 0,
        'klebsiella__abst__AbST': 1,
        'klebsiella__cbst__CbST': 1,
        'klebsiella__smst__SmST': 0,
        'klebsiella__rmst__RmST': 0,
        'klebsiella__ybst__Yersiniabactin': '-',
        'klebsiella__abst__Aerobactin': 'iuc 1',
        'klebsiella__cbst__Colibactin': 'clb 1',
        'klebsiella__smst__Salmochelin': '-',
        'klebsiella__rmst__RmpADC': '-',
        'klebsiella__rmst__RmpADC_status': '-'
    }
    result = get_results(None, None, None, previous_results)
    assert result['virulence_score'] == '5'


def test_get_results_8():
    previous_results = {
        'klebsiella__ybst__spurious_YbST': '-',
        'klebsiella__abst__spurious_AbST': '-',
        'klebsiella__cbst__spurious_CbST': '-',
        'klebsiella__rmst__spurious_RmST': '-',
        'klebsiella__smst__spurious_SmST': '-',
        'klebsiella__ybst__YbST': 1,
        'klebsiella__abst__AbST': 1,
        'klebsiella__cbst__CbST': 1,
        'klebsiella__smst__SmST': 0,
        'klebsiella__rmst__RmST': 0,
        'klebsiella__ybst__Yersiniabactin': 'ybt 1',
        'klebsiella__abst__Aerobactin': 'iuc 1',
        'klebsiella__cbst__Colibactin': 'clb 1',
        'klebsiella__smst__Salmochelin': '-',
        'klebsiella__rmst__RmpADC': '-',
        'klebsiella__rmst__RmpADC_status': '-'
    }
    result = get_results(None, None, None, previous_results)
    assert result['virulence_score'] == '5'