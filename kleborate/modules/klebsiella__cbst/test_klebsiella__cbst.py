"""
This file contains tests for Kleborate. To run all tests, go the repo's root directory and run:
  python3 -m pytest

To get code coverage stats:
  coverage run --source . -m pytest && coverage report -m

Copyright 2023 Kat Holt, Ryan Wick (rrwick@gmail.com), Mary Maranga (gathonimaranga@gmail.com)
https://github.com/katholt/KleborateModular/

This file is part of Kleborate. Kleborate is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Kleborate is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Kleborate. If
not, see <https://www.gnu.org/licenses/>.
"""

import collections
import pytest

from .klebsiella__cbst import *


def get_test_genome_dir():
    return pathlib.Path(__file__).parents[3] / 'test' / 'test_genomes'


def test_prerequisite_modules():
    assert prerequisite_modules() == []


def test_check_cli_options_1():
    Args = collections.namedtuple('Args', ['klebsiella__cbst_min_identity', 'klebsiella__cbst_min_coverage',
                                           'klebsiella__cbst_min_spurious_identity', 'klebsiella__cbst_min_spurious_coverage',
                                           'klebsiella__cbst_required_exact_matches'])
    check_cli_options(Args(klebsiella__cbst_min_identity=90.0, klebsiella__cbst_min_coverage=90.0,
                           klebsiella__cbst_min_spurious_identity=80.0, klebsiella__cbst_min_spurious_coverage=40.0,
                           klebsiella__cbst_required_exact_matches=3))



def test_check_cli_options_2():
    Args = collections.namedtuple('Args', ['klebsiella__cbst_min_identity', 'klebsiella__cbst_min_coverage',
                                           'klebsiella__cbst_required_exact_matches'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(klebsiella__cbst_min_identity=0.90, klebsiella__cbst_min_coverage=90.0,
                               klebsiella__cbst_required_exact_matches=3))


def test_check_cli_options_3():
    Args = collections.namedtuple('Args', ['klebsiella__cbst_min_identity', 'klebsiella__cbst_min_coverage',
                                           'klebsiella__cbst_required_exact_matches'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(klebsiella__cbst_min_identity=-90.0, klebsiella__cbst_min_coverage=0.90,
                               klebsiella__cbst_required_exact_matches=3))


def test_check_cli_options_4():
    Args = collections.namedtuple('Args', ['klebsiella__cbst_min_identity', 'klebsiella__cbst_min_coverage',
                                           'klebsiella__cbst_required_exact_matches'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(klebsiella__cbst_min_identity=-10.0, klebsiella__cbst_min_coverage=90.0,
                               klebsiella__cbst_required_exact_matches=3))


def test_check_cli_options_5():
    Args = collections.namedtuple('Args', ['klebsiella__cbst_min_identity', 'klebsiella__cbst_min_coverage',
                                           'klebsiella__cbst_required_exact_matches'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(klebsiella__cbst_min_identity=90.0, klebsiella__cbst_min_coverage=120.0,
                               klebsiella__cbst_required_exact_matches=3))


def test_check_cli_options_6():
    Args = collections.namedtuple('Args', ['klebsiella__cbst_min_identity', 'klebsiella__cbst_min_coverage',
                                           'klebsiella__cbst_min_spurious_identity', 'klebsiella__cbst_min_spurious_coverage',
                                           'klebsiella__cbst_required_exact_matches'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(klebsiella__cbst_min_identity=90.0, klebsiella__cbst_min_coverage=90.0,
                               klebsiella__cbst_min_spurious_identity=80.0, klebsiella__cbst_min_spurious_coverage=40.0,
                               klebsiella__cbst_required_exact_matches=-2))


def test_check_external_programs_1(mocker):
    # Tests the good case where minimap2 is found.
    mocker.patch(
        'shutil.which',
        side_effect=lambda x: {'minimap2': '/usr/bin/minimap2'}[x],
    )
    assert check_external_programs() == ['minimap2']


def test_check_external_programs_2(mocker):
    # Tests the bad case where minimap2 is missing.
    mocker.patch(
        'shutil.which',
        side_effect=lambda x: {'minimap2': None}[x],
    )
    with pytest.raises(SystemExit):
        check_external_programs()


def test_get_results_1():
    Args = collections.namedtuple('Args', ['klebsiella__cbst_min_identity', 'klebsiella__cbst_min_coverage',
                                           'klebsiella__cbst_min_spurious_identity', 'klebsiella__cbst_min_spurious_coverage',
                                           'klebsiella__cbst_required_exact_matches'])
    results = get_results(get_test_genome_dir() / 'GCF_000968155.1.fna.gz', None,
                          Args(klebsiella__cbst_min_identity=90.0, klebsiella__cbst_min_coverage=80.0,
                               klebsiella__cbst_min_spurious_identity=80.0, klebsiella__cbst_min_spurious_coverage=40.0,
                               klebsiella__cbst_required_exact_matches=3), {})

    assert results['CbST'] == '9'
    assert results['Colibactin'] == 'clb 1'
    assert results['clbA'] == '1'
    assert results['clbB'] == '4'
    assert results['clbC'] == '1'
    assert results['clbD'] == '1'
    assert results['clbE'] == '1'
    assert results['clbF'] == '1'
    assert results['clbG'] == '1'
    assert results['clbH'] == '5'
    assert results['clbI'] == '1'
    assert results['clbL'] == '1'
    assert results['clbM'] == '1'
    assert results['clbN'] == '1'
    assert results['clbO'] == '1'
    assert results['clbP'] == '1'
    assert results['clbQ'] == '1'


def test_get_results_2():
    # Tests an E. coli without the iro locus, so no ST should be assigned.
    Args = collections.namedtuple('Args', ['klebsiella__cbst_min_identity', 'klebsiella__cbst_min_coverage',
                                           'klebsiella__cbst_min_spurious_identity', 'klebsiella__cbst_min_spurious_coverage',
                                           'klebsiella__cbst_required_exact_matches'])
    results = get_results(get_test_genome_dir() / 'GCF_000008865.2.fna.gz', None,
                          Args(klebsiella__cbst_min_identity=90.0, klebsiella__cbst_min_coverage=80.0,
                               klebsiella__cbst_min_spurious_identity=80.0, klebsiella__cbst_min_spurious_coverage=40.0,
                               klebsiella__cbst_required_exact_matches=3), {})

    assert results['CbST'] == 0
    assert results['Colibactin'] == '-'
