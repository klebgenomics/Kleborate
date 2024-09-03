"""
This file contains tests for Kleborate. To run all tests, go the repo's root directory and run:
  python3 -m pytest

To get code coverage stats:
  coverage run --source . -m pytest && coverage report -m

Copyright 2023 Kat Holt
Copyright 2023 Ryan Wick (rrwick@gmail.com)
https://github.com/katholt/Kleborate/

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

from .escherichia__mlst_achtman import *


def get_test_genome_dir():
    return pathlib.Path(__file__).parents[3] / 'test' / 'test_genomes'


def test_prerequisite_modules():
    assert prerequisite_modules() == []


def test_check_cli_options_1():
    Args = collections.namedtuple('Args', ['escherichia_mlst_achtman_min_identity',
                                           'escherichia_mlst_achtman_min_coverage',
                                           'escherichia_mlst_achtman_required_exact_matches'])
    check_cli_options(Args(escherichia_mlst_achtman_min_identity=90.0,
                           escherichia_mlst_achtman_min_coverage=90.0,
                           escherichia_mlst_achtman_required_exact_matches=3))


def test_check_cli_options_2():
    Args = collections.namedtuple('Args', ['escherichia_mlst_achtman_min_identity',
                                           'escherichia_mlst_achtman_min_coverage',
                                           'escherichia_mlst_achtman_required_exact_matches'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(escherichia_mlst_achtman_min_identity=0.90,
                               escherichia_mlst_achtman_min_coverage=90.0,
                               escherichia_mlst_achtman_required_exact_matches=3))


def test_check_cli_options_3():
    Args = collections.namedtuple('Args', ['escherichia_mlst_achtman_min_identity',
                                           'escherichia_mlst_achtman_min_coverage',
                                           'escherichia_mlst_achtman_required_exact_matches'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(escherichia_mlst_achtman_min_identity=-90.0,
                               escherichia_mlst_achtman_min_coverage=0.90,
                               escherichia_mlst_achtman_required_exact_matches=3))


def test_check_cli_options_4():
    Args = collections.namedtuple('Args', ['escherichia_mlst_achtman_min_identity',
                                           'escherichia_mlst_achtman_min_coverage',
                                           'escherichia_mlst_achtman_required_exact_matches'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(escherichia_mlst_achtman_min_identity=-10.0,
                               escherichia_mlst_achtman_min_coverage=90.0,
                               escherichia_mlst_achtman_required_exact_matches=3))


def test_check_cli_options_5():
    Args = collections.namedtuple('Args', ['escherichia_mlst_achtman_min_identity',
                                           'escherichia_mlst_achtman_min_coverage',
                                           'escherichia_mlst_achtman_required_exact_matches'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(escherichia_mlst_achtman_min_identity=90.0,
                               escherichia_mlst_achtman_min_coverage=120.0,
                               escherichia_mlst_achtman_required_exact_matches=3))


def test_check_cli_options_6():
    Args = collections.namedtuple('Args', ['escherichia_mlst_achtman_min_identity',
                                           'escherichia_mlst_achtman_min_coverage',
                                           'escherichia_mlst_achtman_required_exact_matches'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(escherichia_mlst_achtman_min_identity=90.0,
                               escherichia_mlst_achtman_min_coverage=90.0,
                               escherichia_mlst_achtman_required_exact_matches=-2))


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
    Args = collections.namedtuple('Args', ['escherichia_mlst_achtman_min_identity',
                                           'escherichia_mlst_achtman_min_coverage',
                                           'escherichia_mlst_achtman_required_exact_matches'])
    results = get_results(get_test_genome_dir() / 'GCF_000005845.2.fna.gz', None,
                          Args(escherichia_mlst_achtman_min_identity=90.0,
                               escherichia_mlst_achtman_min_coverage=80.0,
                               escherichia_mlst_achtman_required_exact_matches=3), {})
    assert results['ST'] == 'ST10'
    assert results['clonal_complex'] == 'ST10 Cplx'
    assert results['adk'] == '10'
    assert results['fumC'] == '11'
    assert results['gyrB'] == '4'
    assert results['icd'] == '8'
    assert results['mdh'] == '8'
    assert results['purA'] == '8'
    assert results['recA'] == '2'


def test_get_results_2():
    Args = collections.namedtuple('Args', ['escherichia_mlst_achtman_min_identity',
                                           'escherichia_mlst_achtman_min_coverage',
                                           'escherichia_mlst_achtman_required_exact_matches'])
    results = get_results(get_test_genome_dir() / 'GCF_000008865.2.fna.gz', None,
                          Args(escherichia_mlst_achtman_min_identity=90.0,
                               escherichia_mlst_achtman_min_coverage=80.0,
                               escherichia_mlst_achtman_required_exact_matches=3), {})
    assert results['ST'] == 'ST11'
    assert results['clonal_complex'] == 'ST11 Cplx'
    assert results['adk'] == '12'
    assert results['fumC'] == '12'
    assert results['gyrB'] == '8'
    assert results['icd'] == '12'
    assert results['mdh'] == '15'
    assert results['purA'] == '2'
    assert results['recA'] == '2'


def test_get_results_3():
    # Tests a Klebsiella pneumoniae using the Escherichia scheme, so no ST should be assigned.
    Args = collections.namedtuple('Args', ['escherichia_mlst_achtman_min_identity',
                                           'escherichia_mlst_achtman_min_coverage',
                                           'escherichia_mlst_achtman_required_exact_matches'])
    results = get_results(get_test_genome_dir() / 'GCF_000009885.1.fna.gz', None,
                          Args(escherichia_mlst_achtman_min_identity=90.0,
                               escherichia_mlst_achtman_min_coverage=80.0,
                               escherichia_mlst_achtman_required_exact_matches=3), {})
    assert results['ST'] == 'NA'
    assert results['clonal_complex'] == '-'
