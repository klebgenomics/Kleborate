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

from .klebsiella_pneumo_complex__mlst import *


def get_test_genome_dir():
    return pathlib.Path(__file__).parents[3] / 'test' / 'test_genomes'


def test_prerequisite_modules():
    assert prerequisite_modules() == []


def test_check_cli_options_1():
    Args = collections.namedtuple('Args', ['klebsiella_pneumo_complex__mlst_min_identity', 'klebsiella_pneumo_complex__mlst_min_coverage',
                                           'klebsiella_pneumo_complex__mlst_required_exact_matches'])
    check_cli_options(Args(klebsiella_pneumo_complex__mlst_min_identity=90.0, klebsiella_pneumo_complex__mlst_min_coverage=90.0,
                           klebsiella_pneumo_complex__mlst_required_exact_matches=3))


def test_check_cli_options_2():
    Args = collections.namedtuple('Args', ['klebsiella_pneumo_complex__mlst_min_identity', 'klebsiella_pneumo_complex__mlst_min_coverage',
                                           'klebsiella_pneumo_complex__mlst_required_exact_matches'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(klebsiella_pneumo_complex__mlst_min_identity=0.90, klebsiella_pneumo_complex__mlst_min_coverage=90.0,
                               klebsiella_pneumo_complex__mlst_required_exact_matches=3))


def test_check_cli_options_3():
    Args = collections.namedtuple('Args', ['klebsiella_pneumo_complex__mlst_min_identity', 'klebsiella_pneumo_complex__mlst_min_coverage',
                                           'klebsiella_pneumo_complex__mlst_required_exact_matches'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(klebsiella_pneumo_complex__mlst_min_identity=-90.0, klebsiella_pneumo_complex__mlst_min_coverage=0.90,
                               klebsiella_pneumo_complex__mlst_required_exact_matches=3))


def test_check_cli_options_4():
    Args = collections.namedtuple('Args', ['klebsiella_pneumo_complex__mlst_min_identity', 'klebsiella_pneumo_complex__mlst_min_coverage',
                                           'klebsiella_pneumo_complex__mlst_required_exact_matches'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(klebsiella_pneumo_complex__mlst_min_identity=-10.0, klebsiella_pneumo_complex__mlst_min_coverage=90.0,
                               klebsiella_pneumo_complex__mlst_required_exact_matches=3))


def test_check_cli_options_5():
    Args = collections.namedtuple('Args', ['klebsiella_pneumo_complex__mlst_min_identity', 'klebsiella_pneumo_complex__mlst_min_coverage',
                                           'klebsiella_pneumo_complex__mlst_required_exact_matches'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(klebsiella_pneumo_complex__mlst_min_identity=90.0, klebsiella_pneumo_complex__mlst_min_coverage=120.0,
                               klebsiella_pneumo_complex__mlst_required_exact_matches=3))


def test_check_cli_options_6():
    Args = collections.namedtuple('Args', ['klebsiella_pneumo_complex__mlst_min_identity', 'klebsiella_pneumo_complex__mlst_min_coverage',
                                           'klebsiella_pneumo_complex__mlst_required_exact_matches'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(klebsiella_pneumo_complex__mlst_min_identity=90.0, klebsiella_pneumo_complex__mlst_min_coverage=90.0,
                               klebsiella_pneumo_complex__mlst_required_exact_matches=-2))


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
    Args = collections.namedtuple('Args', ['klebsiella_pneumo_complex__mlst_min_identity', 'klebsiella_pneumo_complex__mlst_min_coverage',
                                           'klebsiella_pneumo_complex__mlst_required_exact_matches'])
    results = get_results(get_test_genome_dir() / 'GCF_000968155.1.fna.gz', None,
                          Args(klebsiella_pneumo_complex__mlst_min_identity=90.0, klebsiella_pneumo_complex__mlst_min_coverage=80.0,
                               klebsiella_pneumo_complex__mlst_required_exact_matches=3), {})
    assert results['ST'] == 'ST66'
    assert results['gapA'] == '2'
    assert results['infB'] == '3'
    assert results['mdh'] == '2'
    assert results['pgi'] == '1'
    assert results['phoE'] == '10'
    assert results['rpoB'] == '1'
    assert results['tonB'] == '13'


def test_get_results_2():
    Args = collections.namedtuple('Args', ['klebsiella_pneumo_complex__mlst_min_identity', 'klebsiella_pneumo_complex__mlst_min_coverage',
                                           'klebsiella_pneumo_complex__mlst_required_exact_matches'])
    results = get_results(get_test_genome_dir() / 'GCF_001068035.1.fna.gz', None,
                          Args(klebsiella_pneumo_complex__mlst_min_identity=90.0, klebsiella_pneumo_complex__mlst_min_coverage=80.0,
                               klebsiella_pneumo_complex__mlst_required_exact_matches=3), {})
    assert results['ST'] == 'ST592-1LV'
    assert results['gapA'] == '2'
    assert results['infB'] == '3'
    assert results['mdh'] == '6'
    assert results['pgi'] == '-'
    assert results['phoE'] == '9'
    assert results['rpoB'] == '4'
    assert results['tonB'] == '13'


def test_get_results_3():
    # Tests a Klebsiella oxytoca using the KpSC scheme, so no ST should be assigned.
    Args = collections.namedtuple('Args', ['klebsiella_pneumo_complex__mlst_min_identity', 'klebsiella_pneumo_complex__mlst_min_coverage',
                                           'klebsiella_pneumo_complex__mlst_required_exact_matches'])
    results = get_results(get_test_genome_dir() / 'GCF_000247855.1.fna.gz', None,
                          Args(klebsiella_pneumo_complex__mlst_min_identity=90.0, klebsiella_pneumo_complex__mlst_min_coverage=80.0,
                               klebsiella_pneumo_complex__mlst_required_exact_matches=3), {})
    assert results['ST'] == 'NA'
