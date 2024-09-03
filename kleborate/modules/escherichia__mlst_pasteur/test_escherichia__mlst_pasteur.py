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

from .escherichia__mlst_pasteur import *


def get_test_genome_dir():
    return pathlib.Path(__file__).parents[3] / 'test' / 'test_genomes'


def test_prerequisite_modules():
    assert prerequisite_modules() == []


def test_check_cli_options_1():
    Args = collections.namedtuple('Args', ['escherichia_mlst_pasteur_min_identity',
                                           'escherichia_mlst_pasteur_min_coverage',
                                           'escherichia_mlst_pasteur_required_exact_matches'])
    check_cli_options(Args(escherichia_mlst_pasteur_min_identity=90.0,
                           escherichia_mlst_pasteur_min_coverage=90.0,
                           escherichia_mlst_pasteur_required_exact_matches=3))


def test_check_cli_options_2():
    Args = collections.namedtuple('Args', ['escherichia_mlst_pasteur_min_identity',
                                           'escherichia_mlst_pasteur_min_coverage',
                                           'escherichia_mlst_pasteur_required_exact_matches'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(escherichia_mlst_pasteur_min_identity=0.90,
                               escherichia_mlst_pasteur_min_coverage=90.0,
                               escherichia_mlst_pasteur_required_exact_matches=3))


def test_check_cli_options_3():
    Args = collections.namedtuple('Args', ['escherichia_mlst_pasteur_min_identity',
                                           'escherichia_mlst_pasteur_min_coverage',
                                           'escherichia_mlst_pasteur_required_exact_matches'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(escherichia_mlst_pasteur_min_identity=-90.0,
                               escherichia_mlst_pasteur_min_coverage=0.90,
                               escherichia_mlst_pasteur_required_exact_matches=3))


def test_check_cli_options_4():
    Args = collections.namedtuple('Args', ['escherichia_mlst_pasteur_min_identity',
                                           'escherichia_mlst_pasteur_min_coverage',
                                           'escherichia_mlst_pasteur_required_exact_matches'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(escherichia_mlst_pasteur_min_identity=-10.0,
                               escherichia_mlst_pasteur_min_coverage=90.0,
                               escherichia_mlst_pasteur_required_exact_matches=3))


def test_check_cli_options_5():
    Args = collections.namedtuple('Args', ['escherichia_mlst_pasteur_min_identity',
                                           'escherichia_mlst_pasteur_min_coverage',
                                           'escherichia_mlst_pasteur_required_exact_matches'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(escherichia_mlst_pasteur_min_identity=90.0,
                               escherichia_mlst_pasteur_min_coverage=120.0,
                               escherichia_mlst_pasteur_required_exact_matches=3))


def test_check_cli_options_6():
    Args = collections.namedtuple('Args', ['escherichia_mlst_pasteur_min_identity',
                                           'escherichia_mlst_pasteur_min_coverage',
                                           'escherichia_mlst_pasteur_required_exact_matches'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(escherichia_mlst_pasteur_min_identity=90.0,
                               escherichia_mlst_pasteur_min_coverage=90.0,
                               escherichia_mlst_pasteur_required_exact_matches=-2))


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
    Args = collections.namedtuple('Args', ['escherichia_mlst_pasteur_min_identity',
                                           'escherichia_mlst_pasteur_min_coverage',
                                           'escherichia_mlst_pasteur_required_exact_matches'])
    results = get_results(get_test_genome_dir() / 'GCF_000005845.2.fna.gz', None,
                          Args(escherichia_mlst_pasteur_min_identity=90.0,
                               escherichia_mlst_pasteur_min_coverage=80.0,
                               escherichia_mlst_pasteur_required_exact_matches=3), {})

    assert results['ST'] == 'ST262'
    assert results['dinB'] == '8'
    assert results['icdA'] == '118'
    assert results['pabB'] == '7'
    assert results['polB'] == '3'
    assert results['putP'] == '7'
    assert results['trpA'] == '1'
    assert results['trpB'] == '4'
    assert results['uidA'] == '2'


def test_get_results_2():
    Args = collections.namedtuple('Args', ['escherichia_mlst_pasteur_min_identity',
                                           'escherichia_mlst_pasteur_min_coverage',
                                           'escherichia_mlst_pasteur_required_exact_matches'])
    results = get_results(get_test_genome_dir() / 'GCF_000008865.2.fna.gz', None,
                          Args(escherichia_mlst_pasteur_min_identity=90.0,
                               escherichia_mlst_pasteur_min_coverage=80.0,
                               escherichia_mlst_pasteur_required_exact_matches=3), {})
    assert results['ST'] == 'ST296'
    assert results['dinB'] == '68'
    assert results['icdA'] == '110'
    assert results['pabB'] == '91'
    assert results['polB'] == '63'
    assert results['putP'] == '93'
    assert results['trpA'] == '84'
    assert results['trpB'] == '86'
    assert results['uidA'] == '83'


def test_get_results_3():
    # Tests a Klebsiella pneumoniae using the Escherichia scheme, so no ST should be assigned.
    Args = collections.namedtuple('Args', ['escherichia_mlst_pasteur_min_identity',
                                           'escherichia_mlst_pasteur_min_coverage',
                                           'escherichia_mlst_pasteur_required_exact_matches'])
    results = get_results(get_test_genome_dir() / 'GCF_000009885.1.fna.gz', None,
                          Args(escherichia_mlst_pasteur_min_identity=90.0,
                               escherichia_mlst_pasteur_min_coverage=80.0,
                               escherichia_mlst_pasteur_required_exact_matches=3), {})
    assert results['ST'] == 'NA'
