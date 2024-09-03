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

from .klebsiella__rmst import *


def get_test_genome_dir():
    return pathlib.Path(__file__).parents[3] / 'test' / 'test_genomes'


def get_test_file_dir():
    return pathlib.Path(__file__).parents[0] / 'test_files'


def test_prerequisite_modules():
    assert prerequisite_modules() == []


def test_check_cli_options_1():
    Args = collections.namedtuple('Args', ['klebsiella__rmst_min_identity', 'klebsiella__rmst_min_coverage',
                                           'klebsiella__rmst_min_spurious_identity', 'klebsiella__rmst_min_spurious_coverage',
                                           'klebsiella__rmst_required_exact_matches'])
    check_cli_options(Args(klebsiella__rmst_min_identity=90.0, klebsiella__rmst_min_coverage=90.0,
                           klebsiella__rmst_min_spurious_identity=85.0, klebsiella__rmst_min_spurious_coverage=50.0,
                           klebsiella__rmst_required_exact_matches=3))


def test_check_cli_options_2():
    Args = collections.namedtuple('Args', ['klebsiella__rmst_min_identity', 'klebsiella__rmst_min_coverage',
                                           'klebsiella__rmst_required_exact_matches'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(klebsiella__rmst_min_identity=0.90, klebsiella__rmst_min_coverage=90.0,
                               klebsiella__rmst_required_exact_matches=3))


def test_check_cli_options_3():
    Args = collections.namedtuple('Args', ['klebsiella__rmst_min_identity', 'klebsiella__rmst_min_coverage',
                                           'klebsiella__rmst_required_exact_matches'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(klebsiella__rmst_min_identity=-90.0, klebsiella__rmst_min_coverage=0.90,
                               klebsiella__rmst_required_exact_matches=3))


def test_check_cli_options_4():
    Args = collections.namedtuple('Args', ['klebsiella__rmst_min_identity', 'klebsiella__rmst_min_coverage',
                                           'klebsiella__rmst_required_exact_matches'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(klebsiella__rmst_min_identity=-10.0, klebsiella__rmst_min_coverage=90.0,
                               klebsiella__rmst_required_exact_matches=3))


def test_check_cli_options_5():
    Args = collections.namedtuple('Args', ['klebsiella__rmst_min_identity', 'klebsiella__rmst_min_coverage',
                                           'klebsiella__rmst_required_exact_matches'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(klebsiella__rmst_min_identity=90.0, klebsiella__rmst_min_coverage=120.0,
                               klebsiella__rmst_required_exact_matches=3))


def test_check_cli_options_6():
    Args = collections.namedtuple('Args', ['klebsiella__rmst_min_identity', 'klebsiella__rmst_min_coverage',
                                           'klebsiella__rmst_min_spurious_identity', 'klebsiella__rmst_min_spurious_coverage',
                                           'klebsiella__rmst_required_exact_matches'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(klebsiella__rmst_min_identity=90.0, klebsiella__rmst_min_coverage=90.0,
                               klebsiella__rmst_min_spurious_identity=80.0, klebsiella__rmst_min_spurious_coverage=50.0,
                               klebsiella__rmst_required_exact_matches=-2))



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
    Args = collections.namedtuple('Args', ['klebsiella__rmst_min_identity', 'klebsiella__rmst_min_coverage',
                                           'klebsiella__rmst_min_spurious_identity', 'klebsiella__rmst_min_spurious_coverage',
                                           'klebsiella__rmst_required_exact_matches'])
    results = get_results(get_test_genome_dir() / 'GCF_000968155.1.fna.gz', None,
                          Args(klebsiella__rmst_min_identity=90.0, klebsiella__rmst_min_coverage=80.0,
                               klebsiella__rmst_min_spurious_identity=80.0, klebsiella__rmst_min_spurious_coverage=40.0,
                               klebsiella__rmst_required_exact_matches=3), {})

    assert results['RmST'] == '2'
    assert results['RmpADC'] == 'rmp 2; KpVP-2'
    assert results['rmpA'] == '9'
    assert results['rmpD'] == '32'
    assert results['rmpC'] == '5'


def test_get_results_2():
    Args = collections.namedtuple('Args', ['klebsiella__rmst_min_identity', 'klebsiella__rmst_min_coverage',
                                           'klebsiella__rmst_min_spurious_identity', 'klebsiella__rmst_min_spurious_coverage',
                                           'klebsiella__rmst_required_exact_matches'])
    results = get_results(get_test_genome_dir() / 'GCF_001068035.1.fna.gz', None,
                          Args(klebsiella__rmst_min_identity=90.0, klebsiella__rmst_min_coverage=80.0,
                               klebsiella__rmst_min_spurious_identity=80.0, klebsiella__rmst_min_spurious_coverage=40.0,
                               klebsiella__rmst_required_exact_matches=3), {})
    assert results['RmST'] == '26'
    assert results['RmpADC'] == 'rmp 1; KpVP-1'
    assert results['rmpD'] == '2'
    assert results['rmpA'] == '2'
    assert results['rmpC'] == '2'


def test_get_results_3():
    # Tests an E. coli without the iro locus, so no ST should be assigned.
    Args = collections.namedtuple('Args', ['klebsiella__rmst_min_identity', 'klebsiella__rmst_min_coverage',
                                           'klebsiella__rmst_min_spurious_identity', 'klebsiella__rmst_min_spurious_coverage',
                                           'klebsiella__rmst_required_exact_matches'])
    results = get_results(get_test_genome_dir() / 'GCF_000008865.2.fna.gz', None,
                          Args(klebsiella__rmst_min_identity=90.0, klebsiella__rmst_min_coverage=80.0,
                               klebsiella__rmst_min_spurious_identity=80.0, klebsiella__rmst_min_spurious_coverage=40.0,
                               klebsiella__rmst_required_exact_matches=3), {})
    assert results['RmST'] == 0
    assert results['RmpADC'] == '-'


def test_gcf_000968155():
    Args = collections.namedtuple('Args', ['klebsiella__rmst_min_identity', 'klebsiella__rmst_min_coverage',
                                           'klebsiella__rmst_min_spurious_identity', 'klebsiella__rmst_min_spurious_coverage',
                                           'klebsiella__rmst_required_exact_matches'])
    results = get_results(get_test_file_dir() / 'GCF_000968155.1_rmp_locus.fasta', None,
                          Args(klebsiella__rmst_min_identity=90.0, klebsiella__rmst_min_coverage=80.0,
                               klebsiella__rmst_min_spurious_identity=80.0, klebsiella__rmst_min_spurious_coverage=40.0,
                               klebsiella__rmst_required_exact_matches=3), {})
    assert results['RmST'] == '2'
    assert results['RmpADC'] == 'rmp 2; KpVP-2'
    assert results['rmpA'] == '9'
    assert results['rmpD'] == '32'
    assert results['rmpC'] == '5'


def test_gcf_000968155_incomplete():
    Args = collections.namedtuple('Args', ['klebsiella__rmst_min_identity', 'klebsiella__rmst_min_coverage',
                                           'klebsiella__rmst_min_spurious_identity', 'klebsiella__rmst_min_spurious_coverage',
                                           'klebsiella__rmst_required_exact_matches'])
    results = get_results(get_test_file_dir() / 'GCF_000968155.1_rmp_locus_incomplete.fasta', None,
                          Args(klebsiella__rmst_min_identity=90.0, klebsiella__rmst_min_coverage=80.0,
                               klebsiella__rmst_min_spurious_identity=80.0, klebsiella__rmst_min_spurious_coverage=40.0,
                               klebsiella__rmst_required_exact_matches=3), {})
    assert results['RmST'] == 0
    assert results['RmpADC'] == '-'
    assert results['rmpA'] == '9'
    assert results['rmpD'] == '32'
    assert results['rmpC'] == '-'


def test_gcf_000968155_truncated():
    Args = collections.namedtuple('Args', ['klebsiella__rmst_min_identity', 'klebsiella__rmst_min_coverage',
                                           'klebsiella__rmst_min_spurious_identity', 'klebsiella__rmst_min_spurious_coverage',
                                           'klebsiella__rmst_required_exact_matches'])
    results = get_results(get_test_file_dir() / 'GCF_000968155.1_rmp_locus_truncated.fasta', None,
                          Args(klebsiella__rmst_min_identity=90.0, klebsiella__rmst_min_coverage=80.0,
                               klebsiella__rmst_min_spurious_identity=80.0, klebsiella__rmst_min_spurious_coverage=40.0,
                               klebsiella__rmst_required_exact_matches=3), {})
    assert results['RmST'] == 0
    assert results['RmpADC'] == '-'
    assert results['rmpA'] == '9'
    assert results['rmpD'] == '32'
    assert results['rmpC'] == '5*-28%'
