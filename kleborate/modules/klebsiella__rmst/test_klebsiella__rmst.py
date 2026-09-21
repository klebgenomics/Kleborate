"""
This file contains tests for Kleborate. To run all tests, go the repo's root directory and run:
  python3 -m pytest

To get code coverage stats:
  coverage run --source . -m pytest && coverage report -m

Copyright 2026 Mary Maranga (gathonimaranga@gmail.com)
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
import pathlib
import pytest
import sys

from .klebsiella__rmst import *


def get_test_genome_dir():
    return pathlib.Path(__file__).parents[3] / 'test' / 'test_genomes'


def get_test_file_dir():
    return pathlib.Path(__file__).parents[0] / 'test_files'


def get_file_or_skip(filename):
    """Finds test genome across possible paths and extensions, or skips if not found."""
    candidate_dirs = [get_test_file_dir(), get_test_genome_dir()]
    stem = filename.replace('.fasta', '').replace('.fna', '').replace('.gz', '')
    extensions = ['.fasta', '.fasta.gz', '.fna', '.fna.gz']

    for d in candidate_dirs:
        # direct check
        if (d / filename).is_file():
            return d / filename
        # check extensions
        for ext in extensions:
            if (d / (stem + ext)).is_file():
                return d / (stem + ext)

    pytest.skip(f"Test genome file {filename} not present in test directories.")


def test_prerequisite_modules():
    assert prerequisite_modules() == []


def test_check_cli_options_1():
    Args = collections.namedtuple('Args', ['klebsiella__rmst_min_identity', 'klebsiella__rmst_min_coverage',
                                           'klebsiella__rmst_min_spurious_identity', 'klebsiella__rmst_min_spurious_coverage',
                                           'klebsiella__rmst_required_exact_matches', 'klebsiella__rmst_min_gene_count'])
    check_cli_options(Args(klebsiella__rmst_min_identity=90.0, klebsiella__rmst_min_coverage=80.0,
                           klebsiella__rmst_min_spurious_identity=80.0, klebsiella__rmst_min_spurious_coverage=40.0,
                           klebsiella__rmst_required_exact_matches=1, klebsiella__rmst_min_gene_count=2))


def test_check_cli_options_2():
    Args = collections.namedtuple('Args', ['klebsiella__rmst_min_identity', 'klebsiella__rmst_min_coverage',
                                           'klebsiella__rmst_min_spurious_identity', 'klebsiella__rmst_min_spurious_coverage',
                                           'klebsiella__rmst_required_exact_matches', 'klebsiella__rmst_min_gene_count'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(klebsiella__rmst_min_identity=0.90, klebsiella__rmst_min_coverage=80.0,
                               klebsiella__rmst_min_spurious_identity=80.0, klebsiella__rmst_min_spurious_coverage=40.0,
                               klebsiella__rmst_required_exact_matches=1, klebsiella__rmst_min_gene_count=2))


def test_check_cli_options_3():
    Args = collections.namedtuple('Args', ['klebsiella__rmst_min_identity', 'klebsiella__rmst_min_coverage',
                                           'klebsiella__rmst_min_spurious_identity', 'klebsiella__rmst_min_spurious_coverage',
                                           'klebsiella__rmst_required_exact_matches', 'klebsiella__rmst_min_gene_count'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(klebsiella__rmst_min_identity=-90.0, klebsiella__rmst_min_coverage=0.90,
                               klebsiella__rmst_min_spurious_identity=80.0, klebsiella__rmst_min_spurious_coverage=40.0,
                               klebsiella__rmst_required_exact_matches=1, klebsiella__rmst_min_gene_count=2))


def test_check_cli_options_4():
    Args = collections.namedtuple('Args', ['klebsiella__rmst_min_identity', 'klebsiella__rmst_min_coverage',
                                           'klebsiella__rmst_min_spurious_identity', 'klebsiella__rmst_min_spurious_coverage',
                                           'klebsiella__rmst_required_exact_matches', 'klebsiella__rmst_min_gene_count'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(klebsiella__rmst_min_identity=-10.0, klebsiella__rmst_min_coverage=80.0,
                               klebsiella__rmst_min_spurious_identity=80.0, klebsiella__rmst_min_spurious_coverage=40.0,
                               klebsiella__rmst_required_exact_matches=1, klebsiella__rmst_min_gene_count=2))


def test_check_cli_options_5():
    Args = collections.namedtuple('Args', ['klebsiella__rmst_min_identity', 'klebsiella__rmst_min_coverage',
                                           'klebsiella__rmst_min_spurious_identity', 'klebsiella__rmst_min_spurious_coverage',
                                           'klebsiella__rmst_required_exact_matches', 'klebsiella__rmst_min_gene_count'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(klebsiella__rmst_min_identity=90.0, klebsiella__rmst_min_coverage=120.0,
                               klebsiella__rmst_min_spurious_identity=80.0, klebsiella__rmst_min_spurious_coverage=40.0,
                               klebsiella__rmst_required_exact_matches=1, klebsiella__rmst_min_gene_count=2))


def test_check_cli_options_6():
    Args = collections.namedtuple('Args', ['klebsiella__rmst_min_identity', 'klebsiella__rmst_min_coverage',
                                           'klebsiella__rmst_min_spurious_identity', 'klebsiella__rmst_min_spurious_coverage',
                                           'klebsiella__rmst_required_exact_matches', 'klebsiella__rmst_min_gene_count'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(klebsiella__rmst_min_identity=90.0, klebsiella__rmst_min_coverage=80.0,
                               klebsiella__rmst_min_spurious_identity=-10.0, klebsiella__rmst_min_spurious_coverage=40.0,
                               klebsiella__rmst_required_exact_matches=-1, klebsiella__rmst_min_gene_count=2))


def test_check_cli_options_7():
    Args = collections.namedtuple('Args', ['klebsiella__rmst_min_identity', 'klebsiella__rmst_min_coverage',
                                           'klebsiella__rmst_min_spurious_identity', 'klebsiella__rmst_min_spurious_coverage',
                                           'klebsiella__rmst_required_exact_matches', 'klebsiella__rmst_min_gene_count'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(klebsiella__rmst_min_identity=90.0, klebsiella__rmst_min_coverage=80.0,
                               klebsiella__rmst_min_spurious_identity=80.0, klebsiella__rmst_min_spurious_coverage=40.0,
                               klebsiella__rmst_required_exact_matches=1, klebsiella__rmst_min_gene_count=-1))


def test_check_external_programs_success(mocker):
    # Tests the good case where rammappy is successfully imported.
    mock_module = mocker.MagicMock()
    mocker.patch.dict(sys.modules, {'rammappy': mock_module})

    assert check_external_programs() == ['rammappy']


def test_check_external_programs_import_error(mocker):
    # Tests the bad case where rammappy cannot be imported.
    mocker.patch.dict(sys.modules, {'rammappy': None})

    with pytest.raises(SystemExit) as exc_info:
        check_external_programs()

    assert 'Error: could not import rammappy' in str(exc_info.value)


def test_get_results_1():
    Args = collections.namedtuple('Args', ['klebsiella__rmst_min_identity', 'klebsiella__rmst_min_coverage',
                                           'klebsiella__rmst_min_spurious_identity', 'klebsiella__rmst_min_spurious_coverage',
                                           'klebsiella__rmst_required_exact_matches', 'klebsiella__rmst_min_gene_count'])
    results = get_results(get_test_genome_dir() / 'GCF_000968155.1.fna.gz', None,
                          Args(klebsiella__rmst_min_identity=90.0, klebsiella__rmst_min_coverage=80.0,
                               klebsiella__rmst_min_spurious_identity=80.0, klebsiella__rmst_min_spurious_coverage=40.0,
                               klebsiella__rmst_required_exact_matches=1, klebsiella__rmst_min_gene_count=2), {})

    assert results['RmST'] == '2'
    assert results['RmpADC'] == 'rmp 2; KpVP-2'
    assert results['RmpADC_status'] == 'OFF (reduced)'
    assert results['rmpA'] == '9'
    assert results['rmpD'] == '32'
    assert results['rmpC'] == '5 (OFF)'
    assert results['rmpA_promoter'] == '10T (reduced expression)'
    assert results['argR'] in ('present', ['present'])


def test_get_results_2():
    Args = collections.namedtuple('Args', ['klebsiella__rmst_min_identity', 'klebsiella__rmst_min_coverage',
                                           'klebsiella__rmst_min_spurious_identity', 'klebsiella__rmst_min_spurious_coverage',
                                           'klebsiella__rmst_required_exact_matches', 'klebsiella__rmst_min_gene_count'])
    results = get_results(get_test_genome_dir() / 'GCF_001068035.1.fna.gz', None,
                          Args(klebsiella__rmst_min_identity=90.0, klebsiella__rmst_min_coverage=80.0,
                               klebsiella__rmst_min_spurious_identity=80.0, klebsiella__rmst_min_spurious_coverage=40.0,
                               klebsiella__rmst_required_exact_matches=1, klebsiella__rmst_min_gene_count=2), {})
    assert results['RmST'] == '26'
    assert results['RmpADC'] == 'rmp 1; KpVP-1'
    assert results['RmpADC_status'] == 'ON'
    assert results['rmpD'] == '2'
    assert results['rmpA'] == '2'
    assert results['rmpC'] == '2'
    assert results['rmpA_promoter'] == '11T'
    assert results['argR'] in ('present', ['present'])


def test_get_results_3():
    # Tests an E. coli without the rmp locus, so no ST should be assigned.
    Args = collections.namedtuple('Args', ['klebsiella__rmst_min_identity', 'klebsiella__rmst_min_coverage',
                                           'klebsiella__rmst_min_spurious_identity', 'klebsiella__rmst_min_spurious_coverage',
                                           'klebsiella__rmst_required_exact_matches', 'klebsiella__rmst_min_gene_count'])
    results = get_results(get_test_genome_dir() / 'GCF_000008865.2.fna.gz', None,
                          Args(klebsiella__rmst_min_identity=90.0, klebsiella__rmst_min_coverage=80.0,
                               klebsiella__rmst_min_spurious_identity=80.0, klebsiella__rmst_min_spurious_coverage=40.0,
                               klebsiella__rmst_required_exact_matches=1, klebsiella__rmst_min_gene_count=2), {})
    assert results['RmST'] == 0
    assert results['RmpADC'] == '-'


def test_gcf_000968155():
    Args = collections.namedtuple('Args', ['klebsiella__rmst_min_identity', 'klebsiella__rmst_min_coverage',
                                           'klebsiella__rmst_min_spurious_identity', 'klebsiella__rmst_min_spurious_coverage',
                                           'klebsiella__rmst_required_exact_matches', 'klebsiella__rmst_min_gene_count'])
    results = get_results(get_test_file_dir() / 'GCF_000968155.1_rmp_locus.fasta', None,
                          Args(klebsiella__rmst_min_identity=90.0, klebsiella__rmst_min_coverage=80.0,
                               klebsiella__rmst_min_spurious_identity=80.0, klebsiella__rmst_min_spurious_coverage=40.0,
                               klebsiella__rmst_required_exact_matches=1, klebsiella__rmst_min_gene_count=2), {})
    assert results['RmST'] == '2'
    assert results['RmpADC'] == 'rmp 2; KpVP-2'
    assert results['RmpADC_status'] == 'OFF (reduced), (argR missing)'
    assert results['rmpA'] == '9'
    assert results['rmpD'] == '32'
    assert results['rmpC'] == '5 (OFF)'
    assert results['rmpA_promoter'] == '10T (reduced expression)'
    assert results['argR'] in ('-', ['-'])


def test_gcf_000968155_incomplete():
    Args = collections.namedtuple('Args', ['klebsiella__rmst_min_identity', 'klebsiella__rmst_min_coverage',
                                           'klebsiella__rmst_min_spurious_identity', 'klebsiella__rmst_min_spurious_coverage',
                                           'klebsiella__rmst_required_exact_matches', 'klebsiella__rmst_min_gene_count'])
    results = get_results(get_test_file_dir() / 'GCF_000968155.1_rmp_locus_incomplete.fasta', None,
                          Args(klebsiella__rmst_min_identity=90.0, klebsiella__rmst_min_coverage=80.0,
                               klebsiella__rmst_min_spurious_identity=80.0, klebsiella__rmst_min_spurious_coverage=40.0,
                               klebsiella__rmst_required_exact_matches=1, klebsiella__rmst_min_gene_count=2), {})
    assert results['RmST'] == '2-1LV'
    assert results['RmpADC'] == 'rmp 2; KpVP-2 (partial)'
    assert results['RmpADC_status'] == '-'
    assert results['rmpA'] == '9'
    assert results['rmpD'] == '32'
    assert results['rmpC'] == '-'
    assert results['rmpA_promoter'] == '10T (reduced expression)'
    assert results['argR'] in ('-', ['-'])


def test_gcf_000968155_truncated():
    Args = collections.namedtuple('Args', ['klebsiella__rmst_min_identity', 'klebsiella__rmst_min_coverage',
                                           'klebsiella__rmst_min_spurious_identity', 'klebsiella__rmst_min_spurious_coverage',
                                           'klebsiella__rmst_required_exact_matches', 'klebsiella__rmst_min_gene_count'])
    results = get_results(get_test_file_dir() / 'GCF_000968155.1_rmp_locus_truncated.fasta', None,
                          Args(klebsiella__rmst_min_identity=90.0, klebsiella__rmst_min_coverage=80.0,
                               klebsiella__rmst_min_spurious_identity=80.0, klebsiella__rmst_min_spurious_coverage=40.0,
                               klebsiella__rmst_required_exact_matches=1, klebsiella__rmst_min_gene_count=2), {})
    
    assert results['RmST'] == '2-1LV'
    assert results['RmpADC'] == 'rmp 2; KpVP-2 (partial)'
    assert results['RmpADC_status'] == 'OFF (reduced), (argR missing)'
    assert results['rmpA'] == '9'
    assert results['rmpD'] == '32'
    assert results['rmpC'] == '5*-28% (OFF)'
    assert results['rmpA_promoter'] == '10T (reduced expression)'
    assert results['argR'] in ('-', ['-'])


def test_get_results_KP01():
    # Test for strain KP01
    Args = collections.namedtuple('Args', ['klebsiella__rmst_min_identity', 'klebsiella__rmst_min_coverage',
                                           'klebsiella__rmst_min_spurious_identity', 'klebsiella__rmst_min_spurious_coverage',
                                           'klebsiella__rmst_required_exact_matches', 'klebsiella__rmst_min_gene_count'])
    results = get_results(get_test_genome_dir() / 'KP01.fasta', None,
                          Args(klebsiella__rmst_min_identity=90.0, klebsiella__rmst_min_coverage=80.0,
                               klebsiella__rmst_min_spurious_identity=80.0, klebsiella__rmst_min_spurious_coverage=40.0,
                               klebsiella__rmst_required_exact_matches=1, klebsiella__rmst_min_gene_count=2), {})

    assert results['RmST'] == '40'
    assert results['RmpADC'] == 'rmp 1; KpVP-1'
    assert results['RmpADC_status'] == 'ON'
    assert results['rmpA'] == '2'
    assert results['rmpD'] == '3'
    assert results['rmpC'] == '2'
    assert results['rmpA_promoter'] == '11T'
    assert results['argR'] in ('present', ['present'])


def test_get_results_KP02():
    # Test for strain KP02
    Args = collections.namedtuple('Args', ['klebsiella__rmst_min_identity', 'klebsiella__rmst_min_coverage',
                                           'klebsiella__rmst_min_spurious_identity', 'klebsiella__rmst_min_spurious_coverage',
                                           'klebsiella__rmst_required_exact_matches', 'klebsiella__rmst_min_gene_count'])
    results = get_results(get_test_genome_dir() / 'KP02.fasta', None,
                          Args(klebsiella__rmst_min_identity=90.0, klebsiella__rmst_min_coverage=80.0,
                               klebsiella__rmst_min_spurious_identity=80.0, klebsiella__rmst_min_spurious_coverage=40.0,
                               klebsiella__rmst_required_exact_matches=1, klebsiella__rmst_min_gene_count=2), {})

    assert results['RmST'] == '40'
    assert results['RmpADC'] == 'rmp 1; KpVP-1'
    assert results['RmpADC_status'] == 'ON (reduced)'
    assert results['rmpA'] == '2'
    assert results['rmpD'] == '3'
    assert results['rmpC'] == '2'
    assert results['rmpA_promoter'] == '10T (reduced expression)'
    assert results['argR'] in ('present', ['present'])


def test_get_results_DRR389032_argR_box():
    Args = collections.namedtuple('Args', ['klebsiella__rmst_min_identity', 'klebsiella__rmst_min_coverage',
                                           'klebsiella__rmst_min_spurious_identity', 'klebsiella__rmst_min_spurious_coverage',
                                           'klebsiella__rmst_required_exact_matches', 'klebsiella__rmst_min_gene_count'])
    target_file = get_file_or_skip('DRR389032.fasta')
    results = get_results(target_file, None,
                          Args(klebsiella__rmst_min_identity=90.0, klebsiella__rmst_min_coverage=80.0,
                               klebsiella__rmst_min_spurious_identity=80.0, klebsiella__rmst_min_spurious_coverage=40.0,
                               klebsiella__rmst_required_exact_matches=1, klebsiella__rmst_min_gene_count=2), {})

    assert results['RmST'] == '26'
    assert results['RmpADC'] == 'rmp 1; KpVP-1'
    assert results['RmpADC_status'] == 'ON (ARG box lost)'
    assert results['rmpA'] == '2'
    assert results['rmpD'] == '2'
    assert results['rmpC'] == '2'
    assert results['rmpA_promoter'] == '12T, ARG-box lost'
    assert results['argR'] in ('present', ['present'])


def test_get_results_NK_H1_084_argR_box():
    Args = collections.namedtuple('Args', ['klebsiella__rmst_min_identity', 'klebsiella__rmst_min_coverage',
                                           'klebsiella__rmst_min_spurious_identity', 'klebsiella__rmst_min_spurious_coverage',
                                           'klebsiella__rmst_required_exact_matches', 'klebsiella__rmst_min_gene_count'])
    target_file = get_file_or_skip('NK_H1_084.fasta')
    results = get_results(target_file, None,
                          Args(klebsiella__rmst_min_identity=90.0, klebsiella__rmst_min_coverage=80.0,
                               klebsiella__rmst_min_spurious_identity=80.0, klebsiella__rmst_min_spurious_coverage=40.0,
                               klebsiella__rmst_required_exact_matches=1, klebsiella__rmst_min_gene_count=2), {})

    assert results['RmST'] == '40'
    assert results['RmpADC'] == 'rmp 1; KpVP-1'
    assert results['RmpADC_status'] == 'ON (ARG box lost)'
    assert results['rmpA'] == '2'
    assert results['rmpD'] == '3'
    assert results['rmpC'] == '2'
    assert results['rmpA_promoter'] == '11T (ARG-box lost)'
    assert results['argR'] in ('present', ['present'])