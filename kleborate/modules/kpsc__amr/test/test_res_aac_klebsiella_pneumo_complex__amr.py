
"""
This file contains tests for Kleborate. To run all tests, go the repo's root directory and run:
  python3 -m pytest

To get code coverage stats:
  coverage run --source . -m pytest && coverage report -m

Copyright 2024 Kat Holt
Copyright 2024 Ryan Wick (rrwick@gmail.com)
Copyright 2024 (gathonimaranga@gmail.com)
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
import pathlib


from kleborate.shared.resMinimap import read_class_file, get_res_headers, resminimap_assembly
from kleborate.modules.klebsiella_pneumo_complex__amr.klebsiella_pneumo_complex__amr import get_headers, get_results

def get_test_genome_dir():
    return pathlib.Path(__file__).parents[4] / 'test' / 'test_res_aac'

def test_get_results_1():
    Args = collections.namedtuple('Args', ['klebsiella_pneumo_complex__amr_min_identity', 'klebsiella_pneumo_complex__amr_min_coverage','klebsiella_pneumo_complex__amr_min_spurious_identity', 'klebsiella_pneumo_complex__amr_min_spurious_coverage'])
    results = get_results(get_test_genome_dir() / '01.fasta', None,
                          Args(klebsiella_pneumo_complex__amr_min_identity=90.0, klebsiella_pneumo_complex__amr_min_coverage=80.0,klebsiella_pneumo_complex__amr_min_spurious_identity=80.0, klebsiella_pneumo_complex__amr_min_spurious_coverage=40.0), {})
    assert results['AGly_acquired'] == '-'


def test_get_results_2():
    Args = collections.namedtuple('Args', ['klebsiella_pneumo_complex__amr_min_identity', 'klebsiella_pneumo_complex__amr_min_coverage','klebsiella_pneumo_complex__amr_min_spurious_identity', 'klebsiella_pneumo_complex__amr_min_spurious_coverage'])
    results = get_results(get_test_genome_dir() / '02.fasta', None,
                          Args(klebsiella_pneumo_complex__amr_min_identity=90.0, klebsiella_pneumo_complex__amr_min_coverage=80.0,klebsiella_pneumo_complex__amr_min_spurious_identity=80.0, klebsiella_pneumo_complex__amr_min_spurious_coverage=40.0), {})
    assert results['AGly_acquired'] == "aac(6')-31"


def test_get_results_3():
    Args = collections.namedtuple('Args', ['klebsiella_pneumo_complex__amr_min_identity', 'klebsiella_pneumo_complex__amr_min_coverage','klebsiella_pneumo_complex__amr_min_spurious_identity', 'klebsiella_pneumo_complex__amr_min_spurious_coverage'])
    results = get_results(get_test_genome_dir() / '03.fasta', None,
                          Args(klebsiella_pneumo_complex__amr_min_identity=90.0, klebsiella_pneumo_complex__amr_min_coverage=80.0,klebsiella_pneumo_complex__amr_min_spurious_identity=80.0, klebsiella_pneumo_complex__amr_min_spurious_coverage=40.0), {})
    assert results['AGly_acquired'] == "aac(6')-31*"

def test_get_results_4():

    """
    This is a tricky one (and the one which first highlighted the problem). The best bit-score
    match is for a partial hit to Ant3''Ih-Aac6-IId, but the correct answer is Aac6Ib-cr, which
    has a full coverage hit
    """
    Args = collections.namedtuple('Args', ['klebsiella_pneumo_complex__amr_min_identity', 'klebsiella_pneumo_complex__amr_min_coverage','klebsiella_pneumo_complex__amr_min_spurious_identity', 'klebsiella_pneumo_complex__amr_min_spurious_coverage'])
    results = get_results(get_test_genome_dir() / '04.fasta', None,
                          Args(klebsiella_pneumo_complex__amr_min_identity=90.0, klebsiella_pneumo_complex__amr_min_coverage=80.0,klebsiella_pneumo_complex__amr_min_spurious_identity=80.0, klebsiella_pneumo_complex__amr_min_spurious_coverage=40.0), {})
    assert results['AGly_acquired'] == "aac(6')-Ib-cr.v2"


def test_get_results_5():
    """
    Same as test_get_results_4, but with the hit on the other strand
    """
    Args = collections.namedtuple('Args', ['klebsiella_pneumo_complex__amr_min_identity', 'klebsiella_pneumo_complex__amr_min_coverage','klebsiella_pneumo_complex__amr_min_spurious_identity', 'klebsiella_pneumo_complex__amr_min_spurious_coverage'])
    results = get_results(get_test_genome_dir() / '05.fasta', None,
                          Args(klebsiella_pneumo_complex__amr_min_identity=90.0, klebsiella_pneumo_complex__amr_min_coverage=80.0,klebsiella_pneumo_complex__amr_min_spurious_identity=80.0, klebsiella_pneumo_complex__amr_min_spurious_coverage=40.0), {})
    assert results['AGly_acquired'] == "aac(6')-Ib-cr.v2"