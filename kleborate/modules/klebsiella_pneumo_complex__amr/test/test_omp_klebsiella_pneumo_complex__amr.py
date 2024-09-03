
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
    return pathlib.Path(__file__).parents[4] / 'test' / 'test_res_omp' 

"""
Tests calling of carbapenem resistance via the OmpK35/OmpK36 genes.
"""
def test_get_results_1():
    Args = collections.namedtuple('Args', ['klebsiella_pneumo_complex__amr_min_identity', 'klebsiella_pneumo_complex__amr_min_coverage','klebsiella_pneumo_complex__amr_min_spurious_identity', 'klebsiella_pneumo_complex__amr_min_spurious_coverage'])
    results = get_results(get_test_genome_dir() / 'test_res_omp_1.fasta', None,
                          Args(klebsiella_pneumo_complex__amr_min_identity=90.0, klebsiella_pneumo_complex__amr_min_coverage=80.0,klebsiella_pneumo_complex__amr_min_spurious_identity=80.0, klebsiella_pneumo_complex__amr_min_spurious_coverage=40.0), {})
    assert results['Omp_mutations'] == '-'


def test_get_results_2():
    #A frameshift in OmpK35 should cause an early stop and lead to a carbapenem resistance call.
    Args = collections.namedtuple('Args', ['klebsiella_pneumo_complex__amr_min_identity', 'klebsiella_pneumo_complex__amr_min_coverage', 'klebsiella_pneumo_complex__amr_min_spurious_identity', 'klebsiella_pneumo_complex__amr_min_spurious_coverage'])
    results = get_results(get_test_genome_dir() / 'test_res_omp_2.fasta', None,
                          Args(klebsiella_pneumo_complex__amr_min_identity=90.0, klebsiella_pneumo_complex__amr_min_coverage=80.0, klebsiella_pneumo_complex__amr_min_spurious_identity=80.0, klebsiella_pneumo_complex__amr_min_spurious_coverage=40.0), {})
    assert results['Omp_mutations'] == 'OmpK35-36%;OmpK36_c25t'


def test_get_results_3():
    #This tests an early stop mutation (without a frameshift) in OmpK35.
    Args = collections.namedtuple('Args', ['klebsiella_pneumo_complex__amr_min_identity', 'klebsiella_pneumo_complex__amr_min_coverage','klebsiella_pneumo_complex__amr_min_spurious_identity', 'klebsiella_pneumo_complex__amr_min_spurious_coverage'])
    results = get_results(get_test_genome_dir() / 'test_res_omp_3.fasta', None,
                          Args(klebsiella_pneumo_complex__amr_min_identity=90.0, klebsiella_pneumo_complex__amr_min_coverage=80.0, klebsiella_pneumo_complex__amr_min_spurious_identity=80.0, klebsiella_pneumo_complex__amr_min_spurious_coverage=40.0), {})
    assert results['Omp_mutations'] == 'OmpK35-10%'


def test_get_results_4():
    Args = collections.namedtuple('Args', ['klebsiella_pneumo_complex__amr_min_identity', 'klebsiella_pneumo_complex__amr_min_coverage','klebsiella_pneumo_complex__amr_min_spurious_identity', 'klebsiella_pneumo_complex__amr_min_spurious_coverage'])
    results = get_results(get_test_genome_dir() / 'test_res_omp_4.fasta', None,
                          Args(klebsiella_pneumo_complex__amr_min_identity=90.0, klebsiella_pneumo_complex__amr_min_coverage=80.0,klebsiella_pneumo_complex__amr_min_spurious_identity=80.0, klebsiella_pneumo_complex__amr_min_spurious_coverage=40.0), {})
    assert results['Omp_mutations'] == '-'

def test_get_results_5():
    Args = collections.namedtuple('Args', ['klebsiella_pneumo_complex__amr_min_identity', 'klebsiella_pneumo_complex__amr_min_coverage','klebsiella_pneumo_complex__amr_min_spurious_identity', 'klebsiella_pneumo_complex__amr_min_spurious_coverage'])
    results = get_results(get_test_genome_dir() / 'test_res_omp_5.fasta', None,
                          Args(klebsiella_pneumo_complex__amr_min_identity=90.0, klebsiella_pneumo_complex__amr_min_coverage=80.0, klebsiella_pneumo_complex__amr_min_spurious_identity=80.0, klebsiella_pneumo_complex__amr_min_spurious_coverage=40.0), {})
    assert results['Omp_mutations'] == 'OmpK36GD'

def test_get_results_6():
    Args = collections.namedtuple('Args', ['klebsiella_pneumo_complex__amr_min_identity', 'klebsiella_pneumo_complex__amr_min_coverage','klebsiella_pneumo_complex__amr_min_spurious_identity', 'klebsiella_pneumo_complex__amr_min_spurious_coverage'])
    results = get_results(get_test_genome_dir() / 'test_res_omp_6.fasta', None,
                          Args(klebsiella_pneumo_complex__amr_min_identity=90.0, klebsiella_pneumo_complex__amr_min_coverage=80.0, klebsiella_pneumo_complex__amr_min_spurious_identity=80.0, klebsiella_pneumo_complex__amr_min_spurious_coverage=40.0), {})
    assert results['Omp_mutations'] == 'OmpK36TD'

