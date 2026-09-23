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


from kleborate.modules.klebsiella_pneumo_complex__amr.shv_mutations import*
from kleborate.shared.resMinimap import read_class_file, get_res_headers, resminimap_assembly
from kleborate.modules.klebsiella_pneumo_complex__amr.klebsiella_pneumo_complex__amr import get_headers, get_results


def get_test_genome_dir():
    return pathlib.Path(__file__).parents[4] / 'test' / 'test_shv'


def test_get_results_1():
    #This test has an exact match for SHV-1.
    Args = collections.namedtuple('Args', ['klebsiella_pneumo_complex__amr_min_identity', 'klebsiella_pneumo_complex__amr_min_coverage', 'klebsiella_pneumo_complex__amr_min_spurious_identity', 'klebsiella_pneumo_complex__amr_min_spurious_coverage'])
    results = get_results(get_test_genome_dir() / '01.fasta', None,
                          Args(klebsiella_pneumo_complex__amr_min_identity=90.0, klebsiella_pneumo_complex__amr_min_coverage=80.0, klebsiella_pneumo_complex__amr_min_spurious_identity=80.0, klebsiella_pneumo_complex__amr_min_spurious_coverage=40.0), {})
    assert results['Bla_chr'] == 'SHV-1'
    assert results['SHV_mutations'] == '-'

def test_get_results_2():
    """
    This test has a match for SHV-1 with a mutation at site 238 (G -> Y). This changes the
    class to ESBL, so the mutation is included in the
    """
    Args = collections.namedtuple('Args', ['klebsiella_pneumo_complex__amr_min_identity', 'klebsiella_pneumo_complex__amr_min_coverage', 'klebsiella_pneumo_complex__amr_min_spurious_identity', 'klebsiella_pneumo_complex__amr_min_spurious_coverage'])
    results = get_results(get_test_genome_dir() / '02.fasta', None,
                          Args(klebsiella_pneumo_complex__amr_min_identity=90.0, klebsiella_pneumo_complex__amr_min_coverage=80.0, klebsiella_pneumo_complex__amr_min_spurious_identity=80.0, klebsiella_pneumo_complex__amr_min_spurious_coverage=40.0 ), {})
    assert results['Bla_ESBL_acquired'] == 'SHV-1* +238Y'
    assert results['SHV_mutations'] == '238Y'

def test_get_results_3():
    """
    Same as test 2, but the gene is on the reverse strand.
    """
    Args = collections.namedtuple('Args', ['klebsiella_pneumo_complex__amr_min_identity', 'klebsiella_pneumo_complex__amr_min_coverage','klebsiella_pneumo_complex__amr_min_spurious_identity', 'klebsiella_pneumo_complex__amr_min_spurious_coverage'])
    results = get_results(get_test_genome_dir() / '03.fasta', None,
                          Args(klebsiella_pneumo_complex__amr_min_identity=90.0, klebsiella_pneumo_complex__amr_min_coverage=80.0, klebsiella_pneumo_complex__amr_min_spurious_identity=80.0, klebsiella_pneumo_complex__amr_min_spurious_coverage=40.0), {})
    assert results['Bla_ESBL_acquired'] == 'SHV-1* +238Y'
    assert results['SHV_mutations'] == '238Y'

def test_get_results_4():
    """
    This test has a match for SHV-1 with a mutation at site 50 (G -> Y). This doesn't change
    resistance and so won't be reported.
    """
    Args = collections.namedtuple('Args', ['klebsiella_pneumo_complex__amr_min_identity', 'klebsiella_pneumo_complex__amr_min_coverage','klebsiella_pneumo_complex__amr_min_spurious_identity', 'klebsiella_pneumo_complex__amr_min_spurious_coverage'])
    results = get_results(get_test_genome_dir() / '04.fasta', None,
                          Args(klebsiella_pneumo_complex__amr_min_identity=90.0, klebsiella_pneumo_complex__amr_min_coverage=80.0, klebsiella_pneumo_complex__amr_min_spurious_identity=80.0, klebsiella_pneumo_complex__amr_min_spurious_coverage=40.0), {})
    assert results['Bla_chr'] == 'SHV-1*'
    assert results['SHV_mutations'] == '-'

def test_get_results_5():
    """
    This test has an exact match for SHV-29.
    """
    Args = collections.namedtuple('Args', ['klebsiella_pneumo_complex__amr_min_identity', 'klebsiella_pneumo_complex__amr_min_coverage', 'klebsiella_pneumo_complex__amr_min_spurious_identity', 'klebsiella_pneumo_complex__amr_min_spurious_coverage'])
    results = get_results(get_test_genome_dir() / '05.fasta', None,
                          Args(klebsiella_pneumo_complex__amr_min_identity=90.0, klebsiella_pneumo_complex__amr_min_coverage=80.0, klebsiella_pneumo_complex__amr_min_spurious_identity=80.0, klebsiella_pneumo_complex__amr_min_spurious_coverage=40.0), {})
    assert results['Bla_ESBL_acquired'] == 'SHV-29'
    assert results['SHV_mutations'] == '238A;35Q'

def test_get_results_6():
    """
    This test has SHV-29 plus an inhibition mutation.
    """
    Args = collections.namedtuple('Args', ['klebsiella_pneumo_complex__amr_min_identity', 'klebsiella_pneumo_complex__amr_min_coverage','klebsiella_pneumo_complex__amr_min_spurious_identity', 'klebsiella_pneumo_complex__amr_min_spurious_coverage'])
    results = get_results(get_test_genome_dir() / '06.fasta', None,
                          Args(klebsiella_pneumo_complex__amr_min_identity=90.0, klebsiella_pneumo_complex__amr_min_coverage=80.0, klebsiella_pneumo_complex__amr_min_spurious_identity=80.0, klebsiella_pneumo_complex__amr_min_spurious_coverage=40.0), {})
    assert results['Bla_ESBL_inhR_acquired'] == 'SHV-29* +234Y'
    assert results['SHV_mutations'] ==  '234Y;238A;35Q'

def test_get_results_7():
    """
    This test has SHV-1 with position 238 deleted. Since it's not in the omega loop, this isn't
    reported and doesn't have an effect.
    """
    Args = collections.namedtuple('Args', ['klebsiella_pneumo_complex__amr_min_identity', 'klebsiella_pneumo_complex__amr_min_coverage','klebsiella_pneumo_complex__amr_min_spurious_identity', 'klebsiella_pneumo_complex__amr_min_spurious_coverage'])
    results = get_results(get_test_genome_dir() / '07.fasta', None,
                          Args(klebsiella_pneumo_complex__amr_min_identity=90.0, klebsiella_pneumo_complex__amr_min_coverage=80.0,klebsiella_pneumo_complex__amr_min_spurious_identity=80.0, klebsiella_pneumo_complex__amr_min_spurious_coverage=40.0), {})
    assert results['Bla_chr'] == 'SHV-1*'
    assert results['SHV_mutations'] == '-'

def test_get_results_8():
    """
    This test has SHV-1 with a synonymous mutation in the omega loop (so not reported).
    """
    Args = collections.namedtuple('Args', ['klebsiella_pneumo_complex__amr_min_identity', 'klebsiella_pneumo_complex__amr_min_coverage','klebsiella_pneumo_complex__amr_min_spurious_identity', 'klebsiella_pneumo_complex__amr_min_spurious_coverage'])
    results = get_results(get_test_genome_dir() / '08.fasta', None,
                          Args(klebsiella_pneumo_complex__amr_min_identity=90.0, klebsiella_pneumo_complex__amr_min_coverage=80.0, klebsiella_pneumo_complex__amr_min_spurious_identity=80.0, klebsiella_pneumo_complex__amr_min_spurious_coverage=40.0), {})
    assert results['Bla_chr'] == 'SHV-1^'
    assert results['SHV_mutations'] == '-'

def test_get_results_9():
    """
    This test has SHV-1 with a nonsynonymous mutation in the omega loop (so it is reported).
    """
    Args = collections.namedtuple('Args', ['klebsiella_pneumo_complex__amr_min_identity', 'klebsiella_pneumo_complex__amr_min_coverage','klebsiella_pneumo_complex__amr_min_spurious_identity', 'klebsiella_pneumo_complex__amr_min_spurious_coverage'])
    results = get_results(get_test_genome_dir() / '09.fasta', None,
                          Args(klebsiella_pneumo_complex__amr_min_identity=90.0, klebsiella_pneumo_complex__amr_min_coverage=80.0, klebsiella_pneumo_complex__amr_min_spurious_identity=80.0, klebsiella_pneumo_complex__amr_min_spurious_coverage=40.0), {})
    print(results)
    assert results['Bla_chr'] == 'SHV-1*'
    assert results['SHV_mutations'] == '174R;omega-loop=RWETELNEALRGDARD'

def test_bla_class_01():
    assert get_new_bla_class(False, False) == 'Bla_chr'

def test_bla_class_02():
    assert get_new_bla_class(True, False) == 'Bla_ESBL'

def test_bla_class_03():
    assert get_new_bla_class(False, True) =='Bla_inhR'

def test_bla_class_04():
    assert get_new_bla_class(True, True) == 'Bla_ESBL_inhR'

def test_get_class_changing_mutations_01():
    mutations = get_class_changing_mutations('Bla_chr', 'Bla_chr', ['A'], ['B'])
    assert mutations == []

def test_get_class_changing_mutations_02():
    mutations = get_class_changing_mutations('Bla_ESBL', 'Bla_ESBL', ['A'], ['B'])
    assert mutations == []

def test_get_class_changing_mutations_03():
    mutations = get_class_changing_mutations('Bla_chr', 'Bla_ESBL', ['A'], ['B'])
    assert mutations == ['A']

def test_get_class_changing_mutations_04():
    mutations = get_class_changing_mutations('Bla_ESBL', 'Bla_chr', ['A'], ['B'])
    assert mutations == ['A']

def test_get_class_changing_mutations_05():
    mutations = get_class_changing_mutations('Bla_chr', 'Bla_inhR', ['A'], ['B'])
    assert mutations == ['B']

def test_get_class_changing_mutations_06():
    mutations = get_class_changing_mutations('Bla_inhR', 'Bla_chr', ['A'], ['B'])
    assert mutations == ['B']

def test_get_class_changing_mutations_07():
    mutations = get_class_changing_mutations('Bla_chr', 'Bla_ESBL_inhR', ['A'], ['B'])
    assert mutations == ['A', 'B']

def test_get_class_changing_mutations_08():
    mutations = get_class_changing_mutations('Bla_ESBL_inhR', 'Bla_chr', ['A'], ['B'])
    assert mutations ==['A', 'B']





