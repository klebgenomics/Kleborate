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

import pytest

from kleborate.shared.alignment import *


def test_bad_paf():
    with pytest.raises(SystemExit) as e:
        Alignment('not_a_paf_line')
    assert 'PAF format' in str(e.value)


def test_repr():
    a = Alignment('A\t1000\t50\t150\t+\tC\t1000\t60\t160\t100\t100\tAS:i:100\tcg:Z:100=')
    assert str(a) == 'A:50-150(+), C:60-160 (100.000%)'


def test_is_exact():
    a = Alignment('A\t100\t0\t100\t+\tC\t1000\t60\t160\t100\t100\tAS:i:100\tcg:Z:100=')
    assert a.is_exact()

    b = Alignment('A\t100\t0\t100\t+\tC\t1000\t60\t160\t90\t100\tAS:i:100\tcg:Z:100=')
    assert not b.is_exact()  # identity < 100%

    a = Alignment('A\t100\t0\t90\t+\tC\t1000\t60\t160\t90\t90\tAS:i:100\tcg:Z:100=')
    assert not a.is_exact()  # coverage < 100%


def test_get_expanded_cigar():
    assert get_expanded_cigar('5=') == '====='
    assert get_expanded_cigar('3=1I4=2D2=1X4=') == '===I====DD==X===='
    assert get_expanded_cigar('') == ''


def test_sequences_1():
    alignments = align_query_to_ref('test/test_alignment/query.fasta',
                                    'test/test_alignment/forward_hit.fasta')
    assert len(alignments) == 1
    a = alignments[0]
    assert a.strand == '+'
    assert a.percent_identity == pytest.approx(100.0)
    assert a.query_cov == pytest.approx(100.0)
    assert a.ref_cov == pytest.approx(10.0)
    assert len(a.query_seq) == 1000
    assert len(a.ref_seq) == 1000
    assert a.query_seq.startswith('CTTCCACAACCCTCCCAAATGTCCC')
    assert a.ref_seq.startswith('CTTCCACAACCCTCCCAAATGTCCC')
    assert a.query_seq.endswith('ATGCGCGTTAGCTGCCTGACAGCTG')
    assert a.ref_seq.endswith('ATGCGCGTTAGCTGCCTGACAGCTG')


def test_sequences_2():
    alignments = align_query_to_ref('test/test_alignment/query.fasta',
                                    'test/test_alignment/reverse_hit.fasta')
    assert len(alignments) == 1
    a = alignments[0]
    assert a.strand == '-'
    assert a.percent_identity == pytest.approx(100.0)
    assert a.query_cov == pytest.approx(100.0)
    assert a.ref_cov == pytest.approx(10.0)
    assert len(a.query_seq) == 1000
    assert len(a.ref_seq) == 1000
    assert a.query_seq.startswith('CTTCCACAACCCTCCCAAATGTCCC')
    assert a.ref_seq.startswith('CTTCCACAACCCTCCCAAATGTCCC')
    assert a.query_seq.endswith('ATGCGCGTTAGCTGCCTGACAGCTG')
    assert a.ref_seq.endswith('ATGCGCGTTAGCTGCCTGACAGCTG')


def test_sequences_3():
    alignments = align_query_to_ref('test/test_alignment/query.fasta',
                                    'test/test_alignment/imperfect_hit.fasta')
    assert len(alignments) == 1
    a = alignments[0]
    assert a.strand == '+'
    assert a.percent_identity < 100.0
    assert a.query_cov == pytest.approx(100.0)
    assert a.ref_cov > 10.0
    assert len(a.query_seq) == 1000
    assert len(a.ref_seq) == 1001
    assert a.query_seq.startswith('CTTCCACAACCCTCCCAAATGTCCC')
    assert a.ref_seq.startswith('CTTCCACAACCCTCCCAAATGTCCC')
    assert a.query_seq.endswith('ATGCGCGTTAGCTGCCTGACAGCTG')
    assert a.ref_seq.endswith('ATGCGCGTTAGCTGCCTGACAGCTG')
