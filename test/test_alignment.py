"""
This file contains tests for Kleborate. To run all tests, go the repo's root directory and run:
  python3 -m pytest

To get code coverage stats:
  coverage run --source . -m pytest && coverage report -m

Copyright 2026 Kat Holt
Copyright 2026 Ryan Wick (rrwick@gmail.com)
https://github.com/katholt/Kleborate/

This file is part of Kleborate. Kleborate is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Kleborate is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Kleborate. If
not, see <https://www.gnu.org/licenses/>.
"""

import inspect
import sys
from types import SimpleNamespace
import pytest

from kleborate.shared.alignment import *


def create_alignment(paf_line):
    """
    Helper to construct an Alignment object across both legacy
    (single paf_line argument) and new (rammappy mapping object)
    signatures.
    """
    parts = paf_line.strip().split('\t')
    if len(parts) < 11:
        sys.exit('Error: line was not in PAF format')

    sig = inspect.signature(Alignment.__init__)
    num_params = len([p for p in sig.parameters.values() if p.name != 'self'])

    if num_params >= 3:
        query_name = parts[0]
        try:
            query_length = int(parts[1])
        except ValueError:
            query_length = 0

        # Parse tags (AS:i:..., cg:Z:..., etc.)
        score = None
        cigar = None
        for tag in parts[11:]:
            if tag.startswith('AS:i:'):
                score = int(tag.split(':')[-1])
            elif tag.startswith('cg:Z:'):
                cigar = tag.split(':')[-1].encode()

        strand_val = Strand.Forward if parts[4] == '+' else Strand.Reverse

        mapping = SimpleNamespace(
            query_start=int(parts[2]),
            query_end=int(parts[3]),
            strand=strand_val,
            target_name=parts[5],
            target_len=int(parts[6]),
            target_start=int(parts[7]),
            target_end=int(parts[8]),
            matches=int(parts[9]),
            block_len=int(parts[10]),
            score=score,
            cigar=cigar,
        )
        return Alignment(mapping, query_name, query_length)

    return Alignment(paf_line)


def test_bad_paf():
    with pytest.raises(SystemExit) as e:
        create_alignment('not_a_paf_line')
    assert 'PAF format' in str(e.value)


def test_repr():
    a = create_alignment('A\t1000\t50\t150\t+\tC\t1000\t60\t160\t100\t100\tAS:i:100\tcg:Z:100=')
    assert str(a) == 'A:50-150(+), C:60-160 (100.000%)'


def test_is_exact():
    a = create_alignment('A\t100\t0\t100\t+\tC\t1000\t60\t160\t100\t100\tAS:i:100\tcg:Z:100=')
    assert a.is_exact()

    b = create_alignment('A\t100\t0\t100\t+\tC\t1000\t60\t160\t90\t100\tAS:i:100\tcg:Z:100=')
    assert not b.is_exact()  # identity < 100%

    a = create_alignment('A\t100\t0\t90\t+\tC\t1000\t60\t160\t90\t90\tAS:i:100\tcg:Z:100=')
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