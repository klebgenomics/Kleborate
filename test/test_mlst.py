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

import pathlib
import tempfile

from kleborate.shared.mlst import *
from kleborate.shared.alignment import Alignment


def test_load_st_profiles_1():
    # Tests a profile without a final extra-info column.
    gene_names = ['abcD', 'efgH', 'ijkL']
    with tempfile.TemporaryDirectory() as tmp_dir:
        profile_file = pathlib.Path(tmp_dir) / 'profiles'
        with open(profile_file, 'wt') as f:
            f.write('ST\tabcD\tefgH\tijkL\n')
            f.write('1\t1\t1\t1\n')
            f.write('2\t1\t2\t1\n')
            f.write('3\t1\t1\t2\n')
        profiles = load_st_profiles(profile_file, gene_names, None)
    assert profiles == [(1, [1, 1, 1], None), (2, [1, 2, 1], None), (3, [1, 1, 2], None)]


def test_load_st_profiles_2():
    # Tests a profile with a final extra-info column.
    gene_names = ['abcD', 'efgH', 'ijkL']
    with tempfile.TemporaryDirectory() as tmp_dir:
        profile_file = pathlib.Path(tmp_dir) / 'profiles'
        with open(profile_file, 'wt') as f:
            f.write('ST\tabcD\tefgH\tijkL\textra\n')
            f.write('1\t1\t1\t1\t123\n')
            f.write('2\t1\t2\t1\txyz\n')
            f.write('3\t1\t1\t2\t\n')
        profiles = load_st_profiles(profile_file, gene_names, 'extra')
    assert profiles == [(1, [1, 1, 1], '123'), (2, [1, 2, 1], 'xyz'), (3, [1, 1, 2], '')]


def test_get_best_hits_1():
    # In this test, abcD_2 is the only perfect (100% identity, 100% coverage) hit.
    hits = [Alignment('abcD_1\t1000\t0\t1000\t+\t'
                      'tig\t10000\t3000\t4000\t999\t1000\tAS:i:999\tcg:Z:500=1X499='),
            Alignment('abcD_2\t1000\t0\t1000\t+\t'
                      'tig\t10000\t3000\t4000\t1000\t1000\tAS:i:1000\tcg:Z:1000='),
            Alignment('abcD_3\t1000\t0\t1000\t+\t'
                      'tig\t10000\t3000\t4000\t999\t1000\tAS:i:999\tcg:Z:500=1X499=')]
    best_hits = get_best_hits(hits)
    assert len(best_hits) == 1
    assert best_hits[0].query_name == 'abcD_2'


def test_get_best_hits_2():
    # In this test, abcD_2 and abcD_3 both have perfect (100% identity, 100% coverage) hits, but
    # abcD_3 has a higher alignment score.
    hits = [Alignment('abcD_1\t1000\t0\t1000\t+\t'
                      'tig\t10000\t3000\t4000\t999\t1000\tAS:i:999\tcg:Z:500=1X499='),
            Alignment('abcD_2\t1000\t0\t1000\t+\t'
                      'tig\t10000\t3000\t4000\t1000\t1000\tAS:i:1000\tcg:Z:1000='),
            Alignment('abcD_3\t1100\t0\t1100\t+\t'
                      'tig\t10000\t3000\t4000\t1100\t1100\tAS:i:1100\tcg:Z:1100=')]
    best_hits = get_best_hits(hits)
    assert len(best_hits) == 1
    assert best_hits[0].query_name == 'abcD_3'


def test_get_best_hits_3():
    # Tests a tie between three equally-good hits.
    hits = [Alignment('abcD_1\t1000\t0\t1000\t+\t'
                      'tig\t10000\t3000\t4000\t999\t1000\tAS:i:1000\tcg:Z:1000='),
            Alignment('abcD_2\t1000\t0\t1000\t+\t'
                      'tig\t10000\t3000\t4000\t999\t1000\tAS:i:1000\tcg:Z:1000='),
            Alignment('abcD_3\t1000\t0\t1000\t+\t'
                      'tig\t10000\t3000\t4000\t999\t1000\tAS:i:1000\tcg:Z:1000=')]
    best_hits = get_best_hits(hits)
    assert len(best_hits) == 3


def test_get_best_hits_4():
    # Tests an empty list.
    assert get_best_hits([]) == []


def test_number_from_hit_1():
    hit = Alignment('abcD_2\t1000\t1\t1000\t+\t'
                    'tig\t10000\t3000\t4000\t1000\t1000\tAS:i:1000\tcg:Z:1000=')
    assert number_from_hit(hit) == 2


def test_number_from_hit_2():
    hit = Alignment('abcD_323\t1000\t1\t1000\t+\t'
                    'tig\t10000\t3000\t4000\t1000\t1000\tAS:i:1000\tcg:Z:1000=')
    assert number_from_hit(hit) == 323


def test_number_from_hit_3():
    hit = Alignment('tonB648\t1000\t1\t1000\t+\t'
                    'tig\t10000\t3000\t4000\t1000\t1000\tAS:i:1000\tcg:Z:1000=')
    assert number_from_hit(hit) == 648


def test_number_from_hit_4():
    hit = Alignment('tonB\t1000\t1\t1000\t+\t'
                    'tig\t10000\t3000\t4000\t1000\t1000\tAS:i:1000\tcg:Z:1000=')
    assert number_from_hit(hit) == 0


def test_number_from_hit_5():
    assert number_from_hit(None) == 0


def test_number_from_hit_6():
    hit = Alignment('abc123_45\t1000\t1\t1000\t+\t'
                    'tig\t10000\t3000\t4000\t1000\t1000\tAS:i:1000\tcg:Z:1000=')
    assert number_from_hit(hit) == 45


def test_get_best_matching_profile_1():
    # Tests a perfect match to ST1.
    profiles = [(1, [1, 1, 1], None), (2, [1, 2, 1], None), (3, [1, 1, 2], None)]
    gene_names = ['abcD', 'efgH', 'ijkL']
    best_hits_per_gene = {'abcD': [Alignment('abcD_1\t100\t0\t100\t+\t'
                                             'tig\t100\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')],
                          'efgH': [Alignment('efgH_1\t100\t0\t100\t+\t'
                                             'tig\t100\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')],
                          'ijkL': [Alignment('ijkL_1\t100\t0\t100\t+\t'
                                             'tig\t100\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')]}
    st, alleles, extra_info = get_best_matching_profile(profiles, gene_names, best_hits_per_gene)
    assert st == 1
    assert alleles == [1, 1, 1]

    best_hit_per_gene = get_best_hit_per_gene(gene_names, best_hits_per_gene, alleles)
    assert best_hit_per_gene['abcD'].query_name == 'abcD_1'
    assert best_hit_per_gene['efgH'].query_name == 'efgH_1'
    assert best_hit_per_gene['ijkL'].query_name == 'ijkL_1'


def test_get_best_matching_profile_2():
    # Tests a perfect match to ST2.
    profiles = [(1, [1, 1, 1], None), (2, [1, 2, 1], None), (3, [1, 1, 2], None)]
    gene_names = ['abcD', 'efgH', 'ijkL']
    best_hits_per_gene = {'abcD': [Alignment('abcD_1\t100\t0\t100\t+\t'
                                             'tig\t100\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')],
                          'efgH': [Alignment('efgH_2\t100\t0\t100\t+\t'
                                             'tig\t100\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')],
                          'ijkL': [Alignment('ijkL_1\t100\t0\t100\t+\t'
                                             'tig\t100\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')]}
    st, alleles, extra_info = get_best_matching_profile(profiles, gene_names, best_hits_per_gene)
    assert st == 2
    assert alleles == [1, 2, 1]

    best_hit_per_gene = get_best_hit_per_gene(gene_names, best_hits_per_gene, alleles)
    assert best_hit_per_gene['abcD'].query_name == 'abcD_1'
    assert best_hit_per_gene['efgH'].query_name == 'efgH_2'
    assert best_hit_per_gene['ijkL'].query_name == 'ijkL_1'


def test_get_best_matching_profile_3():
    # This test doesn't match any ST, but it's closest to ST3.
    profiles = [(1, [1, 1, 1], None), (2, [1, 2, 1], None), (3, [1, 1, 2], None)]
    gene_names = ['abcD', 'efgH', 'ijkL']
    best_hits_per_gene = {'abcD': [Alignment('abcD_1\t100\t0\t100\t+\t'
                                             'tig\t100\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')],
                          'efgH': [Alignment('efgH_3\t100\t0\t100\t+\t'
                                             'tig\t100\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')],
                          'ijkL': [Alignment('ijkL_2\t100\t0\t100\t+\t'
                                             'tig\t100\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')]}
    st, alleles, extra_info = get_best_matching_profile(profiles, gene_names, best_hits_per_gene)
    assert st == 3
    assert alleles == [1, 1, 2]

    best_hit_per_gene = get_best_hit_per_gene(gene_names, best_hits_per_gene, alleles)
    assert best_hit_per_gene['abcD'].query_name == 'abcD_1'
    assert best_hit_per_gene['efgH'].query_name == 'efgH_3'
    assert best_hit_per_gene['ijkL'].query_name == 'ijkL_2'


def test_get_best_matching_profile_4():
    # Test no hits at all.
    profiles = [(1, [1, 1, 1], None), (2, [1, 2, 1], None), (3, [1, 1, 2], None)]
    gene_names = ['abcD', 'efgH', 'ijkL']
    best_hits_per_gene = {'abcD': [], 'efgH': [], 'ijkL': []}
    st, alleles, extra_info = get_best_matching_profile(profiles, gene_names, best_hits_per_gene)
    assert st == 0
    assert alleles == [0, 0, 0]

    best_hit_per_gene = get_best_hit_per_gene(gene_names, best_hits_per_gene, alleles)
    assert best_hit_per_gene['abcD'] is None
    assert best_hit_per_gene['efgH'] is None
    assert best_hit_per_gene['ijkL'] is None


def test_get_best_matching_profile_5():
    # This test has two hits for efgH, so it matches equally well to ST1 and ST2. The correct
    # answer is ST1 because that's earlier in the profile.
    profiles = [(1, [1, 1, 1], None), (2, [1, 2, 1], None), (3, [1, 1, 2], None)]
    gene_names = ['abcD', 'efgH', 'ijkL']
    best_hits_per_gene = {'abcD': [Alignment('abcD_1\t100\t0\t100\t+\t'
                                             'tig\t100\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')],
                          'efgH': [Alignment('efgH_2\t100\t0\t100\t+\t'
                                             'tig\t100\t0\t100\t99\t100\tAS:i:99\tcg:Z:99='),
                                   Alignment('efgH_1\t100\t0\t100\t+\t'
                                             'tig\t100\t0\t100\t99\t100\tAS:i:99\tcg:Z:99=')],
                          'ijkL': [Alignment('ijkL_1\t100\t0\t100\t+\t'
                                             'tig\t100\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')]}
    st, alleles, extra_info = get_best_matching_profile(profiles, gene_names, best_hits_per_gene)
    assert st == 1
    assert alleles == [1, 1, 1]

    best_hit_per_gene = get_best_hit_per_gene(gene_names, best_hits_per_gene, alleles)
    assert best_hit_per_gene['abcD'].query_name == 'abcD_1'
    assert best_hit_per_gene['efgH'].query_name == 'efgH_1'
    assert best_hit_per_gene['ijkL'].query_name == 'ijkL_1'


def test_get_best_matching_profile_6():
    # This test has two hits for efgH and two hits for ijkL, so it matches equally well to ST2 and
    # ST3. The correct answer is ST2 because that's earlier in the profile.
    profiles = [(1, [1, 1, 1], None), (2, [1, 2, 1], None), (3, [1, 1, 2], None)]
    gene_names = ['abcD', 'efgH', 'ijkL']
    best_hits_per_gene = {'abcD': [Alignment('abcD_1\t100\t0\t100\t+\t'
                                             'tig\t100\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')],
                          'efgH': [Alignment('efgH_2\t100\t0\t100\t+\t'
                                             'tig\t100\t0\t100\t99\t100\tAS:i:99\tcg:Z:99='),
                                   Alignment('efgH_3\t100\t0\t100\t+\t'
                                             'tig\t100\t0\t100\t99\t100\tAS:i:99\tcg:Z:99=')],
                          'ijkL': [Alignment('ijkL_2\t100\t0\t100\t+\t'
                                             'tig\t100\t0\t100\t99\t100\tAS:i:99\tcg:Z:99='),
                                   Alignment('ijkL_3\t100\t0\t100\t+\t'
                                             'tig\t100\t0\t100\t99\t100\tAS:i:99\tcg:Z:99=')]}
    st, alleles, extra_info = get_best_matching_profile(profiles, gene_names, best_hits_per_gene)
    assert st == 2
    assert alleles == [1, 2, 1]

    best_hit_per_gene = get_best_hit_per_gene(gene_names, best_hits_per_gene, alleles)
    assert best_hit_per_gene['abcD'].query_name == 'abcD_1'
    assert best_hit_per_gene['efgH'].query_name == 'efgH_2'
    assert best_hit_per_gene['ijkL'].query_name == 'ijkL_2'
