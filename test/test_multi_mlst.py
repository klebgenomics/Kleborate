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


from kleborate.shared.multi_mlst import *
from kleborate.shared.alignment import Alignment


def test_cluster_hits_by_contig_1():
    gene_names = ['abcD', 'efgH', 'ijkL']
    hits_per_gene = {'abcD': [], 'efgH': [], 'ijkL': []}
    hits_by_contig = cluster_hits_by_contig(hits_per_gene, gene_names)
    assert len(hits_by_contig) == 0


def test_cluster_hits_by_contig_2():
    gene_names = ['abcD', 'efgH', 'ijkL']
    hits_per_gene = {'abcD': [Alignment('abcD_1\t100\t0\t100\t+\t'
                                        'tig_1\t100\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')],
                     'efgH': [Alignment('efgH_1\t100\t0\t100\t+\t'
                                        'tig_1\t100\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')],
                     'ijkL': [Alignment('ijkL_1\t100\t0\t100\t+\t'
                                        'tig_1\t100\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')]}
    hits_by_contig = cluster_hits_by_contig(hits_per_gene, gene_names)
    assert len(hits_by_contig) == 1
    assert len(hits_by_contig['tig_1']['abcD']) == 1
    assert len(hits_by_contig['tig_1']['efgH']) == 1
    assert len(hits_by_contig['tig_1']['ijkL']) == 1


def test_cluster_hits_by_contig_3():
    gene_names = ['abcD', 'efgH', 'ijkL']
    hits_per_gene = {'abcD': [Alignment('abcD_1\t100\t0\t100\t+\t'
                                        'tig_1\t100\t0\t100\t100\t100\tAS:i:100\tcg:Z:100='),
                              Alignment('abcD_2\t100\t0\t100\t+\t'
                                        'tig_1\t100\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')],
                     'efgH': [Alignment('efgH_1\t100\t0\t100\t+\t'
                                        'tig_1\t100\t0\t100\t100\t100\tAS:i:100\tcg:Z:100='),
                              Alignment('efgH_2\t100\t0\t100\t+\t'
                                        'tig_2\t100\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')],
                     'ijkL': [Alignment('ijkL_1\t100\t0\t100\t+\t'
                                        'tig_1\t100\t0\t100\t100\t100\tAS:i:100\tcg:Z:100='),
                              Alignment('ijkL_2\t100\t0\t100\t+\t'
                                        'tig_3\t100\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')]}
    hits_by_contig = cluster_hits_by_contig(hits_per_gene, gene_names)
    assert len(hits_by_contig) == 3
    assert len(hits_by_contig['tig_1']['abcD']) == 2
    assert len(hits_by_contig['tig_1']['efgH']) == 1
    assert len(hits_by_contig['tig_1']['ijkL']) == 1
    assert len(hits_by_contig['tig_2']['abcD']) == 0
    assert len(hits_by_contig['tig_2']['efgH']) == 1
    assert len(hits_by_contig['tig_2']['ijkL']) == 0
    assert len(hits_by_contig['tig_3']['abcD']) == 0
    assert len(hits_by_contig['tig_3']['efgH']) == 0
    assert len(hits_by_contig['tig_3']['ijkL']) == 1


def test_find_full_set_contigs_1():
    hits_by_contig = {'tig_1': {'abcD': [], 'efgH': [], 'ijkL': []},
                      'tig_2': {'abcD': [], 'efgH': [], 'ijkL': []},
                      'tig_3': {'abcD': [], 'efgH': [], 'ijkL': []}}
    assert find_full_set_contigs(hits_by_contig) == []


def test_find_full_set_contigs_2():
    abcd_1_hit = Alignment('abcD_1\t100\t0\t100\t+\t'
                           'tig_1\t100\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')
    efgh_1_hit = Alignment('efgH_1\t100\t0\t100\t+\t'
                           'tig_1\t100\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')
    ijkl_1_hit = Alignment('ijkL_1\t100\t0\t100\t+\t'
                           'tig_1\t100\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')
    hits_by_contig = {'tig_1': {'abcD': [abcd_1_hit], 'efgH': [], 'ijkL': []},
                      'tig_2': {'abcD': [], 'efgH': [efgh_1_hit], 'ijkL': []},
                      'tig_3': {'abcD': [], 'efgH': [], 'ijkL': [ijkl_1_hit]}}
    assert find_full_set_contigs(hits_by_contig) == []


def test_find_full_set_contigs_3():
    abcd_1_hit = Alignment('abcD_1\t100\t0\t100\t+\t'
                           'tig_1\t100\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')
    abcd_2_hit = Alignment('abcD_2\t100\t0\t100\t+\t'
                           'tig_2\t1000\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')
    efgh_1_hit = Alignment('efgH_1\t100\t0\t100\t+\t'
                           'tig_1\t100\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')
    efgh_2_hit = Alignment('efgH_2\t100\t0\t100\t+\t'
                           'tig_2\t1000\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')
    efgh_3_hit = Alignment('efgH_3\t100\t0\t100\t+\t'
                           'tig_3\t10000\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')
    ijkl_1_hit = Alignment('ijkL_1\t100\t0\t100\t+\t'
                           'tig_1\t100\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')
    ijkl_3_hit = Alignment('ijkL_3\t100\t0\t100\t+\t'
                           'tig_3\t10000\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')
    hits_by_contig = {'tig_1': {'abcD': [abcd_1_hit], 'efgH': [efgh_1_hit], 'ijkL': [ijkl_1_hit]},
                      'tig_2': {'abcD': [abcd_2_hit], 'efgH': [efgh_2_hit], 'ijkL': []},
                      'tig_3': {'abcD': [], 'efgH': [efgh_3_hit], 'ijkL': [ijkl_3_hit]}}
    assert find_full_set_contigs(hits_by_contig) == ['tig_1']


def test_find_full_set_contigs_4():
    abcd_1_hit = Alignment('abcD_1\t100\t0\t100\t+\t'
                           'tig_1\t100\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')
    abcd_2_hit = Alignment('abcD_2\t100\t0\t100\t+\t'
                           'tig_2\t1000\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')
    abcd_3_hit = Alignment('abcD_3\t100\t0\t100\t+\t'
                           'tig_3\t10\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')
    efgh_1_hit = Alignment('efgH_1\t100\t0\t100\t+\t'
                           'tig_1\t100\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')
    efgh_2_hit = Alignment('efgH_2\t100\t0\t100\t+\t'
                           'tig_2\t1000\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')
    efgh_3_hit = Alignment('efgH_3\t100\t0\t100\t+\t'
                           'tig_3\t10\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')
    ijkl_1_hit = Alignment('ijkL_1\t100\t0\t100\t+\t'
                           'tig_1\t100\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')
    ijkl_2_hit = Alignment('ijkL_2\t100\t0\t100\t+\t'
                           'tig_2\t1000\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')
    ijkl_3_hit = Alignment('ijkL_3\t100\t0\t100\t+\t'
                           'tig_3\t10\t0\t100\t100\t100\tAS:i:100\tcg:Z:100=')
    hits_by_contig = {'tig_1': {'abcD': [abcd_1_hit], 'efgH': [efgh_1_hit], 'ijkL': [ijkl_1_hit]},
                      'tig_2': {'abcD': [abcd_2_hit], 'efgH': [efgh_2_hit], 'ijkL': [ijkl_2_hit]},
                      'tig_3': {'abcD': [abcd_3_hit], 'efgH': [efgh_3_hit], 'ijkL': [ijkl_3_hit]}}
    assert find_full_set_contigs(hits_by_contig) == ['tig_2', 'tig_1', 'tig_3']


def test_combine_results_1():
    full_set_contigs = ['tig_1', 'tig_2']
    gene_names = ['abcD', 'efgH', 'ijkL']
    contig_results = {'tig_1': ('1', 'lineage1', {'abcD': '1', 'efgH': '2', 'ijkL': '3'}),
                      'tig_2': ('2', 'lineage2', {'abcD': '4', 'efgH': '5', 'ijkL': '6'})}
    combined_st, combined_extra_info, combined_alleles = \
        combine_results(full_set_contigs, contig_results, gene_names)
    assert combined_st == '1,2'
    assert combined_extra_info == 'lineage1,lineage2'
    assert combined_alleles == {'abcD': '1,4', 'efgH': '2,5', 'ijkL': '3,6'}


def test_combine_results_2():
    full_set_contigs = ['tig_2', 'tig_1']
    gene_names = ['abcD', 'efgH', 'ijkL']
    contig_results = {'tig_1': ('1', 'lineage1', {'abcD': '1', 'efgH': '2', 'ijkL': '3'}),
                      'tig_2': ('2', 'lineage2', {'abcD': '4', 'efgH': '5', 'ijkL': '6'})}
    combined_st, combined_extra_info, combined_alleles = \
        combine_results(full_set_contigs, contig_results, gene_names)
    assert combined_st == '2,1'
    assert combined_extra_info == 'lineage2,lineage1'
    assert combined_alleles == {'abcD': '4,1', 'efgH': '5,2', 'ijkL': '6,3'}
