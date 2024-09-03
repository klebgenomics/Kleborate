"""
This file contains tests for Kleborate. To run all tests, go the repo's root directory and run:
  python3 -m pytest

To get code coverage stats:
  coverage run --source . -m pytest && coverage report -m

Copyright 2023 Kat Holt, Ryan Wick (rrwick@gmail.com), Mary Maranga (gathonimaranga@gmail.com)
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

from .general__contig_stats import *


def get_file_dir():
    # Returns the path of the directory with the files for these tests.
    return pathlib.Path(__file__).parents[0] / 'test_files'


def test_prerequisite_modules():
    assert prerequisite_modules() == ['enterobacterales__species']


def test_get_headers():
    # stdout_headers must be a subset of full_headers.
    full_headers, stdout_headers = get_headers()
    assert all(h in full_headers for h in stdout_headers)


def test_empty_functions():
    # Tests the functions that aren't used in this module.
    assert add_cli_options(None) is None
    assert check_cli_options(None) is None
    assert check_external_programs() == []


def test_count_1():
    contig_count, _, _, _, _ = get_contig_stats(get_file_dir() / 'test_1.fasta')
    assert contig_count == 4


def test_count_2():
    contig_count, _, _, _, _ = get_contig_stats(get_file_dir() / 'test_2.fasta')
    assert contig_count == 3


def test_n50_1():
    _, N50, _, _, _ = get_contig_stats(get_file_dir() / 'test_1.fasta')
    assert N50 == 40


def test_n50_2():
    _, N50, _, _, _ = get_contig_stats(get_file_dir() / 'test_2.fasta')
    assert N50 == 200


def test_longest_1():
    _, _, longest_contig, _, _ = get_contig_stats(get_file_dir() / 'test_1.fasta')
    assert longest_contig == 45


def test_longest_2():
    _, _, longest_contig, _, _ = get_contig_stats(get_file_dir() / 'test_2.fasta')
    assert longest_contig == 200


def test_ambiguous_bases_1():
    _, _, _, _, ambiguous = get_contig_stats(get_file_dir() / 'test_1.fasta')
    assert ambiguous == 'no'


def test_ambiguous_bases_2():
    _, _, _, _, ambiguous = get_contig_stats(get_file_dir() / 'test_2.fasta')
    assert ambiguous == 'yes (1)'


def test_ambiguous_bases_3():
    _, _, _, _, ambiguous = get_contig_stats(get_file_dir() / 'test_3.fasta')
    assert ambiguous == 'no'


def test_ambiguous_bases_4():
    _, _, _, _, ambiguous = get_contig_stats(get_file_dir() / 'test_4.fasta')
    assert ambiguous == 'yes (4)'


def test_total_size_1():
    _, _, _, total_size, _ = get_contig_stats(get_file_dir() / 'test_1.fasta')
    assert total_size == 115


def test_total_size_2():
    _, _, _, total_size, _ = get_contig_stats(get_file_dir() / 'test_2.fasta')
    assert total_size == 260


def test_total_size_3():
    _, _, _, total_size, _ = get_contig_stats(get_file_dir() / 'test_3.fasta')
    assert total_size == 260


def test_total_size_4():
    _, _, _, total_size, _ = get_contig_stats(get_file_dir() / 'test_4.fasta')
    assert total_size == 260


def test_qc_warnings_1():
    # Define the species specifications for the test
    species_specification_dict = {
        'Klebsiella pneumoniae': {
            'min_genome_size': 5000000,
            'max_genome_size': 6500000
        }
    }
    
    previous_results = {
        'enterobacterales__species__species': 'Klebsiella pneumoniae'
    }
    
    # A perfectly nice assembly - yields no warnings.
    warnings = get_qc_warnings(
        total_size=5000000,  # Total size within acceptable range
        N50=20000,        
        ambiguous_bases='no',  # No ambiguous bases
        species=previous_results['enterobacterales__species__species'],
        species_specification_dict=species_specification_dict
    )
    
    assert warnings == '-'

def test_qc_warnings_2():
    # Define the species specifications for the test
    species_specification_dict = {
        'Klebsiella pneumoniae': {
            'min_genome_size': 5000000,
            'max_genome_size': 6500000
        }
    }
    
    previous_results = {
        'enterobacterales__species__species': 'Klebsiella pneumoniae'
    }
    
    # Small N50 test case with all required arguments
    warnings = get_qc_warnings(
        total_size=5000000,  
        N50=500,          
        ambiguous_bases='no', 
        species=previous_results['enterobacterales__species__species'],
        species_specification_dict=species_specification_dict
    )
    
    assert warnings == 'N50'


def test_empty_file_1():
    previous_results = {
        'enterobacterales__species__species': 'Klebsiella pneumoniae'
    }
    contig_count, N50, longest_contig, total_size, ambiguous = \
    get_contig_stats(get_file_dir() / 'empty.fasta')
    assert contig_count == 0
    assert N50 == 0
    assert longest_contig == 0
    assert total_size == 0
    assert ambiguous == 'no'


def test_get_results():
    # Final results are all in string format.
    previous_results = {
        'enterobacterales__species__species': 'Klebsiella pneumoniae'
    }
    # results = get_results(get_file_dir() / 'test_1.fasta', None, None, {})
    results = get_results(get_file_dir() / 'test_1.fasta', None, None, previous_results)
    assert results['contig_count'] == '4'
    assert results['N50'] == '40'
    assert results['largest_contig'] == '45'
    assert results['total_size'] == '115'
    assert results['ambiguous_bases'] == 'no'
    assert results['QC_warnings'] == 'total_size,N50'
