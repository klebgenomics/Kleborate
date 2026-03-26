"""
This file contains tests for Kleborate. To run all tests, go the repo's root directory and run:
  python3 -m pytest

To get code coverage stats:
  coverage run --source . -m pytest && coverage report -m

Copyright 2025 Kat Holt
https://github.com/klebgenomics/Kleborate/

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

from .acineto__mlst_oxford import *


def get_test_genome_dir():
    return pathlib.Path(__file__).parents[3] / 'test' / 'test_genomes'


def test_prerequisite_modules():
    assert prerequisite_modules() == []


def test_check_cli_options_1():
    Args = collections.namedtuple('Args', ['acinetobacter_mlst_oxford_min_identity', 'acinetobacter_mlst_oxford_min_coverage',
                                           'acinetobacter_mlst_oxford_required_exact_matches'])
    check_cli_options(Args(acinetobacter_mlst_oxford_min_identity=90.0, acinetobacter_mlst_oxford_min_coverage=90.0,
                           acinetobacter_mlst_oxford_required_exact_matches=3))


def test_check_cli_options_2():
    Args = collections.namedtuple('Args', ['acinetobacter_mlst_oxford_min_identity', 'acinetobacter_mlst_oxford_min_coverage',
                                           'acinetobacter_mlst_oxford_required_exact_matches'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(acinetobacter_mlst_oxford_min_identity=0.90, acinetobacter_mlst_oxford_min_coverage=90.0,
                               acinetobacter_mlst_oxford_required_exact_matches=3))


def test_check_cli_options_3():
    Args = collections.namedtuple('Args', ['acinetobacter_mlst_oxford_min_identity', 'acinetobacter_mlst_oxford_min_coverage',
                                           'acinetobacter_mlst_oxford_required_exact_matches'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(acinetobacter_mlst_oxford_min_identity=-90.0, acinetobacter_mlst_oxford_min_coverage=0.90,
                               acinetobacter_mlst_oxford_required_exact_matches=3))


def test_check_cli_options_4():
    Args = collections.namedtuple('Args', ['acinetobacter_mlst_oxford_min_identity', 'acinetobacter_mlst_oxford_min_coverage',
                                           'acinetobacter_mlst_oxford_required_exact_matches'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(acinetobacter_mlst_oxford_min_identity=-10.0, acinetobacter_mlst_oxford_min_coverage=90.0,
                               acinetobacter_mlst_oxford_required_exact_matches=3))


def test_check_cli_options_5():
    Args = collections.namedtuple('Args', ['acinetobacter_mlst_oxford_min_identity', 'acinetobacter_mlst_oxford_min_coverage',
                                           'acinetobacter_mlst_oxford_required_exact_matches'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(acinetobacter_mlst_oxford_min_identity=90.0, acinetobacter_mlst_oxford_min_coverage=120.0,
                               acinetobacter_mlst_oxford_required_exact_matches=3))


def test_check_cli_options_6():
    Args = collections.namedtuple('Args', ['acinetobacter_mlst_oxford_min_identity', 'acinetobacter_mlst_oxford_min_coverage',
                                           'acinetobacter_mlst_oxford_required_exact_matches'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(acinetobacter_mlst_oxford_min_identity=90.0, acinetobacter_mlst_oxford_min_coverage=90.0,
                               acinetobacter_mlst_oxford_required_exact_matches=-2))


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

