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

import argparse
import collections
import pathlib
import pytest

from .template import *


def get_file_dir():
    # Returns the path of the directory with the files for these tests.
    return pathlib.Path(__file__).parents[0] / 'test_files'


def test_prerequisite_modules():
    assert prerequisite_modules() == []


def test_get_results():
    # Final results are all in string format.
    results = get_results(get_file_dir() / 'test_file', None, None, {})
    assert results['header_a'] == 'result_a'
    assert results['header_b'] == 'result_b'
    assert results['header_c'] == 'result_c'
    assert sorted(results.keys()) == sorted(get_headers()[0])


def test_check_cli_options_1():
    Args = collections.namedtuple('Args', ['template_opt1', 'template_opt2', 'template_flag'])
    check_cli_options(Args(template_opt1='abc', template_opt2=5, template_flag=True))


def test_check_cli_options_2():
    Args = collections.namedtuple('Args', ['template_opt1', 'template_opt2', 'template_flag'])
    with pytest.raises(SystemExit):
        check_cli_options(Args(template_opt1='abc', template_opt2=15, template_flag=True))


def test_check_external_programs_1(mocker):
    # Tests the good case where both programs are found.
    mocker.patch(
        'shutil.which',
        side_effect=lambda x: {'mash': '/usr/bin/mash', 'minimap2': '/usr/bin/minimap2'}[x],
    )
    assert check_external_programs() == ['mash', 'minimap2']


def test_check_external_programs_2(mocker):
    # Tests the case where mash is missing.
    mocker.patch(
        'shutil.which',
        side_effect=lambda x: {'mash': None, 'minimap2': '/usr/bin/minimap2'}[x],
    )
    with pytest.raises(SystemExit):
        check_external_programs()


def test_check_external_programs_3(mocker):
    # Tests the case where minimap2 is missing.
    mocker.patch(
        'shutil.which',
        side_effect=lambda x: {'mash': '/usr/bin/mash', 'minimap2': None}[x],
    )
    with pytest.raises(SystemExit):
        check_external_programs()


def test_add_cli_options():
    parser = argparse.ArgumentParser()
    add_cli_options(parser)
    assert '--template_opt1' in parser.format_help()
    assert '--template_opt2' in parser.format_help()
    assert '--template_flag' in parser.format_help()
