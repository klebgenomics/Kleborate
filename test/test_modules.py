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

from kleborate.__main__ import import_modules


def test_descriptions():
    # Makes sure that a description is provided for each module.
    module_names, modules = import_modules()
    for module_name in module_names:
        description = modules[module_name].description()
        assert len(description) > 0


def test_header_subset():
     # For each module, tests that stdout_headers are a subset of full_headers.
     module_names, modules = import_modules()
     for module_name in module_names:
         full_headers, stdout_headers = modules[module_name].get_headers()
         assert all(h in full_headers for h in stdout_headers)

def test_no_header_duplicates():
    # No duplicate header names are allowed within a module.
    module_names, modules = import_modules()
    for module_name in module_names:
        full_headers, _ = modules[module_name].get_headers()
        assert len(full_headers) == len(set(full_headers))


