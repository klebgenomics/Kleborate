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

import collections
import pathlib
import pytest
import re
import tempfile

import kleborate.__main__


def test_get_version():
    assert re.match(r'\d+\.\d+\.\d+', kleborate.__main__.get_version())


def test_get_all_module_names():
    all_module_names = kleborate.__main__.get_all_module_names()
    assert 'general__contig_stats' in all_module_names
    assert 'template' not in all_module_names


def test_get_headers():
    _, modules = kleborate.__main__.import_modules()
    # top_headers, full_headers, stdout_headers = \
    full_headers, stdout_headers = \
        kleborate.__main__.get_headers(['general__contig_stats'], modules)
    # assert top_headers[0] == ''
    # assert top_headers[1] == 'general__contig_stats'
    assert full_headers[0] == 'strain'
    assert stdout_headers[0] == 'strain'
    assert all(h in full_headers for h in stdout_headers)


# def test_check_modules():
#     _, modules = kleborate.__main__.import_modules()
#     Args = collections.namedtuple('Args', ['assemblies'])
#     kleborate.__main__.check_modules(Args(assemblies=['test/test_main/test.fasta']),
#                                      modules, ['general__contig_stats'])

def test_check_modules():
    _, modules = kleborate.__main__.import_modules()
    Args = collections.namedtuple('Args', ['assemblies'])
    args = Args(assemblies=['test/test_main/test.fasta'])
    
    preset_check_modules = []  
    preset_pass_modules = []  
    
    # Calling the check_modules function with all required arguments
    kleborate.__main__.check_modules(args, modules, ['general__contig_stats'], preset_check_modules, preset_pass_modules)



def test_check_assemblies_1():
    Args = collections.namedtuple('Args', ['assemblies'])
    with pytest.raises(SystemExit) as e:
        kleborate.__main__.check_assemblies(Args(assemblies=['test/test_main']))
    assert 'is a directory' in str(e.value)


def test_check_assemblies_2():
    Args = collections.namedtuple('Args', ['assemblies'])
    with pytest.raises(SystemExit) as e:
        kleborate.__main__.check_assemblies(Args(assemblies=['test/test_main/does_not_exist']))
    assert 'could not find' in str(e.value)


def test_check_assemblies_3():
    Args = collections.namedtuple('Args', ['assemblies'])
    with pytest.raises(SystemExit) as e:
        kleborate.__main__.check_assemblies(Args(assemblies=['test/test_main/bad_format.fasta']))
    assert 'invalid' in str(e.value)


def test_check_assemblies_4():
    Args = collections.namedtuple('Args', ['assemblies'])
    with pytest.raises(SystemExit) as e:
        kleborate.__main__.check_assemblies(Args(assemblies=['test/test_main/empty_seq.fasta']))
    assert 'zero-length sequence' in str(e.value)


def test_check_assemblies_5():
    Args = collections.namedtuple('Args', ['assemblies'])
    kleborate.__main__.check_assemblies(Args(assemblies=['test/test_main/test.fasta']))
    kleborate.__main__.check_assemblies(Args(assemblies=['test/test_main/test.fasta.gz']))


def test_decompress_file():
    with tempfile.TemporaryDirectory() as tmp_dir:
        in_file = 'test/test_main/test.fasta.gz'
        out_file = pathlib.Path(tmp_dir) / 'temp.fasta'
        ref_file = 'test/test_main/test.fasta'
        kleborate.__main__.decompress_file(in_file, out_file)
        assert open(out_file, 'rt').read() == open(ref_file, 'rt').read()


def test_gunzip_assembly_if_necessary_1():
    with tempfile.TemporaryDirectory() as tmp_dir:
        assembly = 'test/test_main/test.fasta'
        assert kleborate.__main__.gunzip_assembly_if_necessary(assembly, tmp_dir) == assembly


def test_gunzip_assembly_if_necessary_2():
    with tempfile.TemporaryDirectory() as tmp_dir:
        assembly = 'test/test_main/test.fasta.gz'
        assert kleborate.__main__.gunzip_assembly_if_necessary(assembly, tmp_dir) != assembly


def test_output_headers(capfd):
    full_headers = ['header_a', 'header_b', 'header_c']
    stdout_headers = ['header_a', 'header_b']
    with tempfile.TemporaryDirectory() as tmp_dir:
        out_file = pathlib.Path(tmp_dir) / 'out.txt'
        kleborate.__main__.output_headers(full_headers, stdout_headers, out_file)
        out, err = capfd.readouterr()
        assert out == 'header_a\theader_b\n'
        assert open(out_file, 'rt').read() == 'header_a\theader_b\theader_c'


def test_output_results_1():
    full_headers = ['header_a', 'header_b']
    stdout_headers = ['header_a']
    results = {'header_a': 'result_a', 'header_b': 'result_b', 'header_c': 'result_c'}
    with tempfile.TemporaryDirectory() as tmp_dir:
        out_file = pathlib.Path(tmp_dir) / 'out.txt'
        with pytest.raises(SystemExit) as e:
            kleborate.__main__.output_results(full_headers, stdout_headers, out_file, results)
        assert 'not covered by the output headers' in str(e.value)


def test_get_presets():
    presets = kleborate.__main__.get_presets()
    all_module_names = kleborate.__main__.get_all_module_names()
    for preset, modules in presets.items():
        # Ensure that 'check' and 'pass' keys exist
        assert 'check' in modules
        assert 'pass' in modules

        # Check 'check' modules
        for check_module in modules['check']:
            assert check_module[0] in all_module_names  # Only check the module name, not the condition

        # Check 'pass' modules
        for pass_module in modules['pass']:
            assert pass_module in all_module_names

def test_get_used_module_names_1():
    all_module_names = ['a', 'b', 'c', 'd', 'e']
    presets = {'1': {'check': [['a'], ['b'], ['c']], 'pass': []}, '2': {'check': [['c'], ['d'], ['e']], 'pass': []}}
    Args = collections.namedtuple('Args', ['modules', 'preset'])
    module_names, check_modules, pass_modules = kleborate.__main__.get_used_module_names(Args(modules='b,c,d', preset=None),
                                                                                         all_module_names, presets)
    assert module_names == ['b', 'c', 'd']
    assert check_modules == []
    assert pass_modules == []


# def test_get_used_module_names_1():
#     all_module_names = ['a', 'b', 'c', 'd', 'e']
#     presets = {'1': ['a', 'b', 'c'], '2': ['c', 'd', 'e']}
#     Args = collections.namedtuple('Args', ['modules', 'preset'])
#     modules = kleborate.__main__.get_used_module_names(Args(modules='b,c,d', preset=None),
#                                                        all_module_names, presets)
#     assert modules == ['b', 'c', 'd']


def test_get_used_module_names_2():
    all_module_names = ['a', 'b', 'c', 'd', 'e']
    presets = {'1': {'check': [['a'], ['b'], ['c']], 'pass': []}, '2': {'check': [['c'], ['d'], ['e']], 'pass': []}}
    Args = collections.namedtuple('Args', ['modules', 'preset'])
    module_names, check_modules, pass_modules = kleborate.__main__.get_used_module_names(Args(modules='c,b,a', preset=None),
                                                                                         all_module_names, presets)
    assert module_names == ['c', 'b', 'a']
    assert check_modules == []
    assert pass_modules == []

def test_get_used_module_names_3():
    all_module_names = ['a', 'b', 'c', 'd', 'e']
    presets = {'1': {'check': [['a'], ['b'], ['c']], 'pass': []}, '2': {'check': [['c'], ['d'], ['e']], 'pass': []}}
    Args = collections.namedtuple('Args', ['modules', 'preset'])
    module_names, check_modules, pass_modules = kleborate.__main__.get_used_module_names(Args(modules=None, preset='2'),
                                                                                         all_module_names, presets)
    assert module_names == ['c', 'd', 'e']
    assert check_modules == ['c', 'd', 'e']
    assert pass_modules == []

# def test_get_used_module_names_4():
#     all_module_names = ['a', 'b', 'c', 'd', 'e']
#     presets = {'1': ['a', 'b', 'c'], '2': ['c', 'd', 'e']}
#     Args = collections.namedtuple('Args', ['modules', 'preset'])
#     modules = kleborate.__main__.get_used_module_names(Args(modules='b,c', preset='2'),
#                                                        all_module_names, presets)
#     assert modules == ['c', 'd', 'e', 'b']

def test_get_used_module_names_4():
    all_module_names = ['a', 'b', 'c', 'd', 'e']
    presets = {'1': {'check': [['a'], ['b'], ['c']], 'pass': []}, '2': {'check': [['c'], ['d'], ['e']], 'pass': []}}
    Args = collections.namedtuple('Args', ['modules', 'preset'])
    module_names, check_modules, pass_modules = kleborate.__main__.get_used_module_names(Args(modules='b,c', preset='2'),
                                                                       all_module_names, presets)
    assert module_names == ['c', 'd', 'e', 'b']  
    assert check_modules == ['c', 'd', 'e']
    assert pass_modules == []  



def test_get_used_module_names_5():
    all_module_names = ['a', 'b', 'c', 'd', 'e']
    presets = {'1': ['a', 'b', 'c'], '2': ['c', 'd', 'e']}
    Args = collections.namedtuple('Args', ['modules', 'preset'])
    with pytest.raises(SystemExit) as e:
        kleborate.__main__.get_used_module_names(Args(modules=None, preset='3'),
                                                 all_module_names, presets)
    assert '3 is not a valid preset' in str(e.value)


def test_get_used_module_names_6():
    all_module_names = ['a', 'b', 'c', 'd', 'e']
    presets = {'1': ['a', 'b', 'c'], '2': ['c', 'd', 'e']}
    Args = collections.namedtuple('Args', ['modules', 'preset'])
    with pytest.raises(SystemExit) as e:
        kleborate.__main__.get_used_module_names(Args(modules='a,b,f', preset=None),
                                                 all_module_names, presets)
    assert 'f is not a valid module name' in str(e.value)


def test_get_used_module_names_7():
    all_module_names = ['a', 'b', 'c', 'd', 'e']
    presets = {'1': ['a', 'b', 'c'], '2': ['c', 'd', 'e']}
    Args = collections.namedtuple('Args', ['modules', 'preset'])
    with pytest.raises(SystemExit) as e:
        kleborate.__main__.get_used_module_names(Args(modules=None, preset=None),
                                                 all_module_names, presets)
    assert 'either --preset or --modules is required' in str(e.value)


def test_paper_refs():
    papers = kleborate.__main__.paper_refs()
    assert 'Lam MMC, et al.' in papers
    assert 'Wyres KL, et al.' in papers
    assert '\n\n' in papers


def remove_formatting(text):
    return re.sub('\033.*?m', '', text)


def test_parse_arguments_1(capsys):
    all_module_names, modules = kleborate.__main__.import_modules()
    with pytest.raises(SystemExit):
        kleborate.__main__.parse_arguments(['--help'], all_module_names, modules)
    out, err = capsys.readouterr()
    out = remove_formatting(out)
    assert 'Kleborate:' in out
    assert 'Input/output:' in out
    assert 'enterobacterales__species module:' not in out
    


def test_parse_arguments_2(capsys):
    all_module_names, modules = kleborate.__main__.import_modules()
    with pytest.raises(SystemExit):
        kleborate.__main__.parse_arguments(['--help_all'], all_module_names, modules)
    out, err = capsys.readouterr()
    out = remove_formatting(out)
    print(out)
    assert 'Kleborate:' in out
    assert 'Input/output:' in out
    assert 'enterobacterales__species module:' in out


def test_parse_arguments_3(capsys):
    all_module_names, modules = kleborate.__main__.import_modules()
    with pytest.raises(SystemExit):
        kleborate.__main__.parse_arguments(['--helpall'], all_module_names, modules)
    out, err = capsys.readouterr()
    out = remove_formatting(out)
    print(out)
    assert 'Kleborate:' in out
    assert 'Input/output:' in out
    assert 'enterobacterales__species module:' in out


def test_parse_arguments_4(capsys):
    all_module_names, modules = kleborate.__main__.import_modules()
    with pytest.raises(SystemExit):
        kleborate.__main__.parse_arguments([], all_module_names, modules)
    out, err = capsys.readouterr()
    err = remove_formatting(err)
    assert 'Kleborate:' in err
    assert 'Input/output:' in err
    assert 'enterobacterales__species module:' not in err


def test_get_run_order_1():
    dependency_graph = {'a': [], 'b': [], 'c': []}
    assert kleborate.__main__.get_run_order(dependency_graph) == ['a', 'b', 'c']

    dependency_graph = {'b': [], 'c': [], 'a': []}
    assert kleborate.__main__.get_run_order(dependency_graph) == ['b', 'c', 'a']

    dependency_graph = {'c': [], 'a': [], 'b': []}
    assert kleborate.__main__.get_run_order(dependency_graph) == ['c', 'a', 'b']


def test_get_run_order_2():
    dependency_graph = {'a': ['b'], 'b': [], 'c': []}
    assert sorted(kleborate.__main__.get_run_order(dependency_graph)) == ['a', 'b', 'c']

    dependency_graph = {'a': [], 'b': ['a'], 'c': []}
    assert sorted(kleborate.__main__.get_run_order(dependency_graph)) == ['a', 'b', 'c']

    dependency_graph = {'a': ['b'], 'b': ['a'], 'c': []}
    with pytest.raises(SystemExit) as e:
        assert kleborate.__main__.get_run_order(dependency_graph)
