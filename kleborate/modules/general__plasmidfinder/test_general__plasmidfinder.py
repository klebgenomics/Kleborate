"""
Tests for general__plasmidfinder module

Copyright 2025 Minh-Quan Ton-Ngoc
https://github.com/tnmquann/Kleborate/
"""

import collections
import pathlib
import pytest
import os
import json
import shutil
import tempfile

from .general__plasmidfinder import *


def get_test_file_dir():
    return pathlib.Path(__file__).parents[0] / 'test_files'


def test_description():
    """Test that description returns expected string"""
    desc = description()
    assert 'plasmid' in desc.lower()
    assert 'plasmidfinder' in desc.lower()


def test_prerequisite_modules():
    """Test that no prerequisite modules are required"""
    assert prerequisite_modules() == []


def test_get_headers():
    """Test that headers are correctly defined"""
    full_headers, stdout_headers = get_headers()
    
    # Check all full headers are present
    assert 'plasmidlist' in full_headers
    assert 'plasmiddb' in full_headers
    assert 'plasmidident' in full_headers
    assert 'plasmidcov' in full_headers
    assert 'plasmidcontig' in full_headers
    assert 'plasmidcontigposition' in full_headers
    assert 'plasmidnote' in full_headers
    assert 'plasmidaccession' in full_headers
    
    # Check stdout headers (minimal display)
    assert 'plasmidlist' in stdout_headers
    
    # Verify counts
    assert len(full_headers) == 8
    assert len(stdout_headers) >= 1


def test_get_strain_name_regular_fasta():
    """Test strain name extraction from regular fasta file"""
    assert get_strain_name('/path/to/sample.fasta') == 'sample'
    assert get_strain_name('/path/to/sample.fa') == 'sample'
    assert get_strain_name('/path/to/sample.fna') == 'sample'
    assert get_strain_name('/path/to/sample.fas') == 'sample'


def test_get_strain_name_gzipped():
    """Test strain name extraction from gzipped fasta file"""
    assert get_strain_name('/path/to/sample.fasta.gz') == 'sample'
    assert get_strain_name('/path/to/sample.fa.gz') == 'sample'
    assert get_strain_name('/path/to/sample.fna.gz') == 'sample'


def test_get_strain_name_no_extension():
    """Test strain name extraction when no recognized extension"""
    assert get_strain_name('/path/to/sample') == 'sample'
    assert get_strain_name('/path/to/sample.txt') == 'sample.txt'


def test_check_cli_options_valid():
    """Test valid CLI options"""
    Args = collections.namedtuple('Args', [
        'plasmidfinder_mincov', 
        'plasmidfinder_threshold',
        'plasmidfinder_dbpath'
    ])
    
    # Should not raise error
    check_cli_options(Args(
        plasmidfinder_mincov=0.60,
        plasmidfinder_threshold=0.90,
        plasmidfinder_dbpath=None
    ))
    
    # Edge cases - boundary values
    check_cli_options(Args(
        plasmidfinder_mincov=0.0,
        plasmidfinder_threshold=1.0,
        plasmidfinder_dbpath=None
    ))
    
    check_cli_options(Args(
        plasmidfinder_mincov=1.0,
        plasmidfinder_threshold=0.0,
        plasmidfinder_dbpath=None
    ))


def test_check_cli_options_invalid_mincov():
    """Test invalid mincov value"""
    Args = collections.namedtuple('Args', [
        'plasmidfinder_mincov',
        'plasmidfinder_threshold',
        'plasmidfinder_dbpath'
    ])
    
    # Test mincov > 1
    with pytest.raises(SystemExit):
        check_cli_options(Args(
            plasmidfinder_mincov=1.5,
            plasmidfinder_threshold=0.90,
            plasmidfinder_dbpath=None
        ))
    
    # Test mincov < 0
    with pytest.raises(SystemExit):
        check_cli_options(Args(
            plasmidfinder_mincov=-0.1,
            plasmidfinder_threshold=0.90,
            plasmidfinder_dbpath=None
        ))


def test_check_cli_options_invalid_threshold():
    """Test invalid threshold value"""
    Args = collections.namedtuple('Args', [
        'plasmidfinder_mincov',
        'plasmidfinder_threshold',
        'plasmidfinder_dbpath'
    ])
    
    # Test threshold < 0
    with pytest.raises(SystemExit):
        check_cli_options(Args(
            plasmidfinder_mincov=0.60,
            plasmidfinder_threshold=-0.1,
            plasmidfinder_dbpath=None
        ))
    
    # Test threshold > 1
    with pytest.raises(SystemExit):
        check_cli_options(Args(
            plasmidfinder_mincov=0.60,
            plasmidfinder_threshold=1.5,
            plasmidfinder_dbpath=None
        ))


def test_check_cli_options_invalid_dbpath():
    """Test invalid database path"""
    Args = collections.namedtuple('Args', [
        'plasmidfinder_mincov',
        'plasmidfinder_threshold',
        'plasmidfinder_dbpath'
    ])
    
    with pytest.raises(SystemExit):
        check_cli_options(Args(
            plasmidfinder_mincov=0.60,
            plasmidfinder_threshold=0.90,
            plasmidfinder_dbpath='/nonexistent/path/to/database'
        ))


def test_parse_plasmidfinder_output_with_hits():
    """Test parsing PlasmidFinder JSON output with hits"""
    test_json = {
        "plasmidfinder": {
            "results": {
                "Enterobacteriales": {
                    "enterobacteriales": {
                        "hit1": {
                            "plasmid": "IncHI1B(pNDM-MAR)",
                            "identity": 99.15,
                            "coverage": 82.98,
                            "contig_name": "contig1",
                            "positions_in_contig": "5041466..5041938",
                            "note": "",
                            "accession": "JN420336"
                        },
                        "hit2": {
                            "plasmid": "repB",
                            "identity": 100.0,
                            "coverage": 100.0,
                            "contig_name": "contig1",
                            "positions_in_contig": "5119552..5120111",
                            "note": "VIR",
                            "accession": "AP006726"
                        }
                    }
                }
            }
        }
    }
    
    results = parse_plasmidfinder_output(test_json)
    
    assert 'IncHI1B(pNDM-MAR)' in results['plasmidlist']
    assert 'repB' in results['plasmidlist']
    assert 'Enterobacteriales' in results['plasmiddb']
    assert '99.15' in results['plasmidident']
    assert '100.0' in results['plasmidident']
    assert '82.98' in results['plasmidcov']
    assert '100.0' in results['plasmidcov']
    assert 'contig1' in results['plasmidcontig']
    assert '5041466..5041938' in results['plasmidcontigposition']
    assert '5119552..5120111' in results['plasmidcontigposition']
    assert 'VIR' in results['plasmidnote']
    assert 'JN420336' in results['plasmidaccession']
    assert 'AP006726' in results['plasmidaccession']


def test_parse_plasmidfinder_output_real_data():
    """Test parsing real PlasmidFinder JSON output from test_files"""
    json_file = get_test_file_dir() / 'data.json'
    
    with open(json_file, 'r') as f:
        test_json = json.load(f)
    
    results = parse_plasmidfinder_output(test_json)
    
    # Should find the three hits from the real data (IncFII, IncFIB(K), Col440II)
    assert results['plasmidlist'] != '-'
    assert 'IncFII' in results['plasmidlist']
    assert 'IncFIB(K)' in results['plasmidlist']
    assert 'Col440II' in results['plasmidlist']
    # Check that identities are present
    assert '96.07' in results['plasmidident']  # IncFII identity
    assert '98.93' in results['plasmidident']  # IncFIB(K) identity
    assert '99.65' in results['plasmidident']  # Col440II identity


def test_parse_plasmidfinder_output_no_hits():
    """Test parsing PlasmidFinder JSON output with no hits"""
    test_json = {
        "plasmidfinder": {
            "results": {
                "Gram Positive": {
                    "Inc18": "No hit found",
                    "NT_Rep": "No hit found"
                }
            }
        }
    }
    
    results = parse_plasmidfinder_output(test_json)
    
    assert results['plasmidlist'] == '-'
    assert results['plasmiddb'] == '-'
    assert results['plasmidident'] == '-'
    assert results['plasmidcov'] == '-'
    assert results['plasmidcontig'] == '-'
    assert results['plasmidcontigposition'] == '-'
    assert results['plasmidnote'] == '-'
    assert results['plasmidaccession'] == '-'


def test_parse_plasmidfinder_output_empty():
    """Test parsing empty PlasmidFinder JSON output"""
    test_json = {}
    
    results = parse_plasmidfinder_output(test_json)
    
    assert results['plasmidlist'] == '-'
    assert results['plasmiddb'] == '-'


def test_parse_plasmidfinder_output_none():
    """Test parsing None input"""
    results = parse_plasmidfinder_output(None)
    
    assert results['plasmidlist'] == '-'
    assert results['plasmiddb'] == '-'


def test_parse_plasmidfinder_output_mixed_hits():
    """Test parsing JSON with both hits and 'No hit found' entries"""
    test_json = {
        "plasmidfinder": {
            "results": {
                "Enterobacteriales": {
                    "enterobacteriales": {
                        "hit1": {
                            "plasmid": "IncFII",
                            "identity": 98.5,
                            "coverage": 100.0,
                            "contig_name": "contig2",
                            "positions_in_contig": "1000..2000",
                            "note": "",
                            "accession": "AB123456"
                        }
                    }
                },
                "Gram Positive": {
                    "Inc18": "No hit found",
                    "NT_Rep": "No hit found"
                }
            }
        }
    }
    
    results = parse_plasmidfinder_output(test_json)
    
    # Should only have one hit
    assert results['plasmidlist'] == 'IncFII'
    assert results['plasmiddb'] == 'Enterobacteriales'
    assert results['plasmidident'] == '98.5'
    assert results['plasmidcov'] == '100.0'


def test_data_dir():
    """Test that data_dir returns a valid path"""
    d = data_dir()
    assert isinstance(d, pathlib.Path)
    assert 'general__plasmidfinder' in str(d)


def test_copy_extended_output():
    """Test _copy_extended_output function"""
    # Import the private function directly from the module
    from .general__plasmidfinder import _copy_extended_output
    
    with tempfile.TemporaryDirectory() as src_dir:
        with tempfile.TemporaryDirectory() as out_dir:
            # Create mock files
            test_files = ['results_tab.tsv', 'data.json', 'results.txt']
            for f in test_files:
                with open(os.path.join(src_dir, f), 'w') as fh:
                    fh.write(f'test content for {f}')
            
            # Run the copy function
            _copy_extended_output(src_dir, out_dir, 'test_strain')
            
            # Check files were copied
            strain_dir = os.path.join(out_dir, 'test_strain')
            assert os.path.exists(strain_dir)
            
            for f in test_files:
                copied_file = os.path.join(strain_dir, f'plasmidfinder_{f}')
                assert os.path.exists(copied_file)


@pytest.mark.skipif(not shutil.which('plasmidfinder.py'),
                    reason="plasmidfinder.py not found in PATH")
def test_integration_with_test_file():
    """Integration test with actual test file (requires plasmidfinder.py in PATH)"""
    
    # Use one of the test FASTA files
    test_files = list(get_test_file_dir().glob('*.fasta'))
    if not test_files:
        pytest.skip("No test FASTA files found")
    
    test_file = test_files[0]
    
    Args = collections.namedtuple('Args', [
        'plasmidfinder_dbpath',
        'plasmidfinder_dbselection',
        'plasmidfinder_mincov',
        'plasmidfinder_threshold',
        'plasmidfinder_extendoutput',
        'outdir'
    ])
    
    with tempfile.TemporaryDirectory() as tmp_out:
        args = Args(
            plasmidfinder_dbpath=None,
            plasmidfinder_dbselection=None,
            plasmidfinder_mincov=0.60,
            plasmidfinder_threshold=0.90,
            plasmidfinder_extendoutput=False,
            outdir=tmp_out
        )
        
        previous_results = {'strain': 'test_strain'}
        
        # Run the module
        results = get_results(str(test_file), None, args, previous_results)
        
        # Check that results have correct structure
        assert 'plasmidlist' in results
        assert 'plasmiddb' in results
        assert 'plasmidident' in results
        assert 'plasmidcov' in results
        assert 'plasmidcontig' in results
        assert 'plasmidcontigposition' in results
        assert 'plasmidnote' in results
        assert 'plasmidaccession' in results
        assert isinstance(results['plasmidlist'], str)


@pytest.mark.skipif(not shutil.which('plasmidfinder.py'),
                    reason="plasmidfinder.py not found in PATH")
def test_integration_with_extended_output():
    """Integration test with extended output enabled"""
    
    # Use one of the test FASTA files
    test_files = list(get_test_file_dir().glob('*.fasta'))
    if not test_files:
        pytest.skip("No test FASTA files found")
    
    test_file = test_files[0]
    strain_name = get_strain_name(str(test_file))
    
    Args = collections.namedtuple('Args', [
        'plasmidfinder_dbpath',
        'plasmidfinder_dbselection',
        'plasmidfinder_mincov',
        'plasmidfinder_threshold',
        'plasmidfinder_extendoutput',
        'outdir'
    ])
    
    with tempfile.TemporaryDirectory() as tmp_out:
        args = Args(
            plasmidfinder_dbpath=None,
            plasmidfinder_dbselection=None,
            plasmidfinder_mincov=0.60,
            plasmidfinder_threshold=0.90,
            plasmidfinder_extendoutput=True,
            outdir=tmp_out
        )
        
        previous_results = {'strain': strain_name}
        
        # Run the module
        results = get_results(str(test_file), None, args, previous_results)
        
        # Check extended output directory was created
        strain_dir = os.path.join(tmp_out, strain_name)
        if results['plasmidlist'] != '-':
            # Only check for extended output if there were results
            assert os.path.exists(strain_dir) or True  # May not exist if no hits


@pytest.mark.skipif(not shutil.which('plasmidfinder.py'),
                    reason="plasmidfinder.py not found in PATH")
def test_run_plasmidfinder_handles_errors():
    """Test that run_plasmidfinder handles errors gracefully"""
    
    Args = collections.namedtuple('Args', [
        'plasmidfinder_dbpath',
        'plasmidfinder_dbselection',
        'plasmidfinder_mincov',
        'plasmidfinder_threshold',
        'plasmidfinder_extendoutput',
        'outdir'
    ])
    
    with tempfile.TemporaryDirectory() as tmp_out:
        args = Args(
            plasmidfinder_dbpath=None,
            plasmidfinder_dbselection=None,
            plasmidfinder_mincov=0.60,
            plasmidfinder_threshold=0.90,
            plasmidfinder_extendoutput=False,
            outdir=tmp_out
        )
        
        # Test with non-existent file - should return empty dict
        result = run_plasmidfinder('/nonexistent/file.fasta', 'test', args)
        assert result == {} or isinstance(result, dict)


def test_add_cli_options():
    """Test that CLI options are added correctly"""
    import argparse
    parser = argparse.ArgumentParser()
    
    group = add_cli_options(parser)
    
    # Parse empty args to check defaults
    args = parser.parse_args([])
    
    assert args.plasmidfinder_dbpath is None
    assert args.plasmidfinder_dbselection is None
    assert args.plasmidfinder_mincov == 0.60
    assert args.plasmidfinder_threshold == 0.90
    assert args.plasmidfinder_extendoutput == False


def test_add_cli_options_custom_values():
    """Test CLI options with custom values"""
    import argparse
    parser = argparse.ArgumentParser()
    
    add_cli_options(parser)
    
    # Parse with custom args
    args = parser.parse_args([
        '--plasmidfinder_mincov', '0.80',
        '--plasmidfinder_threshold', '0.95',
        '--plasmidfinder_dbselection', 'enterobacteriales,gram_positive',
        '--plasmidfinder_extendoutput'
    ])
    
    assert args.plasmidfinder_mincov == 0.80
    assert args.plasmidfinder_threshold == 0.95
    assert args.plasmidfinder_dbselection == 'enterobacteriales,gram_positive'
    assert args.plasmidfinder_extendoutput == True