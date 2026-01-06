"""
Module for detecting plasmid replicons using PlasmidFinder database

This module integrates PlasmidFinder into Kleborate for detecting plasmid replicons.
It is an optional module that must be explicitly called by the user.

Key features:
- Uses PlasmidFinder from conda environment
- Supports multiple input files via Kleborate's -a option
- Custom parameters for PlasmidFinder configuration
- Handles PlasmidFinder's default output behavior safely using temp directories

Copyright 2025 Minh-Quan Ton-Ngoc
https://github.com/tnmquann/Kleborate/

This file is part of Kleborate. Kleborate is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Kleborate is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Kleborate. If
not, see <https://www.gnu.org/licenses/>.
"""

import os
import pathlib
import shutil
import sys
import subprocess
import json
import tempfile
import gzip


def description():
    """
    Returns a string describing what this module does.
    """
    return 'detects plasmid replicons using PlasmidFinder database'


def prerequisite_modules():
    """
    Returns a list of modules that must be run before this module.
    This module is optional and has no prerequisites.
    """
    return []


def get_headers():
    """
    Returns two lists of column headers:
    1. full_headers: all columns this module will output
    2. stdout_headers: subset of columns to display in terminal
    
    Output columns description:
    - plasmidlist: List of plasmid replicon names (semicolon-separated)
    - plasmiddb: Database category for each plasmid (semicolon-separated)
    - plasmidident: Identity percentage for each plasmid (semicolon-separated)
    - plasmidcov: Coverage percentage for each plasmid (semicolon-separated)
    - plasmidcontig: Contig name for each plasmid hit (semicolon-separated)
    - plasmidcontigposition: Position in contig for each hit (semicolon-separated)
    - plasmidnote: Additional notes for each plasmid (semicolon-separated)
    - plasmidaccession: Accession number for each plasmid (semicolon-separated)
    """
    full_headers = [
        'plasmidlist',
        'plasmiddb',
        'plasmidident',
        'plasmidcov',
        'plasmidcontig',
        'plasmidcontigposition',
        'plasmidnote',
        'plasmidaccession'
    ]
    stdout_headers = ['plasmidlist']
    return full_headers, stdout_headers


def add_cli_options(parser):
    """
    Adds module-specific command-line options.
    
    These options allow customization of PlasmidFinder behavior:
    - --plasmidfinder_dbpath: Custom database path
    - --plasmidfinder_dbselection: Select specific databases
    - --plasmidfinder_mincov: Minimum coverage threshold
    - --plasmidfinder_threshold: Minimum identity threshold
    - --plasmidfinder_extendoutput: Generate extended output files
    """
    module_name = os.path.basename(__file__)[:-3]
    group = parser.add_argument_group(f'{module_name} module')
    
    group.add_argument('--plasmidfinder_dbpath', type=str, default=None,
                       help='Path to custom PlasmidFinder database directory')
    group.add_argument('--plasmidfinder_dbselection', type=str, default=None,
                       help='Comma-separated list of databases to use (default: all)')
    group.add_argument('--plasmidfinder_mincov', type=float, default=0.60,
                       help='Minimum coverage for plasmid replicon match (default: 0.60)')
    group.add_argument('--plasmidfinder_threshold', type=float, default=0.90,
                       help='Minimum identity threshold for plasmid replicon match (default: 0.90)')
    group.add_argument('--plasmidfinder_extendoutput', action='store_true',
                       help='Generate extended output with alignment files and TSV results '
                            '(output to {outdir}/{strain}/)')
    
    return group


def check_cli_options(args):
    """
    Validates module-specific command-line options.
    
    Raises SystemExit if:
    - mincov is not between 0 and 1
    - threshold is not between 0 and 1
    - dbpath is specified but doesn't exist
    """
    if args.plasmidfinder_mincov < 0 or args.plasmidfinder_mincov > 1:
        sys.exit('Error: --plasmidfinder_mincov must be between 0 and 1')
    if args.plasmidfinder_threshold < 0 or args.plasmidfinder_threshold > 1:
        sys.exit('Error: --plasmidfinder_threshold must be between 0 and 1')
    if args.plasmidfinder_dbpath and not os.path.exists(args.plasmidfinder_dbpath):
        sys.exit(f'Error: PlasmidFinder database path does not exist: {args.plasmidfinder_dbpath}')


def check_external_programs():
    """
    Returns a list of external programs required by this module.
    
    Checks that plasmidfinder.py is available in PATH.
    This is typically installed via conda: `conda install -c bioconda plasmidfinder`
    """
    if not shutil.which('plasmidfinder.py'):
        sys.exit('Error: could not find plasmidfinder.py. Please install via conda: '
                 'conda install -c bioconda plasmidfinder')
    return ['plasmidfinder.py']


def data_dir():
    """
    Returns the path to this module's data directory.
    """
    return pathlib.Path(__file__).parents[0] / 'data'


def get_strain_name(assembly_path):
    """
    Extracts strain name from assembly filename.
    Handles both regular and gzipped files.
    
    Examples:
    - /path/to/sample.fasta -> sample
    - /path/to/sample.fasta.gz -> sample
    """
    basename = os.path.basename(assembly_path)
    
    # Handle gzipped files
    if basename.endswith('.gz'):
        basename = basename[:-3]
    
    # Remove common extensions
    for ext in ['.fasta', '.fa', '.fna', '.fas']:
        if basename.endswith(ext):
            basename = basename[:-len(ext)]
            break
    
    return basename


def get_results(assembly, minimap2_index, args, previous_results):
    """
    Main function that performs PlasmidFinder analysis.
    
    This function is called by Kleborate's main loop for each assembly.
    It handles multiple input files automatically as Kleborate iterates
    through assemblies provided with the -a option.
    
    Parameters:
    - assembly: path to assembly FASTA file (may be gzipped)
    - minimap2_index: pre-built minimap2 index (not used by this module)
    - args: command-line arguments namespace
    - previous_results: dict of results from prerequisite modules
    
    Returns:
    - dict mapping column headers to results
    """
    # Get strain name from previous results or derive from assembly filename
    strain = previous_results.get('strain', get_strain_name(assembly))
    
    # Run PlasmidFinder and get JSON results
    plasmidfinder_results = run_plasmidfinder(assembly, strain, args)
    
    # Parse and format results into standard columns
    results = parse_plasmidfinder_output(plasmidfinder_results)
    
    return results


def run_plasmidfinder(assembly, strain, args):
    """
    Runs PlasmidFinder on the assembly and returns the parsed JSON output.
    
    Key design decisions to handle PlasmidFinder's default behavior:
    1. Always use a temporary directory for PlasmidFinder output to avoid
       writing 'data.json' in the current working directory.
    2. Use the temp directory for both -o and --tmp_dir options.
    3. Only copy extended output files to user-specified location if requested.
    
    Parameters:
    - assembly: path to assembly FASTA file
    - strain: strain name for output directory naming
    - args: command-line arguments namespace
    
    Returns:
    - dict: parsed JSON data from PlasmidFinder
    """
    # Create temporary directory for PlasmidFinder output
    # This prevents PlasmidFinder from creating 'data.json' in CWD
    with tempfile.TemporaryDirectory() as tmp_dir:
        # Build PlasmidFinder command
        cmd = [
            'plasmidfinder.py',
            '-i', assembly,
            '-o', tmp_dir,  # Output to temp directory
            '-l', str(args.plasmidfinder_mincov),
            '-t', str(args.plasmidfinder_threshold),
            '--tmp_dir', tmp_dir  # Use same temp dir for intermediate files
        ]
        
        # Add optional database path
        if args.plasmidfinder_dbpath:
            cmd.extend(['-p', args.plasmidfinder_dbpath])
        
        # Add database selection
        if args.plasmidfinder_dbselection:
            cmd.extend(['-d', args.plasmidfinder_dbselection])
        
        # Add extended output flag if requested
        if args.plasmidfinder_extendoutput:
            cmd.append('-x')
        
        # Run PlasmidFinder
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
        except subprocess.CalledProcessError as e:
            # Return empty results on error instead of exiting
            # This allows Kleborate to continue processing other assemblies
            print(f'Warning: PlasmidFinder failed for {strain}: {e.stderr}', file=sys.stderr)
            return {}
        except FileNotFoundError:
            print(f'Warning: plasmidfinder.py not found in PATH', file=sys.stderr)
            return {}
        
        # Read JSON output
        json_file = os.path.join(tmp_dir, 'data.json')
        if not os.path.exists(json_file):
            print(f'Warning: PlasmidFinder did not produce output for {strain}', file=sys.stderr)
            return {}
        
        with open(json_file, 'r') as f:
            json_data = json.load(f)
        
        # Copy extended output files to user-specified output directory if requested
        if args.plasmidfinder_extendoutput and hasattr(args, 'outdir'):
            _copy_extended_output(tmp_dir, args.outdir, strain)
        
        return json_data


def _copy_extended_output(src_dir, outdir, strain):
    """
    Copies extended output files from PlasmidFinder to the user's output directory.
    
    Creates a subdirectory named after the strain:
    {outdir}/{strain}/plasmidfinder_*
    
    Files copied:
    - results_tab.tsv -> plasmidfinder_results_tab.tsv
    - Hit_in_genome_seq.fsa -> plasmidfinder_Hit_in_genome_seq.fsa
    - Plasmid_seqs.fsa -> plasmidfinder_Plasmid_seqs.fsa
    - results.txt -> plasmidfinder_results.txt
    - data.json -> plasmidfinder_data.json
    
    Parameters:
    - src_dir: temporary directory containing PlasmidFinder output
    - outdir: user-specified output directory
    - strain: strain name for subdirectory naming
    """
    strain_output_dir = os.path.join(outdir, strain)
    os.makedirs(strain_output_dir, exist_ok=True)
    
    # List of files to copy (source_name, dest_prefix)
    extended_files = [
        'results_tab.tsv',
        'Hit_in_genome_seq.fsa',
        'Plasmid_seqs.fsa',
        'results.txt',
        'data.json'
    ]
    
    for filename in extended_files:
        src = os.path.join(src_dir, filename)
        if os.path.exists(src):
            dst = os.path.join(strain_output_dir, f'plasmidfinder_{filename}')
            shutil.copy2(src, dst)


def parse_plasmidfinder_output(json_data):
    """
    Parses PlasmidFinder JSON output and extracts relevant information.
    
    The JSON structure from PlasmidFinder is:
    {
        "plasmidfinder": {
            "results": {
                "Database_Category": {
                    "database_name": {
                        "hit_id": {
                            "plasmid": "name",
                            "identity": 99.15,
                            "coverage": 82.98,
                            "contig_name": "contig1",
                            "positions_in_contig": "5041466..5041938",
                            "note": "",
                            "accession": "JN420336"
                        }
                    }
                }
            }
        }
    }
    
    Parameters:
    - json_data: parsed JSON dict from PlasmidFinder
    
    Returns:
    - dict with standardized column names, values semicolon-separated for multiple hits
    """
    # Initialize with default values (dash indicates no results)
    results = {
        'plasmidlist': '-',
        'plasmiddb': '-',
        'plasmidident': '-',
        'plasmidcov': '-',
        'plasmidcontig': '-',
        'plasmidcontigposition': '-',
        'plasmidnote': '-',
        'plasmidaccession': '-'
    }
    
    # Return defaults if no plasmidfinder results
    if not json_data or 'plasmidfinder' not in json_data:
        return results
    
    pf_results = json_data['plasmidfinder'].get('results', {})
    
    # Collect all hits across all databases
    plasmids = []
    databases = []
    identities = []
    coverages = []
    contigs = []
    positions = []
    notes = []
    accessions = []
    
    # Iterate through database categories (e.g., "Enterobacteriales", "Gram Positive")
    for db_category, db_data in pf_results.items():
        # db_data can be a dict of databases or nested structure
        if not isinstance(db_data, dict):
            continue
            
        # Iterate through databases or hits within the category
        for db_name, hits in db_data.items():
            # Skip if no hits found (string value like "No hit found")
            if isinstance(hits, str):
                continue
            
            # Skip if not a dict (unexpected format)
            if not isinstance(hits, dict):
                continue
            
            # Check if this level contains hit data directly or is another nested level
            # Hit data has keys like 'plasmid', 'identity', etc.
            if 'plasmid' in hits:
                # This is a single hit at this level
                hit_info = hits
                plasmids.append(str(hit_info.get('plasmid', '-')) or '-')
                databases.append(db_category or '-')
                identities.append(str(hit_info.get('identity', '-')) or '-')
                coverages.append(str(hit_info.get('coverage', '-')) or '-')
                contigs.append(str(hit_info.get('contig_name', '-')) or '-')
                positions.append(str(hit_info.get('positions_in_contig', '-')) or '-')
                notes.append(str(hit_info.get('note', '-')) or '-')
                accessions.append(str(hit_info.get('accession', '-')) or '-')
            else:
                # This is a collection of hits
                for hit_id, hit_info in hits.items():
                    if not isinstance(hit_info, dict):
                        continue
                    
                    plasmids.append(str(hit_info.get('plasmid', '-')) or '-')
                    databases.append(db_category or '-')
                    identities.append(str(hit_info.get('identity', '-')) or '-')
                    coverages.append(str(hit_info.get('coverage', '-')) or '-')
                    contigs.append(str(hit_info.get('contig_name', '-')) or '-')
                    positions.append(str(hit_info.get('positions_in_contig', '-')) or '-')
                    notes.append(str(hit_info.get('note', '-')) or '-')
                    accessions.append(str(hit_info.get('accession', '-')) or '-')
    
    # Format results with semicolon separator if there are hits
    if plasmids:
        results['plasmidlist'] = ';'.join(plasmids)
        results['plasmiddb'] = ';'.join(databases)
        results['plasmidident'] = ';'.join(identities)
        results['plasmidcov'] = ';'.join(coverages)
        results['plasmidcontig'] = ';'.join(contigs)
        results['plasmidcontigposition'] = ';'.join(positions)
        results['plasmidnote'] = ';'.join(notes)
        results['plasmidaccession'] = ';'.join(accessions)
    
    return results