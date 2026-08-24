"""
This module contains classes for interacting with bacterial genome assemblies and contigs and a pipeline
to type them.

Copyright 2026 Mary Maranga, Tom Stanton
https://github.com/klebgenomics/Kleborate/
https://github.com/klebgenomics/Kaptive

This file is part of Kaptive. Kaptive is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Kaptive is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Kaptive.
If not, see <https://www.gnu.org/licenses/>.
"""

import os
os.environ['KMP_WARNINGS'] = '0'
import sys
from pathlib import Path

from kaptive.db import Database
from kaptive.db.manager import DatabaseManager
from kaptive.core.genome import GenomeAssembly
from kaptive.serotyping import Serotyper


def description():
    return 'In silico serotyping of K and L locus for the Klebsiella pneumoniae species complex'


def prerequisite_modules():
    return []


def get_headers():
    full_headers = [
        'K_locus', 'K_type', 'K_locus_confidence', 'K_locus_problems', 'K_locus_identity',
        'K_Missing_expected_genes', 'K_Database_name', 'K_Database_version',
        'O_locus', 'O_type', 'O_locus_confidence', 'O_locus_problems', 
        'O_locus_identity', 'O_Missing_expected_genes', 'O_Database_name', 'O_Database_version',
        'Kaptive version'
    ]
    stdout_headers = []
    return full_headers, stdout_headers


def add_cli_options(parser):
    module_name = os.path.basename(__file__)[:-3]
    group = parser.add_argument_group(f'{module_name} module')
    group.add_argument('--kpsc-k-db', type=str, default='kpsc_k', metavar='',
                       help="Kaptive database for K-locus typing (default: kpsc_k)")
    group.add_argument('--kpsc-o-db', type=str, default='kpsc_o', metavar='',
                       help="Kaptive database for O-locus typing (default: kpsc_o)")
    return group


def load_or_install_db(db_input):
    """Loads/downloads a Kaptive database"""
    if isinstance(db_input, Database):
        return db_input
    
    db_str = str(db_input)
    
    if os.path.exists(db_str):
        return Database.load(db_str)

    db_mgr = DatabaseManager()
    
    db_obj = db_mgr.get(db_str)
    
    if isinstance(db_obj, Database):
        return db_obj
    return Database.load(db_obj)



def check_cli_options(args):
    # Reads from module-specific args namespace
    args.k_db = load_or_install_db(args.kpsc_k_db)
    args.o_db = load_or_install_db(args.kpsc_o_db)

    args.k_typer = Serotyper(args.k_db)
    args.o_typer = Serotyper(args.o_db)


def check_external_programs():
    return []


def extract_fields(prefix, result, full_headers):
    """Pull the fields from SerotypingResult into results_dict entries."""
    fields = {}
    if result is None:
        return fields

    fields[f'{prefix}_locus'] = result.best_locus_name
    fields[f'{prefix}_type'] = result.phenotype
    fields[f'{prefix}_locus_confidence'] = 'Typeable' if result.typeable else 'Untypeable'
    fields[f'{prefix}_locus_problems'] = result.problems.to_symbols().decode('utf-8')
    fields[f'{prefix}_locus_identity'] = '%.2f%%' % result.percent_identity
    fields[f'{prefix}_Missing_expected_genes'] = ';'.join(result.missing_expected_genes)
    fields[f'{prefix}_Database_name'] = result.database_name
    fields[f'{prefix}_Database_version'] = result.database_version
    fields['Kaptive version'] = result.kaptive_version

    for h in fields.keys():
        if h not in full_headers:
            sys.exit(f'Error: results contained a value ({h}) that is not covered by the full headers')

    return fields


def get_results(assembly, minimap2_index, args, previous_results):
    full_headers, _ = get_headers()
    assembly_path = Path(assembly)

    genome = GenomeAssembly.ensure(assembly_path)

    results_dict = {}

    k_result = args.k_typer(genome)
    if k_result is not None:
        results_dict.update(extract_fields('K', k_result, full_headers))
    else:
        print("Warning: No gene alignments sufficient for K-locus typing. Skipping k_results processing.")

    o_result = args.o_typer(genome)
    if o_result is not None:
        results_dict.update(extract_fields('O', o_result, full_headers))
    else:
        print("Warning: No gene alignments sufficient for O-locus typing. Skipping o_results processing.")

    results_dict = {k: (v if v else '-') for k, v in results_dict.items()}
    return results_dict

