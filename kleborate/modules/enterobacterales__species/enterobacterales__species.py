"""
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
import os
import pathlib
import shutil
import sys


def description():
    return 'Mash-based species detection for enterobacterales species'


def prerequisite_modules():
    return []


def get_headers():
    full_headers = ['species', 'species_match']
    stdout_headers = ['species']
    return full_headers, stdout_headers


def test_add_cli_options():
    parser = argparse.ArgumentParser()
    add_cli_options(parser)
    assert '--enterobacterales__species_strong' in parser.format_help()
    assert '--enterobacterales__species_weak' in parser.format_help()


def add_cli_options(parser):
    module_name = os.path.basename(__file__)[:-3]
    group = parser.add_argument_group(f'{module_name} module')
    group.add_argument('--enterobacterales__species_strong', type=float, default=0.02,
                       help='Mash distance threshold for a strong species match')
    group.add_argument('--enterobacterales__species_weak', type=float, default=0.04,
                       help='Mash distance threshold for a weak species match')
    return group


def check_cli_options(args):
    if args.enterobacterales__species_strong <= 0.0 or args.enterobacterales__species_strong >= 1.0:
        sys.exit('Error: --enterobacterales__species_strong must be between 0.0 and 1.0')
    if args.enterobacterales__species_weak <= 0.0 or args.enterobacterales__species_weak >= 1.0:
        sys.exit('Error: --enterobacterales__species_weak must be between 0.0 and 1.0')
    if args.enterobacterales__species_weak < args.enterobacterales__species_strong:
        sys.exit('Error: --enterobacterales__species_weak must be greater than '
                 '--enterobacterales__species_strong')


def check_external_programs():
    if not shutil.which('mash'):
        sys.exit('Error: could not find mash')
    return ['mash']


def get_results(assembly, minimap2_index, args, previous_results):
    #print(previous_results)
    sketch_file = pathlib.Path(__file__).parents[0] / 'data' / 'species_mash_sketches.msh'
    species, distance = get_enterobacterales__species(assembly, sketch_file)
    if distance <= args.enterobacterales__species_strong:
        species_hit_strength = 'strong'
    elif distance <= args.enterobacterales__species_weak:
        species_hit_strength = 'weak'
    else:
        species = 'unknown'
        species_hit_strength = ''
    return {'species': species,
            'species_match': species_hit_strength}


def get_enterobacterales__species(assembly, sketch_file):
    f = os.popen('mash dist ' + str(sketch_file) + ' ' + str(assembly))
    best_species, best_distance = None, 1.0
    for line in f:
        line_parts = line.split('\t')
        reference = line_parts[0]
        if len(line_parts) >= 3:
            species = reference.split('/')[0]
            distance = float(line_parts[2])
            if distance < best_distance:
                best_distance = distance
                best_species = species
    f.close()
    best_species = clean_species_name(best_species)
    return best_species, best_distance


def clean_species_name(species):
    if species is None:
        return None
    species = species.replace('Escherichia_coli', 'Escherichia coli / Shigella')
    species = species.replace('Raoultella', 'Klebsiella (Raoultella)')
    species = species.replace('_', ' ')
    species = species.replace(' subsp ', ' subsp. ')
    if species.endswith(' unknown'):
        species = species.replace(' unknown', ' (unknown species)')
    return species
