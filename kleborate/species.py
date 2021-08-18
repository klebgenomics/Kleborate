"""
Copyright 2020 Kat Holt
Copyright 2020 Ryan Wick (rrwick@gmail.com)
https://github.com/katholt/Kleborate/

This file is part of Kleborate. Kleborate is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Kleborate is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Kleborate. If
not, see <http://www.gnu.org/licenses/>.
"""

import os


def get_species_results(contigs, data_folder):
    species, species_hit_strength = get_klebsiella_species(contigs, data_folder)
    return {'species': species,
            'species_match': species_hit_strength}


def get_klebsiella_species(contigs, data_folder):
    f = os.popen('mash dist ' + data_folder + '/species_mash_sketches.msh ' + contigs)

    best_species = None
    best_distance = 1.0

    for line in f:
        line_parts = line.split('\t')
        reference = line_parts[0]
        if len(line_parts) < 4:
            continue
        species = reference.split('/')[0]
        distance = float(line_parts[2])

        # Fix up the species name formatting a bit.
        species = species.replace('Escherichia_coli', 'Escherichia coli / Shigella')
        species = species.replace('Raoultella', 'Klebsiella (Raoultella)')
        species = species.replace('_', ' ')
        species = species.replace(' subsp ', ' subsp. ')
        if species.endswith(' unknown'):
            species = species.replace(' unknown', ' (unknown species)')

        if distance < best_distance:
            best_distance = distance
            best_species = species

    f.close()

    if best_distance <= 0.02:
        return best_species, 'strong'
    elif best_distance <= 0.04:
        return best_species, 'weak'
    else:
        return 'unknown', ''


def is_kp_complex(results):
    """
    Returns True if the species call is in the Kp-complex, otherwise false.
    """
    assert 'species' in results
    species = results['species']
    if species.startswith('Klebsiella pneumoniae'):
        return True
    if species.startswith('Klebsiella quasipneumoniae'):
        return True
    if species.startswith('Klebsiella variicola'):
        return True
    if species.startswith('Klebsiella quasivariicola'):
        return True
    if species.startswith('Klebsiella africana'):
        return True
    return False
