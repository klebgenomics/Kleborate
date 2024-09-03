"""
Copyright 2024 Mary Maranga, Kat Holt, Ryan Wick
https://github.com/klebgenomics/KleborateModular/

This file is part of Kleborate. Kleborate is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Kleborate is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Kleborate. If
not, see <https://www.gnu.org/licenses/>.
"""

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


def is_ko_complex(results):
    """
    Returns True if the species call is in the Ko-complex, otherwise false.
    """
    assert 'species' in results
    species = results['species']
    if species.startswith('Klebsiella oxytoca'):
        return True
    if species.startswith('Klebsiella grimontii'):
        return True
    if species.startswith('Klebsiella michiganensis'):
        return True
    if species.startswith('Klebsiella pasteurii'):
        return True
    if species.startswith('Klebsiella huaxiensis'):
        return True
    if species.startswith('Klebsiella spallanzanii'):
        return True
    return False


def is_escherichia(results):
    """
    Returns True if the species call is in the Escherichia genus, otherwise false.
    """
    assert 'species' in results
    species = results['species']
    if species.startswith('Escherichia'):
        return True
    return False
