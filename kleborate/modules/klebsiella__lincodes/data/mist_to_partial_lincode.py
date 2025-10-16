#! /usr/bin/env python3
"""
Converts the MiST output to a partial LINcode.

The JSON input file should be generated with the MiST software and the Institut Pasteur Klebsiella scgMLST629_S scheme
(with cgMLST profiles!).

Example command:
mist call --fasta my_kleb_genome.fasta --db kleb_scgmlst_s --output my_kleb_genome.json
"""
import json
import sys

MISMATCH_CUTOFFS = [629, 610, 585, 190, 43, 10, 7, 4, 2, 1]
METADATA = [
    {'key': 'Phylogroup', 'bin_min': 2},
    {'key': 'Sublineage', 'bin_min': 3},
    {'key': 'Clonal group', 'bin_min': 4}
]


def determine_bin(nb_diff) -> int:
    """
    Determine the deepest bin index in the LINcode that can be reliably assigned.
    :param nb_diff: Number of differences
    :return: Bin number
    """
    for i, threshold in enumerate(MISMATCH_CUTOFFS):
        if i == len(MISMATCH_CUTOFFS) - 1:
            return i
        if threshold >= nb_diff >= MISMATCH_CUTOFFS[i+1]:
            return i
    raise ValueError('Cannot determine bin')

def create_lin_code(lin_code_in: str, bin_: int) -> str:
    """
    Construct the partial LINcode, inserting wildcard symbols for bins that cannot be reliably determined.
    :param lin_code_in: Input LINcode
    :param bin_: Bin number
    :return: Partial LINcode
    """
    lin_code_parts = lin_code_in.split('_')
    return '_'.join(lin_code_parts[:bin_] + (['*'] * (9 - bin_)))


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise ValueError('Usage: mist_to_partial_lincode.py <input_file.json>')
    with open(sys.argv[1]) as handle:
        data_in = json.load(handle)

    # Check if the input contains profile information
    if data_in['profile'] is None:
        raise ValueError('Profile information is missing')

    # Log the results
    best_matching_st = data_in['profile']['name']
    print(f'Best matching: scgST-{best_matching_st}')

    nb_loci = len(data_in['profile']['alleles'])
    try:
        nb_matches = data_in['profile']['nb_matches']
    except KeyError:
        nb_matches = int(data_in['profile']['pct_match'] * nb_loci / 100)
    print(f'Number of matches: {nb_matches}/{nb_loci}')

    # Parse the ST metadata
    metadata_st = {k: v for k, v in data_in['profile']['metadata']}
    lin_code = metadata_st['LINcode']
    print(f'LINcode for scgMST-{best_matching_st}: {lin_code}')

    # Construct the partial LINcode
    bin_nb = determine_bin(nb_loci - nb_matches)
    partial_lin_code = create_lin_code(lin_code, bin_nb)
    print(f'Partial LINcode for input strain: {partial_lin_code}')

    # Report the associated metadata
    for metadata_col in METADATA:
        value = data_in['profile']
        if metadata_col['bin_min'] > bin_nb:
            value = 'Cannot be determined'
        else:
            value = metadata_st[metadata_col['key']]
        print(f"{metadata_col['key']}: {value}")
