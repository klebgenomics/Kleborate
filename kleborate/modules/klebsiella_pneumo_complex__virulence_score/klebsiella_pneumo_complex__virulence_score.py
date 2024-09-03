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


def description():
    return 'virulence score (0-5) for the Klebsiella pneumoniae species complex, based on the ' \
           'results of the abst, cbst and ybst modules'


def prerequisite_modules():
    return ['klebsiella__abst', 'klebsiella__cbst', 'klebsiella__ybst','klebsiella__rmst', 'klebsiella__smst']


def get_headers():
    full_headers = ['virulence_score', 'spurious_virulence_hits']
    stdout_headers = ['virulence_score']
    return full_headers, stdout_headers


def add_cli_options(parser):
    pass


def check_cli_options(args):
    pass


def check_external_programs():
    return []


def get_results(assembly, minimap2_index, args, previous_results):
    # spurious hits
    
    ybt = previous_results['klebsiella__ybst__spurious_ybt_hits']
    abst = previous_results['klebsiella__abst__spurious_abst_hits']
    clb = previous_results['klebsiella__cbst__spurious_clb_hits']
    rmst = previous_results['klebsiella__rmst__spurious_rmst_hits']
    smst = previous_results['klebsiella__smst__spurious_smst_hits']
    
    # Concatenate all lists
    all_hits = ybt + abst + clb + rmst + smst

    all_hits = [s for s in all_hits if s != '-']

    # Check if the resulting list is empty
    if not all_hits:
        all_hits = '-'
    else:
        all_hits = ''.join(all_hits)

    # virulence score
    has_ybt = (previous_results['klebsiella__ybst__YbST'] != 0)
    has_aero = (previous_results['klebsiella__abst__AbST'] != 0)
    has_coli = (previous_results['klebsiella__cbst__CbST'] != 0)


    # Calculate virulence score
    if has_coli and has_aero:
        virulence_score = '5'
    elif has_aero and has_ybt:
        virulence_score = '4'
    elif has_aero:
        virulence_score = '3'
    elif has_coli:
        virulence_score = '2'
    elif has_ybt:
        virulence_score = '1'
    else:
        virulence_score = '0'

    # Return both spurious hits and virulence score
    return {
        'spurious_virulence_hits': all_hits,
        'virulence_score': virulence_score
    }


# def get_results(assembly, minimap2_index, args, previous_results):
#     has_ybt = (previous_results['klebsiella__ybst__YbST'] != 'NA')
#     has_aero = (previous_results['klebsiella__abst__AbST'] != 'NA')
#     has_coli = (previous_results['klebsiella__cbst__CbST'] != 'NA')

#     if has_coli and has_aero:
#         return {'virulence_score': '5'}
#     elif has_aero and has_ybt:
#         return {'virulence_score': '4'}
#     elif has_aero:
#         return {'virulence_score': '3'}
#     elif has_coli:
#         return {'virulence_score': '2'}
#     elif has_ybt:
#         return {'virulence_score': '1'}
#     else:
#         return {'virulence_score': '0'}
