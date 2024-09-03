"""
Copyright 2023 Kat Holt, Ryan Wick (rrwick@gmail.com), Mary Maranga (gathonimaranga@gmail.com)
https://github.com/klebgenomics/KleborateModular/

This file is part of Kleborate. Kleborate is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Kleborate is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Kleborate. If
not, see <https://www.gnu.org/licenses/>.
"""
import re

from ...shared.alignment import align_query_to_ref, truncation_check
from ...shared.mlst import get_best_hits, number_from_hit


def load_st_profiles(database_path, gene_name, extra_info_name=None):
    """
    This function returns a list of WZI ST profiles, where each value is a tuple:
    (ST number, list of allele numbers, extra info)
    """
    profiles, first_line = [], True
    with open(database_path, 'r') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if first_line:
                assert parts[0] == 'ST'
                if extra_info_name is None:
                    assert parts[1] == gene_name
                else:
                    assert parts[1] == gene_name and parts[2] == extra_info_name
                first_line = False
            else:
                st = int(parts[0])
                allele = [int(parts[1])]  # Store allele number as a list
                extra_info = parts[2] if extra_info_name else None
                profiles.append((st, allele, extra_info))
    return profiles


def get_best_matching_profile(profiles, gene_name, best_hits_per_gene):
    """
    This function looks for an ST which best matches the hits. Each ST is scored based on the
    number of genes which have a matching hit (i.e. a hit with the same number). So the score for
    any ST can be 0 to the number of genes in the scheme. STs earlier in the profiles are
    preferred, so if an assembly matches multiple STs equally well, this function will return
    whichever is first in the profiles.
    """
    best_st, best_extra_info, best_matches = 0, '', 0
    best_allele = 0
    for st, alleles, extra_info in profiles:
        matches = 0
        for allele in alleles:
            if any(allele == number_from_hit(h) for h in best_hits_per_gene):
                matches += 1
                if matches > best_matches:
                    best_st, best_allele, best_extra_info, best_matches = st, alleles, extra_info, matches
    return best_st, best_allele, best_extra_info


def get_best_hit_per_gene(gene_name, best_hits_per_gene, alleles):
    """
    This function processes a list of best hits for a single gene. It chooses a single hit based
    on the ST allele provided. If one of the hits matches the ST allele, that's chosen. If none
    of the hits match the ST allele, then the hit with the lowest allele number is chosen.
    """
    if isinstance(alleles, (list, tuple)):
        allele = alleles[0]
    else:
        allele = alleles

    hits_matching_st = [h for h in best_hits_per_gene if number_from_hit(h) == allele]

    if not best_hits_per_gene:
        return None
    elif hits_matching_st:
        return hits_matching_st[0]
    else:
        return sorted(best_hits_per_gene, key=lambda h: number_from_hit(h))[0]


def run_single_mlst(profiles, hits_per_gene, gene_name, required_exact_matches,
                    check_for_truncation=False, report_incomplete=False):

    best_hits_per_gene = get_best_hits(hits_per_gene) 
    st, alleles, extra_info = get_best_matching_profile(profiles, gene_name, best_hits_per_gene)
    best_hit_per_gene = get_best_hit_per_gene(gene_name, best_hits_per_gene, alleles)
    
    exact_matches, lv_count, allele_numbers = 0, 0, {}
    any_truncations, any_missing = False, False
    hit = best_hit_per_gene
    if isinstance(alleles, (list, tuple)):
        st_allele = alleles[0]
    else:
        st_allele = alleles
    hit_allele = number_from_hit(hit)
    if hit is None:
        allele_numbers = '-'
        any_missing = True
    elif hit.is_exact():
        allele_numbers = str(hit_allele)
    else:
        allele_numbers = str(hit_allele) + '*'

        
    if check_for_truncation and hit is not None:
        truncation_suffix = truncation_check(hit)[0]
        allele_numbers += truncation_suffix
        if truncation_suffix:
            any_truncations = True

    if hit is not None and hit.is_exact() and st_allele == hit_allele:
        exact_matches += 1
    else:
        lv_count += 1
            
    if exact_matches < required_exact_matches:
        st, extra_info = 'NA', '-'
    
        
    elif lv_count == 0:
        st = 'ST' + str(st)
        
    else:
        st = 'ST' + str(st) + f'-{lv_count}LV'
        
        
    if extra_info != '-' and any_missing and report_incomplete:
        extra_info += ' (incomplete)'
        if extra_info != '-' and any_truncations:
            extra_info += ' (truncated)'

    return st, extra_info, allele_numbers


def mlst(assembly_path, minimap2_index, profiles_path, allele_path, gene_name, extra_info,
         min_identity, min_coverage, required_exact_matches, check_for_truncation=False):
    """
    This function is different from the mlst function in shared folder as it takes a single gene
    This function takes:
    * assembly_path: a path for an assembly in FASTA format
    * minimap2_index: a path for the assembly's minimap2 index (for faster alignment)
    * profiles_path: a path for the MLST profiles file in TSV format
    * allele_path:takes single allele
    * gene_name:takes single gene
    * extra_info: the name of an additional extra-info column in the ST scheme
    * min_identity: hits with a lower percent identity than this are discarded
    * min_coverage: hits with a lower percent coverage than this are discarded
    * required_exact_matches: at least this many alleles must be an exact match to assign an ST
    * check_for_truncation: if true, truncation strings will be added to the allele numbers

    This function returns:
    * the best matching profile
    * allele numbers
    """
    profiles = load_st_profiles(profiles_path, gene_name, extra_info)
    hits_per_gene = align_query_to_ref(allele_path, assembly_path,
                                           ref_index=minimap2_index, min_identity=min_identity,
                                           min_query_coverage=min_coverage)
    return run_single_mlst(profiles, hits_per_gene, gene_name, required_exact_matches,
                           check_for_truncation)
