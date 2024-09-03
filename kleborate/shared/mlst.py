"""
This file contains code for standard assigning a standard MLST scheme (e.g. 7-gene Klebsiella
pneumoniae ST) to an assembly.

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

import re

from .alignment import align_query_to_ref, truncation_check


def mlst(assembly_path, minimap2_index, profiles_path, allele_paths, gene_names, extra_info,
         min_identity, min_coverage, required_exact_matches, check_for_truncation=False):
    """
    This function takes:
    * assembly_path: a path for an assembly in FASTA format
    * minimap2_index: a path for the assembly's minimap2 index (for faster alignment)
    * profiles_path: a path for the MLST profiles file in TSV format
    * allele_paths: a dictionary {gene name: path for the allele FASTA file}
    * gene_names: a list of the gene names in the MLST scheme
    * extra_info: the name of an additional extra-info column in the MLST scheme (or None)
    * min_identity: hits with a lower percent identity than this are discarded
    * min_coverage: hits with a lower percent coverage than this are discarded
    * required_exact_matches: at least this many alleles must be an exact match to assign an ST
    * check_for_truncation: if true, truncation strings will be added to the allele numbers

    This function returns:
    * the best matching ST profile (e.g. 'ST123', 'ST456-1LV' or 'NA')
    * the extra-info value for the best ST (if used, otherwise an empty string)
    * a dictionary of allele numbers in str format {gene name: allele number}
    """
    profiles = load_st_profiles(profiles_path, gene_names, extra_info)
    hits_per_gene = {g: align_query_to_ref(allele_paths[g], assembly_path,
                                           ref_index=minimap2_index, min_identity=min_identity,
                                           min_query_coverage=min_coverage) for g in gene_names}
    return run_single_mlst(profiles, hits_per_gene, gene_names, required_exact_matches,
                           check_for_truncation)


def run_single_mlst(profiles, hits_per_gene, gene_names, required_exact_matches,
                    check_for_truncation=False, report_incomplete=False):
    """
    This function is factored out because it is also called by the multi_mlst.py file.
    """
    best_hits_per_gene = {gene: get_best_hits(hits_per_gene[gene]) for gene in gene_names}
    st, alleles, extra_info = get_best_matching_profile(profiles, gene_names, best_hits_per_gene)
    best_hit_per_gene = get_best_hit_per_gene(gene_names, best_hits_per_gene, alleles)

    exact_matches, lv_count, allele_numbers = 0, 0, {}
    any_truncations, any_missing = False, False
    for gene_name, st_allele in zip(gene_names, alleles):
        hit = best_hit_per_gene[gene_name]
        hit_allele = number_from_hit(hit)

        if hit is None:
            allele_numbers[gene_name] = '-'
            any_missing = True
        elif hit.is_exact():
            allele_numbers[gene_name] = str(hit_allele)
        else:
            allele_numbers[gene_name] = str(hit_allele) + '*'

        if check_for_truncation and hit is not None:
            truncation_suffix = truncation_check(hit)[0]
            allele_numbers[gene_name] += truncation_suffix
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


def load_st_profiles(database_path, gene_names, extra_info_name):
    """
    This function reads through a tab-delimited MLST database file where the first column is the ST
    number and the subsequent columns are allele numbers, with a final optional extra-info column.
    The first line of the file should be a header: 'ST', followed by gene names, optionally
    followed by the extra-info column name . All other lines should only contain positive integers
    for the ST and gene columns, and the final optional column can contain anything.

    This function returns a list of ST profiles, where each value is a tuple:
    (ST number, list of allele numbers, extra info)
    """
    profiles, first_line = [], True
    with open(database_path, 'r') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if first_line:
                assert parts[0] == 'ST'
                if extra_info_name is None:
                    assert parts[1:] == gene_names
                else:
                    assert parts[1:] == gene_names + [extra_info_name]
                first_line = False
            else:
                st = int(parts[0])
                if extra_info_name is None:
                    alleles = [int(a) for a in parts[1:]]
                    extra_info = None
                else:
                    alleles = [int(a) for a in parts[1:-1]]
                    extra_info = parts[-1]
                profiles.append((st, alleles, extra_info))
    return profiles


def get_best_hits(hits):
    """
    Given a bunch of hits to an allele, this function returns a list of the best hits. 'Best' is
    defined as highest identity. If there is a tie, hits with higher alignment scores are
    preferred. Usually this results in just a single hit, but if there is still a tie (i.e. same
    identity and same alignment score), then multiple hits can be returned.
    """
    if not hits:
        return []
    best_identity = max(h.percent_identity for h in hits)
    hits = [h for h in hits if h.percent_identity == best_identity]
    best_score = max(h.alignment_score for h in hits)
    return [h for h in hits if h.alignment_score == best_score]


def number_from_hit(hit):
    """
    Given an alignment, returns the numerical part of the query name as an int. E.g if the query
    name was gapA_234, it will return 234.
    """
    if hit is None:
        return 0
    try:
        last_part = hit.query_name.split('_')[-1]
        return int(last_part)
    except ValueError:
        pass
    try:
        return int(re.sub('[^0-9]', '', hit.query_name))
    except ValueError:
        return 0


def get_best_matching_profile(profiles, gene_names, best_hits_per_gene):
    """
    This function looks for an ST which best matches the hits. Each ST is scored based on the
    number of genes which have a matching hit (i.e. a hit with the same number). So the score for
    any ST can be 0 to the number of genes in the scheme. STs earlier in the profiles are
    preferred, so if an assembly matches multiple STs equally well, this function will return
    whichever is first in the profiles.
    """
    best_st, best_extra_info, best_matches = 0, '', 0
    best_alleles = [0] * len(gene_names)
    for st, alleles, extra_info in profiles:
        matches = 0
        for gene_name, allele in zip(gene_names, alleles):
            if any(allele == number_from_hit(h) for h in best_hits_per_gene[gene_name]):
                matches += 1
        if matches > best_matches:
            best_st, best_alleles, best_extra_info, best_matches = st, alleles, extra_info, matches
    return best_st, best_alleles, best_extra_info


def get_best_hit_per_gene(gene_names, best_hits_per_gene, alleles):
    """
    This function takes the best_hits_per_gene dict, and for any genes that have multiple hits, it
    chooses a single hit. Which hit is chosen is first based on the ST call, i.e. if one of the
    hits matches the ST allele, that's chosen. If none of the hits match the ST allele, then the
    lowest allele number is chosen.
    """
    best_hit_per_gene = {}
    for gene_name, allele in zip(gene_names, alleles):
        best_hits = best_hits_per_gene[gene_name]
        hits_matching_st = [h for h in best_hits if number_from_hit(h) == allele]
        if not best_hits:
            best_hit_per_gene[gene_name] = None
        elif hits_matching_st:
            best_hit_per_gene[gene_name] = hits_matching_st[0]
        else:
            best_hit_per_gene[gene_name] = sorted(best_hits, key=lambda h: number_from_hit(h))[0]
    return best_hit_per_gene
