"""
This file carries out multi-MLST functions: similar to the regular MLST (found in mlst.py) but
allowing for multiple STs per genome. This is useful for some virulence loci which can appear more
than once per genome (e.g. on the chromosome and on a plasmid).

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

from .alignment import align_query_to_ref
from .mlst import load_st_profiles, run_single_mlst
from .alignment import truncation_check

def multi_mlst(assembly_path, minimap2_index, profiles_path, allele_paths, gene_names, extra_info,
               min_identity, min_coverage, required_exact_matches, check_for_truncation=False,
               report_incomplete=False, min_spurious_identity=None, min_spurious_coverage=None,):
    """
    This function takes and returns the same things as the mlst function in mlst.py. However, it
    will look for cases where multiple contigs have hits for the full set of MLST genes, and in
    that case, MLST is run on each of them. Otherwise, it behaves like normal MLST.
    """
    profiles = load_st_profiles(profiles_path, gene_names, extra_info)
    
    if min_spurious_coverage is not None:
        hits_per_gene = {g: align_query_to_ref(allele_paths[g], assembly_path,
                                               ref_index=minimap2_index, min_identity=min_spurious_identity,
                                               min_query_coverage=min_spurious_coverage) for g in gene_names}

        spurious_hits = {g: [h for h in hits_per_gene[g] 
                     if h.query_cov < min_coverage and h.percent_identity < min_identity] for g in gene_names}
    else:
        hits_per_gene = {g: align_query_to_ref(allele_paths[g], assembly_path,
                                               ref_index=minimap2_index, min_identity=min_identity,
                                               min_query_coverage=min_coverage) for g in gene_names}
        spurious_hits = None


    if spurious_hits is not None:
        spurious_hits = {g: process_spurious_hits(spurious_hits[g]) for g in gene_names}


    hits_by_contig = cluster_hits_by_contig(hits_per_gene, gene_names)
    full_set_contigs = find_full_set_contigs(hits_by_contig)

    # If zero or one contigs have the full set of genes, then this is treated as a non-multi-MLST
    # case, i.e. the same as regular MLST.

    if len(full_set_contigs) < 2:
        return run_single_mlst(profiles, hits_per_gene, gene_names, required_exact_matches,
                               check_for_truncation, report_incomplete), spurious_hits

    # If more than one contig has the full set of genes, then this is treated as a multi-MLST case,
    # where each full-set contig gets an MLST call.
    contig_results = {}
    for contig in full_set_contigs:
        contig_results[contig] = run_single_mlst(profiles, hits_by_contig[contig], gene_names,
                                                 required_exact_matches, check_for_truncation,
                                                 report_incomplete)
         
    return combine_results(full_set_contigs, contig_results, gene_names), spurious_hits


# def multi_mlst(assembly_path, minimap2_index, profiles_path, allele_paths, gene_names, extra_info,
#                min_identity, min_coverage, required_exact_matches, check_for_truncation=False,
#                report_incomplete=False):
#     """
#     This function takes and returns the same things as the mlst function in mlst.py. However, it
#     will look for cases where multiple contigs have hits for the full set of MLST genes, and in
#     that case, MLST is run on each of them. Otherwise, it behaves like normal MLST.
#     """
#     profiles = load_st_profiles(profiles_path, gene_names, extra_info)
#     hits_per_gene = {g: align_query_to_ref(allele_paths[g], assembly_path,
#                                            ref_index=minimap2_index, min_identity=min_identity,
#                                            min_query_coverage=min_coverage) for g in gene_names}
#     hits_by_contig = cluster_hits_by_contig(hits_per_gene, gene_names)
#     full_set_contigs = find_full_set_contigs(hits_by_contig)

#     # If zero or one contigs have the full set of genes, then this is treated as a non-multi-MLST
#     # case, i.e. the same as regular MLST.
#     if len(full_set_contigs) < 2:
#         return run_single_mlst(profiles, hits_per_gene, gene_names, required_exact_matches,
#                                check_for_truncation, report_incomplete)

#     # If more than one contig has the full set of genes, then this is treated as a multi-MLST case,
#     # where each full-set contig gets an MLST call.
#     contig_results = {}
#     for contig in full_set_contigs:
#         contig_results[contig] = run_single_mlst(profiles, hits_by_contig[contig], gene_names,
#                                                  required_exact_matches, check_for_truncation,
#                                                  report_incomplete)
#     return combine_results(full_set_contigs, contig_results, gene_names)


def cluster_hits_by_contig(hits_per_gene, gene_names):
    """
    Takes a dictionary of hits per gene (key: gene name, value: list of hits) and returns the same
    hits grouped by the contig name (key: contig name, value: gene name to list of hits dict).
    """
    all_contig_names = set()
    for gene, hits in hits_per_gene.items():
        all_contig_names.update([h.ref_name for h in hits])

    hits_by_contig = {}
    for contig in sorted(all_contig_names):
        hits_by_contig[contig] = {g: [] for g in gene_names}

    for gene, hits in hits_per_gene.items():
        for h in hits:
            hits_by_contig[h.ref_name][gene].append(h)
    return hits_by_contig


def find_full_set_contigs(hits_by_contig):
    """
    Returns the names of contigs which have a full set of MLST genes. Contigs are sorted by length
    from big to small (so chromosomal contigs come before plasmid contigs in a completed assembly).
    """
    contig_lengths = {}
    full_set_contigs = []
    for contig, hits_per_gene in hits_by_contig.items():
        if all(len(hits) > 0 for hits in hits_per_gene.values()):
            full_set_contigs.append(contig)
        for hits in hits_per_gene.values():
            for h in hits:
                contig_lengths[h.ref_name] = h.ref_length
    return sorted(full_set_contigs, key=lambda c: contig_lengths[c], reverse=True)


def combine_results(full_set_contigs, contig_results, gene_names):
    """
    Takes per-contig MLST results (ST, extra info and alleles) and combines them into a single
    comma-delimited list.
    """
    combined_st, combined_extra_info = [], []
    combined_alleles = {g: [] for g in gene_names}
    for contig in full_set_contigs:
        st, extra_info, alleles = contig_results[contig]
        combined_st.append(st)
        combined_extra_info.append(extra_info)
        for g in gene_names:
            combined_alleles[g].append(alleles[g])
    return ','.join(combined_st), ','.join(combined_extra_info), \
           {g: ','.join(combined_alleles[g]) for g in gene_names}


def get_allele_and_locus(hit):
    """
    Parses the allele name and locus name from the hit's gene ID.
    """
    if '__' in hit.query_name:  # srst2 formatted file
        gene_id_components = hit.query_name.split('__')
        locus = gene_id_components[1]
        allele = gene_id_components[2]
    else:
        allele = hit.query_name
        locus = hit.query_name.split('_')[0]
    return allele, locus

def process_spurious_hits(hits):
    hit_strings = []
    for hit in hits:
        allele, locus = get_allele_and_locus(hit)
        alignment_length = hit.query_end - hit.query_start
        if hit.percent_identity < 100.00 or alignment_length < hit.query_length:
            allele += '*'  # inexact match
        allele += truncation_check(hit)[0]
        hit_strings.append(allele)
    return hit_strings
