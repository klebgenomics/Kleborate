#!/usr/bin/env python3
"""
Reports best match
  * :   best matching allele is not precise match
 -nLV : best matching ST is n-locus variant

If an annotation column is provided (such as clonal complex) in the final column of the profiles
file, this annotation will be reported in column 2 of the output table.

NOTE there is a bug with the culling_limit parameter in older versions of BLAST+. This code has
been tested with BLAST+2.2.30. It does not work with BLAST2.2.25. Not sure about other versions.

Copyright 2018 Kat Holt
Copyright 2018 Ryan Wick (rrwick@gmail.com)
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

from .blastn import run_blastn


def mlst_blast(seqs, database, info_arg, assemblies, min_ident, maxmissing, print_header):
    sts = {}  # key = concatenated string of alleles, value = st
    st_info = {}  # key = st, value = info relating to this ST, eg clonal group
    max_st = 0  # holds the highest current ST, incremented when novel combinations are encountered
    header = []
    info_title = 'info'
    with open(database, 'r') as f:
        for line in f:
            fields = line.rstrip().split('\t')
            if len(header) == 0:
                header = fields
                header.pop(0)  # remove st label
                if info_arg == 'yes':
                    info_title = header.pop()  # remove info label
            else:
                st = fields.pop(0)
                if info_arg == 'yes':
                    info = fields.pop()
                else:
                    info = ''
                sts[','.join(fields)] = st
                if int(st) > max_st:
                    max_st = int(st)
                if info_arg == 'yes':
                    st_info[st] = info

    # In order to call an ST, there needs to be an exact match for half (rounded down) of the
    # relevant alleles.
    required_exact_matches = int(len(header) / 2)

    if print_header:
        if info_arg == 'yes':
            print('\t'.join(['strain', info_title, 'ST'] + header))
        else:
            print('\t'.join(['strain', 'ST'] + header))

    contigs = assemblies[0]
    _, filename = os.path.split(contigs)
    name, _ = os.path.splitext(filename)

    best_score = {}   # key = locus, value = BLAST score for best allele encountered so far
    best_allele = {}  # key = locus, value = best allele (* if imprecise match)

    hits = run_blastn(seqs, contigs, None, min_ident, ungapped=True, culling_limit=2)
    for hit in hits:
        if '__' in hit.gene_id:  # srst2 formatted file
            gene_id_components = hit.gene_id.split('__')
            locus = gene_id_components[1]
            allele = gene_id_components[2]
        else:
            allele = hit.gene_id
            locus = hit.gene_id.split('_')[0]
        if hit.pcid < 100.00 or hit.alignment_length < hit.ref_length:
            allele += '*'  # inexact match
        # store best match for each one locus
        if locus in best_score:
            if hit.score > best_score[locus]:    # update
                best_score[locus] = hit.score
                best_allele[locus] = allele.split('_')[1]  # store number only
        else:  # initialise
            best_score[locus] = hit.score
            best_allele[locus] = allele.split('_')[1]  # store number only

    best_st = []
    best_st_annotated = []

    mismatch_loci_including_snps = 0

    for locus in header:
        if locus in best_allele:
            allele = best_allele[locus]
            allele_number = allele.replace('*', '')
            if allele.endswith('*'):
                mismatch_loci_including_snps += 1
            best_st.append(allele_number)
            best_st_annotated.append(allele)  # will still have character if imperfect match
        else:
            best_st.append('-')
            best_st_annotated.append('-')
            mismatch_loci_including_snps += 1

    # assign ST
    bst = ','.join(best_st)

    if mismatch_loci_including_snps <= maxmissing:
        # only report ST if enough loci are precise matches
        if bst in sts:
            # note may have mismatching alleles due to SNPs, this will be recorded in
            # mismatch_loci_including_snps
            bst = sts[bst]
        else:
            # determine closest ST
            bst, _, mismatch_loci_including_snps = \
                get_closest_locus_variant(best_st, best_st_annotated, sts)
    else:
        bst = '0'

    exact_matches = len(best_st) - mismatch_loci_including_snps
    if exact_matches < required_exact_matches:
        bst = '0'

    # pull info column
    info_final = ''
    if info_arg == 'yes':
        if bst in st_info:
            info_final = st_info[bst]

    if mismatch_loci_including_snps > 0 and bst != '0':
        bst += '-' + str(mismatch_loci_including_snps) + 'LV'

    if info_arg == 'yes':
        return [name, info_final, bst] + best_st_annotated
    else:
        return [name, bst] + best_st_annotated


def get_closest_locus_variant(query, annotated_query, sts):
    annotated_query = list(annotated_query)  # copy the list so we don't change the original
    closest = []
    closest_alleles = {}   # key = st, value = list
    min_dist = len(query)  # number mismatching loci, ignoring SNPs

    for index, item in enumerate(query):
        if item == '-':
            query[index] = '0'

    # get distance from closest ST, ignoring SNPs
    for st in sts:
        d = sum(map(lambda x, y: bool(int(x)-int(y)), st.split(','), query))
        if d == min_dist:
            closest.append(int(sts[st]))
            closest_alleles[sts[st]] = st
        elif d < min_dist:
            # reset
            closest = [int(sts[st])]
            closest_alleles[sts[st]] = st
            min_dist = d  # distance from closest ST, ignoring SNPs

    closest_st = str(min(closest))

    for index, item in enumerate(annotated_query):
        if item == '-' or item.endswith('*'):
            annotated_query[index] = '0'

    # get distance from closest ST, including SNPs
    min_dist_incl_snps = sum(map(lambda x, y: bool(int(x)-int(y)),
                                 closest_alleles[closest_st].split(','), annotated_query))

    return closest_st, min_dist, min_dist_incl_snps
