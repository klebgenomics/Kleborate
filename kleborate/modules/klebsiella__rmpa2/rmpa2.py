"""
Reports best match
  * :    best matching allele is not precise match
  -nLV : best matching ST is n-locus variant

If an annotation column is provided (such as clonal complex) in the final column of the profiles
file, this annotation will be reported in column 2 of the output table

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

from ...shared.alignment import align_query_to_ref, cull_redundant_hits, truncation_check
from ...shared.misc import load_fasta, reverse_complement

def rmpa2_minimap(ref_file, assembly, minimap2_index, min_coverage, min_identity):
    rmpa2_calls = []
    hits = align_query_to_ref(ref_file, assembly,ref_index=minimap2_index,  min_identity=min_identity, min_query_coverage=min_coverage)
    hits = cull_redundant_hits(hits)
    # Get rid of hits that start with 'delete_'
    hits = [h for h in hits if not h.query_name.startswith('delete_')]
    for hit in hits:
        alignment_length = hit.ref_end - hit.ref_start
        if alignment_length > hit.query_length / 2:
            gene_id = hit.query_name
            if hit.percent_identity < 100.00 or alignment_length < hit.query_length:
                gene_id += '*'
            gene_id += truncation_check(hit)[0]
            rmpa2_calls.append(gene_id)
        if len(rmpa2_calls) == 0:
            rmpa2_calls.append('-')
    return ','.join(rmpa2_calls)
