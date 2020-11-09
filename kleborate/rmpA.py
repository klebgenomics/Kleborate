"""
Reports best match
  * :    best matching allele is not precise match
  -nLV : best matching ST is n-locus variant

If an annotation column is provided (such as clonal complex) in the final column of the profiles
file, this annotation will be reported in column 2 of the output table

NOTE there is a bug with the culling_limit parameter in older versions of BLAST+. This code has
been tested with BLAST+2.2.30. It does not work with BLAST2.2.25. Not sure about other versions.

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

from .blastn import run_blastn
from .truncation import truncation_check


def rmpa2_blast(seqs, assemblies, min_cov, min_ident):
    for contigs in assemblies:
        rmpa2_calls = []
        hits = run_blastn(seqs, contigs, min_cov, min_ident)
        for hit in hits:
            if hit.alignment_length > hit.ref_length / 2:
                gene_id = hit.gene_id
                if hit.pcid < 100.00 or hit.alignment_length < hit.ref_length:
                    gene_id += '*'
                gene_id += truncation_check(hit)[0]
                rmpa2_calls.append(gene_id)
        if len(rmpa2_calls) == 0:
            rmpa2_calls.append('-')
        return ','.join(rmpa2_calls)
