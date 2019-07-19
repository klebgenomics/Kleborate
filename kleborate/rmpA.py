"""
Reports best match
  * :    best matching allele is not precise match
  -nLV : best matching ST is n-locus variant

If an annotation column is provided (such as clonal complex) in the final column of the profiles
file, this annotation will be reported in column 2 of the output table

NOTE there is a bug with the culling_limit parameter in older versions of BLAST+. This code has
been tested with BLAST+2.2.30. It does not work with BLAST2.2.25. Not sure about other versions.

Copyright 2017 Kat Holt
Copyright 2017 Ryan Wick (rrwick@gmail.com)
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
import subprocess

from .truncation import truncation_check


def print_header():
    print('\t'.join(['strain', 'rmpA_allele', 'rmpA_lineage', 'rmpA2_allele']))


def rmpa_blast(seqs, database, assemblies, minident):
    if not os.path.exists(seqs + '.nin'):
        with open(os.devnull, 'w') as devnull:
            subprocess.check_call('makeblastdb -dbtype nucl -in ' + seqs, stdout=devnull,
                                  shell=True)

    # read in rmpA database
    st_info = {}  # key = st, value = info relating to this ST, eg clonal group
    header = []
    with open(database, 'r') as f:
        for line in f:
            fields = line.rstrip().split('\t')
            if len(header) == 0:
                header = fields
            else:
                st_info[fields[0]] = fields[1]

    # search input assemblies
    for contigs in assemblies:
        _, file_name = os.path.split(contigs)
        name, ext = os.path.splitext(file_name)

        # blast against all rmpA and rmpA2 alleles
        f = os.popen('blastn -task blastn -db ' + seqs + ' -query ' + contigs +
                     " -outfmt '6 sacc pident slen length score qseq sstrand sstart send' " +
                     '-dust no -evalue 1E-20 -word_size 32 -max_target_seqs 10000 ' +
                     '-culling_limit 1 -perc_identity ' + str(minident))

        rmpa_calls, rmpa2_calls = [], []
        for line in f:
            fields = line.rstrip().split('\t')
            gene_id, pcid, ref_length, length, score, hit_seq, strand, ref_start, ref_end = \
                fields[0], float(fields[1]), int(fields[2]), int(fields[3]), float(fields[4]),\
                fields[5], fields[6], int(fields[7]), int(fields[8])

            if length > (ref_length / 2):
                if pcid < 100.00 or ref_length < length:
                    gene_id += '*'


                gene_id += truncation_check(hit_seq, strand, ref_start, ref_end, ref_length)

                gene = gene_id.split('_')[0]
                if gene == 'rmpA':
                    info = '(' + st_info[gene_id.split('_')[1]] + ')'  # predict from best hit
                    rmpa_calls.append(gene_id + info)
                else:
                    rmpa2_calls.append(gene_id)
        f.close()

        if len(rmpa_calls) == 0:
            rmpa_calls.append('-')
        if len(rmpa2_calls) == 0:
            rmpa2_calls.append('-')

        return [name, ','.join(rmpa_calls), ','.join(rmpa2_calls)]
