"""
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
import subprocess


def run_blastn(db, query, minident, culling_limit=1, ungapped=False):
    build_blast_database_if_needed(db)

    cmd = 'blastn -task blastn -db {} -query {}'.format(db, query)
    cmd += " -outfmt '6 sacc pident slen length score qseq sstrand sstart send'"
    cmd += ' -dust no -evalue 1E-20 -word_size 32 -max_target_seqs 10000'
    cmd += ' -culling_limit {} -perc_identity {}'.format(culling_limit, minident)
    if ungapped:
        cmd += ' -ungapped'

    # TODO: switch this over to subprocess
    blast_hits = []
    f = os.popen(cmd)
    for line in f:
        blast_hits.append(BlastHit(line))
    f.close()

    return blast_hits


class BlastHit(object):
    def __init__(self, line):
        fields = line.rstrip().split('\t')
        self.gene_id = fields[0]                # sacc
        self.pcid = float(fields[1])            # pident
        self.ref_length = int(fields[2])        # slen
        self.alignment_length = int(fields[3])  # length
        self.score = float(fields[4])           # score
        self.hit_seq = fields[5]                # qseq
        self.strand = fields[6]                 # sstrand
        self.ref_start = int(fields[7])         # sstart
        self.ref_end = int(fields[8])           # send


def build_blast_database_if_needed(seqs):
    if not os.path.exists(seqs + '.nin'):
        with open(os.devnull, 'w') as devnull:
            subprocess.check_call('makeblastdb -dbtype nucl -in ' + seqs, stdout=devnull,
                                  shell=True)
