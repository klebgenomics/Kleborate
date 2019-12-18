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


def run_blastn(db, query, min_cov, min_ident, culling_limit=1, ungapped=False):
    build_blast_database_if_needed(db)

    cmd = 'blastn -task blastn -db {} -query {}'.format(db, query)
    cmd += " -outfmt '6 sacc pident slen length bitscore qseq sstrand sstart send" \
           " qacc qstart qend qframe'"
    cmd += ' -dust no -evalue 1E-20 -word_size 32 -max_target_seqs 10000'
    cmd += ' -perc_identity {}'.format(min_ident)
    if ungapped:
        cmd += ' -ungapped'

    # TODO: switch this over to subprocess
    blast_hits = []
    f = os.popen(cmd)
    for line in f:
        blast_hits.append(BlastHit(line))
    f.close()

    # Toss out low identity and low coverage hits.
    if min_ident is not None:
        blast_hits = [h for h in blast_hits if h.pcid * 100 >= min_ident]
    if min_cov is not None:
        blast_hits = [h for h in blast_hits if h.ref_cov * 100 >= min_cov]

    return cull_redundant_hits(blast_hits)


def cull_redundant_hits(blast_hits):
    """
    Cull out redundant hits here (essentially implementing BLAST's -culling_limit 1 feature but
    with our own logic).
    """
    # Sort the hits from best to worst. Hit quality is defined as the product of gene coverage and
    # identity.
    # blast_hits = sorted(blast_hits, key=lambda x: (x.ref_cov * x.pcid), reverse=True)
    # blast_hits = sorted(blast_hits, key=lambda x: (1/x.score, x.gene_id))
    # blast_hits = sorted(blast_hits, key=lambda x: (1/x.pcid, 1/x.ref_hit_len, x.gene_id))
    blast_hits = sorted(blast_hits, key=lambda x: (1/(x.pcid * x.score), x.gene_id))

    filtered_blast_hits = []

    for h in blast_hits:
        if not overlapping(h, filtered_blast_hits):
            filtered_blast_hits.append(h)

    return filtered_blast_hits


def overlapping(hit, existing_hits):
    # Only consider hits in the same reading frame.
    existing_hits = [h for h in existing_hits if
                     h.strand == hit.strand and h.frame == hit.frame and
                     h.contig_name == hit.contig_name]

    for existing_hit in existing_hits:
        if hits_overlap(hit, existing_hit):
            return True

    return False


def hits_overlap(a, b):
    return a.contig_start <= b.contig_end and b.contig_start <= a.contig_end


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
        self.contig_name = fields[9]            # qacc
        self.contig_start = int(fields[10])     # qstart
        self.contig_end = int(fields[11])       # qend
        self.frame = fields[12]                 # qframe

        if self.strand == 'plus':
            assert self.ref_end >= self.ref_start
            self.ref_hit_len = self.ref_end - self.ref_start + 1
        else:
            assert self.strand == 'minus'
            assert self.ref_start >= self.ref_end
            self.ref_hit_len = self.ref_start - self.ref_end + 1
        assert self.contig_end >= self.contig_start

        self.ref_cov = self.ref_hit_len / self.ref_length


def build_blast_database_if_needed(seqs):
    if not os.path.exists(seqs + '.nin'):
        with open(os.devnull, 'w') as devnull:
            subprocess.check_call('makeblastdb -dbtype nucl -in ' + seqs, stdout=devnull,
                                  shell=True)
