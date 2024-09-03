"""
Copyright 2024 Kat Holt
Copyright 2020 Ryan Wick (rrwick@gmail.com)
https://github.com/katholt/Kleborate/

This file is part of Kleborate. Kleborate is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Kleborate is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Kleborate. If
not, see <https://www.gnu.org/licenses/>.
"""

import os
import re
import subprocess
import sys

from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError
from .misc import load_fasta, reverse_complement


class Alignment(object):
    """
    Defines a minimap2 alignment. Each object is created from a single line in a minimap2 PAF file.

    It is assumed that minimap2 was run with its -c option so that full alignment is performed. If
    -c wasn't used, then some pieces will be missing (e.g. CIGAR) and some will be incorrect (e.g.
    percent_identity).

    If dictionaries of the query and reference sequences are also provided (key=name, value=seq),
    then the Alignment object will also contain the
    """

    def __init__(self, paf_line, query_seqs=None, ref_seqs=None):
        self.query_name, self.query_length = None, None
        self.query_start, self.query_end = None, None
        self.strand = None
        self.ref_name, self.ref_length = None, None
        self.ref_start, self.ref_end = None, None
        self.matching_bases, self.num_bases = None, None
        self.percent_identity = None
        self.query_cov, self.ref_cov = None, None
        self.cigar, self.alignment_score = None, None
        self.query_seq, self.ref_seq = None, None

        self.parse_paf_line(paf_line)
        self.set_identity_and_coverages()
        self.set_sequences(query_seqs, ref_seqs)

    def parse_paf_line(self, paf_line):
        line_parts = paf_line.strip().split('\t')
        if len(line_parts) < 11:
            sys.exit('Error: alignment file does not seem to be in PAF format')

        self.query_name = line_parts[0]
        self.query_length = int(line_parts[1])
        self.query_start = int(line_parts[2])
        self.query_end = int(line_parts[3])
        self.strand = line_parts[4]

        self.ref_name = line_parts[5]
        self.ref_length = int(line_parts[6])
        self.ref_start = int(line_parts[7])
        self.ref_end = int(line_parts[8])

        self.matching_bases = int(line_parts[9])
        self.num_bases = int(line_parts[10])

        self.cigar, self.alignment_score = None, None
        for part in line_parts:
            if part.startswith('cg:Z:'):
                self.cigar = part[5:]
            if part.startswith('AS:i:'):
                self.alignment_score = int(part[5:])

    def set_identity_and_coverages(self):
        self.percent_identity = 100.0 * self.matching_bases / self.num_bases
        self.query_cov = 100.0 * (self.query_end - self.query_start) / self.query_length
        self.ref_cov = 100.0 * (self.ref_end - self.ref_start) / self.ref_length

    def set_sequences(self, query_seqs, ref_seqs):
        if query_seqs is not None:
            self.query_seq = query_seqs[self.query_name][self.query_start:self.query_end]
        if ref_seqs is not None:
            self.ref_seq = ref_seqs[self.ref_name][self.ref_start:self.ref_end]
            if self.strand == '-':
                self.ref_seq = reverse_complement(self.ref_seq)

    def __repr__(self):
        return self.query_name + ':' + str(self.query_start) + '-' + str(self.query_end) + \
               '(' + self.strand + '), ' + \
               self.ref_name + ':' + str(self.ref_start) + '-' + str(self.ref_end) + \
               ' (' + ('%.3f' % self.percent_identity) + '%)'

    def get_translated_ref_seq(self):
        nucl_seq = self.ref_seq
        ambiguous_bases = set(b for b in nucl_seq) - {'A', 'C', 'G', 'T'}
        for b in ambiguous_bases:
            nucl_seq = nucl_seq.split(b)[0]  # truncate to first ambiguous base
        nucl_seq = nucl_seq[:len(nucl_seq) // 3 * 3]  # truncate to a multiple of 3
        coding_dna = Seq(nucl_seq)
        return str(coding_dna.translate(table='Bacterial', to_stop=True))

    def is_exact(self):
        """
        Returns True if the alignment covers the entire query with perfect identity.
        """
        return (self.matching_bases == self.num_bases and  # 100% identity
                self.query_end - self.query_start == self.query_length)  # 100% coverage


def align_query_to_ref(query_filename, ref_filename, ref_index=None, preset='map-ont',
                        min_identity=None, min_query_coverage=None):
     """
     Runs minimap2 on two sequence files (FASTA or FASTQ) and returns a list of Alignment objects.
     Optional arguments:
     * ref_index: a minimap2 index for the reference. If provided, this will save a bit of time
                  because minimap2 won't need to make the index.
     * preset: the value for minimap2's preset option (-x)
     * min_identity: if provided, alignments with an identity lower than this are discarded.
                     Expressed as a percentage, so values should be 0-100.
     * min_query_coverage: if provided, alignments with a query coverage lower than this are
                           discarded. Expressed as a percentage, so values should be 0-100.
     """
     query_seqs = dict(load_fasta(query_filename))
     ref_seqs = dict(load_fasta(ref_filename))
     ref = ref_filename if ref_index is None else ref_index
     with open(os.devnull, 'w') as dev_null:
         out = subprocess.check_output(['minimap2','--end-bonus=10','--eqx', '-c', '-x', preset,
                                        str(ref), str(query_filename)], stderr=dev_null)
     alignments = [Alignment(x, query_seqs=query_seqs, ref_seqs=ref_seqs)
                   for x in out.decode().splitlines()]
     if min_identity is not None:
         alignments = [a for a in alignments if a.percent_identity >= min_identity]
     if min_query_coverage is not None:
         alignments = [a for a in alignments if a.query_cov >= min_query_coverage]
     return alignments

def get_expanded_cigar(cigar):
    """
    Takes in a normal CIGAR string and returns an expanded version.
    E.g. 5=1D3= -> =====D===
    """
    expanded_cigar = []
    cigar_parts = re.findall(r'\d+[IDX=M]', cigar)
    for p in cigar_parts:
        size = int(p[:-1])
        letter = p[-1]
        expanded_cigar.append(letter * size)
    return ''.join(expanded_cigar)


def hits_overlap(a, b):
    if a.ref_start <= b.ref_end and b.ref_start <= a.ref_end:  # There is some overlap
        allowed_overlap = 50
        overlap_size = len(range(max(a.ref_start, b.ref_start),
                                 min(a.ref_end, b.ref_end)))
        return overlap_size > allowed_overlap
    else:
        return False


def overlapping(hit, existing_hits):
    existing_hits = [h for h in existing_hits if
                     h.strand == hit.strand and h.ref_name == hit.ref_name]

    for existing_hit in existing_hits:
        if hits_overlap(hit, existing_hit):
            return True

    return False
  

def cull_redundant_hits(minimap_hits):
    
    # Sort the hits from best to worst. Hit quality is defined as the product of gene coverage,identity and score
    
    minimap_hits = sorted(minimap_hits, key=lambda x: (1/(x.percent_identity * x.alignment_score * x.query_cov), x.query_name))

    filtered_minimap_hits = []

    for h in minimap_hits:
        if not overlapping(h, filtered_minimap_hits):
            filtered_minimap_hits.append(h)

    return filtered_minimap_hits


def truncation_check(alignment, cov_threshold=90.0): 
    """
    This function checks to see if a gene alignment is truncated at the amino acid level. It
    assumes that the query sequence is a full coding sequence for a gene and the reference is an
    assembly which may or may not be a complete coding sequence.

    It returns:
    * a string to be appended to the Kleborate result, e.g. '-60%'.
    * the amino acid coverage of the reference sequence, e.g. 60.3.
    """
    # The hit must start at the first base of the gene. If not, the gene is considered 0%.
    if alignment.query_start != 0:
        return '-0%', 0.0,''


    # The assumption is that the reference allele is a full CDS with a stop codon at the end. This
    # isn't always true (the reference sequence is sometimes broken) but will serve to make our
    # denominator for coverage.
    query_aa_length = (alignment.query_length - 3) // 3

    translation = alignment.get_translated_ref_seq()
    coverage = 100.0 * len(translation) / query_aa_length
    
    if coverage >= cov_threshold:
        return '', coverage, translation
    else:
        return '-{:.0f}%'.format(coverage), coverage, translation


def check_for_exact_aa_match(ref_file, hit, contigs):
    
    """
    This function checks to see if an exact amino acid match can be found for a sequence that had
    an inexact nucleotide match. If so, return the gene_id, otherwise None. If multiple references
    have exact amino acid matches, it returns the longest one. If multiple references have
    equally-long exact amino acid matches, it returns the alphabetically first.
    """

    
    # First, we extract the nucleotide sequence from the assembly.
    hit_seq = hit.ref_seq
    assembly_seqs = dict(load_fasta(contigs))
    contig_start, contig_end = hit.ref_start, hit.ref_end  # 0-based indexing
    contig_length = len(assembly_seqs[hit.ref_name])
    gene_nucl_seq = assembly_seqs[hit.ref_name][contig_start:contig_end]
    if hit.strand == '-':
         gene_nucl_seq = reverse_complement(gene_nucl_seq)
    assert hit_seq == gene_nucl_seq
    
    # We also need to check whether the first few or last few bases of the sequence is missing.
    # This is to catch cases where an alternative start/stop codon can lead to an incomplete
    # nucleotide match even when there is an exact amino acid match. If we find that the hit is
    # missing start or end bases (relative to the reference), then we add those back on and will
    # include this augmented sequence in the exact amino acid check.
    
    ref_seqs = load_fasta(ref_file)
    ref_length = len(dict(ref_seqs)[hit.query_name])
    ref_start, ref_end = sorted([hit.query_start, hit.query_end])
    missing_start = ref_start
    missing_end = ref_length - ref_end
    if missing_start == 0 and missing_end == 0:
        augmented_gene_nucl_seq = None
    elif missing_start > 10 and missing_end > 10:  # don't bother with too much missing start/end
        augmented_gene_nucl_seq = None
    else:
        if hit.strand == '+':
            contig_start -= missing_start
            contig_end += missing_end
        elif hit.strand == '-':
            contig_start -= missing_end
            contig_end += missing_start
        else:
            assert False
        contig_start = max(contig_start, 0)
        contig_end = min(contig_end, contig_length)
        augmented_gene_nucl_seq = assembly_seqs[hit.ref_name][contig_start:contig_end]
        if hit.strand == '-':
            augmented_gene_nucl_seq = reverse_complement(augmented_gene_nucl_seq)
            
    # Look for an amino acid match between the assembly sequence and any reference sequence.
    best_match_length = 0
    best_matches = []
    for name, ref_nucl_seq in ref_seqs:
        match = is_exact_aa_match(gene_nucl_seq, ref_nucl_seq)
        if augmented_gene_nucl_seq is not None and \
                is_exact_aa_match(augmented_gene_nucl_seq, ref_nucl_seq):
            match = True
        if match:
            if len(ref_nucl_seq) > best_match_length:
                best_matches = [name]
                best_match_length = len(ref_nucl_seq)
            elif len(ref_nucl_seq) == best_match_length:
                best_matches.append(name)
    if not best_matches:
        return None
    else:
        return sorted(best_matches)[0]



def is_exact_aa_match(gene_nucl_seq_1, ref_nucl_seq):
    # look at the gene nucleotide sequence in all three frames of the forward strand.
    gene_nucl_seq_2 = gene_nucl_seq_1[1:]
    gene_nucl_seq_3 = gene_nucl_seq_1[2:]

    gene_prot_1 = translate_nucl_to_prot(gene_nucl_seq_1)
    gene_prot_2 = translate_nucl_to_prot(gene_nucl_seq_2)
    gene_prot_3 = translate_nucl_to_prot(gene_nucl_seq_3)
    ref_prot = translate_nucl_to_prot(ref_nucl_seq)

    # If the reference protein sequence is contained within any frame of the gene protein sequence,
    # that counts as a match.
    return (ref_prot in gene_prot_1) or (ref_prot in gene_prot_2) or (ref_prot in gene_prot_3)


def translate_nucl_to_prot(nucl_seq):
    # First try to translate as a complete coding sequence. This will allow for alternative start
    # codons (e.g. GTG -> M) if it works. We have to manually add the stop codon (*) here because
    # using cds=True turns that off.
    try:
        return str(Seq(nucl_seq).translate(table='Bacterial', to_stop=False, cds=True)) + '*'
    except TranslationError:
        pass

    # If that failed, we will translate in a more relaxed way using a nucleotide sequence truncated
    # to a multiple-of-three length.
    truncated_nucl_seq = nucl_seq[:len(nucl_seq) // 3 * 3]
    return str(Seq(truncated_nucl_seq).translate(table='Bacterial', to_stop=False, cds=False))

def get_bases_per_ref_pos(alignment):
    aligned_seq1, aligned_seq2 = alignment[0], alignment[1]
    bases_per_ref_pos = {}
    ref_pos = 1
    for i, ref_b in enumerate(aligned_seq1):
        if ref_b == '-' or ref_b == '.':
            continue
        assembly_b = aligned_seq2[i]
        bases_per_ref_pos[ref_pos] = assembly_b
        ref_pos += 1
    return bases_per_ref_pos
