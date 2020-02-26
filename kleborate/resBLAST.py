"""
Blast for resistance genes, summarise by class (one class per column)

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

import collections
import subprocess
import tempfile

from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62

from .blastn import run_blastn
from .truncation import truncation_check


def resblast_one_assembly(contigs, gene_info, qrdr, trunc, omp, seqs, min_cov, min_ident,
                          min_spurious_cov, min_spurious_ident):
    hits_dict = blast_against_all(seqs, min_cov, min_ident, contigs, gene_info,
                                  min_spurious_cov, min_spurious_ident)
    if qrdr:
        check_for_qrdr_mutations(hits_dict, contigs, qrdr, min_ident, 90.0)
    if trunc:
        check_for_mgrb_pmrb_gene_truncations(hits_dict, contigs, trunc, min_ident)
    if omp:
        check_omp_genes(hits_dict, contigs, omp, min_ident, 90.0)
    return hits_dict


def read_class_file(res_class_file):
    gene_info = {}  # key = sequence id (fasta header in seq file), value = (allele,class,Bla_Class)
    res_classes = []
    bla_classes = []

    with open(res_class_file, 'r') as f:
        header = 0
        for line in f:
            if header == 0:
                header = 1
                # clusterid,queryID,class,gene,allele,seqID,accession,positions,size,
                # cluster_contains_multiple_genes,gene_found_in_multiple_clusters,bla_description,
                # bla_class
            else:
                fields = line.rstrip().split(',')

                cluster_id, res_class, gene, allele_symbol, seq_id, bla_class = \
                    fields[0], fields[2], fields[3], fields[4], fields[5], fields[12]
                seq_header = '__'.join([cluster_id, gene + '_' + res_class, allele_symbol, seq_id])

                if res_class == 'Bla' and bla_class == 'NA':
                    bla_class = 'Bla'
                gene_info[seq_header] = (allele_symbol, res_class, bla_class)
                if res_class not in res_classes:
                    res_classes.append(res_class)
                if bla_class not in bla_classes:
                    bla_classes.append(bla_class)

    res_classes.sort()
    if 'Bla' in res_classes:
        res_classes.remove('Bla')
    bla_classes.sort()
    if 'NA' in bla_classes:
        bla_classes.remove('NA')

    if 'Omp' not in res_classes:
        res_classes.append('Omp')

    return gene_info, res_classes, bla_classes


def get_res_headers(res_classes, bla_classes):
    return res_classes + bla_classes


def blast_against_all(seqs, min_cov, min_ident, contigs, gene_info, min_spurious_cov,
                      min_spurious_ident):
    hits_dict = collections.defaultdict(list)  # key = class, value = list
    hits = run_blastn(seqs, contigs, min_spurious_cov, min_spurious_ident)
    for hit in hits:
        coverage = hit.alignment_length / hit.ref_length * 100.0
        if coverage >= min_spurious_cov:
            if hit.pcid < 100.0:
                aa_result = check_for_exact_aa_match(seqs, hit.hit_seq)
                if aa_result is not None:
                    hit.gene_id = aa_result
            else:
                aa_result = None

            hit_allele, hit_class, hit_bla_class = gene_info[hit.gene_id]
            if hit_class == 'Bla':
                hit_class = hit_bla_class

            if aa_result is not None:
                hit_allele += '^'
            else:
                if hit.pcid < 100.0:
                    hit_allele += '*'
                if hit.alignment_length < hit.ref_length:
                    hit_allele += '?'
                hit_allele += truncation_check(hit)[0]

            # If the hit is decent (above the min coverage and identity thresholds), it goes in the
            # column for the class.
            if coverage >= min_cov and hit.pcid >= min_ident:
                hits_dict[hit_class].append(hit_allele)

            # If the hit is bad (below the min coverage and identity thresholds but above the
            # thresholds for spurious hits) then it goes in the spurious hit column.
            else:
                hits_dict['spurious_resistance_hits'].append(hit_allele)

    return hits_dict


def check_for_exact_aa_match(seqs, gene_nucl_seq):
    """
    This function checks to see if an exact amino acid match can be found for a sequence that had
    an inexact nucleotide match. If so, return the gene_id, otherwise None.
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        # Save just the query allele in a temporary FASTA file.
        q_filename = tmp_dir + '/query.fasta'
        with open(q_filename, 'wt') as q_fasta:
            q_fasta.write('>query\n')
            q_fasta.write(gene_nucl_seq)
            q_fasta.write('\n')

        # tblastx: translated query to translated database
        tblastx_cmd = 'tblastx -db ' + seqs + ' -query ' + q_filename + ' -query_gencode 11' + \
                      " -db_gencode 11 -outfmt '6 bitscore sacc pident slen length sstart send'"
        process = subprocess.Popen(tblastx_cmd, stdout=subprocess.PIPE, stderr=None, shell=True)
        blast_output = process.communicate()[0]
        if not isinstance(blast_output, str):
            blast_output = blast_output.decode()

        # Load in the BLAST results.
        blast_results = []
        for line in blast_output.splitlines():
            p = line.strip().split('\t')
            bitscore, gene_id, pcid, length, allele_length, allele_start, allele_end = \
                float(p[0]), p[1], float(p[2]), float(p[3]), float(p[4]), int(p[5]), int(p[6])
            cov = 100.0 * allele_length * 3.0 / length

            # Only keep perfect matches (identity and coverage of 100%) on the forward ref strand.
            if pcid == 100.0 and cov == 100.0 and allele_start < allele_end:
                blast_results.append((pcid, cov, bitscore, gene_id, allele_length, allele_start,
                                      allele_end))
        if not blast_results:
            return None

        # If there are still multiple matches, we break ties first with bitscore (higher is better)
        # and then with gene ID (alphabetically first is better).
        blast_results = sorted(blast_results, key=lambda r: ((1.0 / r[2]), r[3]))

        best_match = blast_results[0]
        gene_id = best_match[3]
        return gene_id


def check_for_qrdr_mutations(hits_dict, contigs, qrdr, min_ident, min_cov):
    qrdr_loci = {'GyrA': [(83, 'S'), (87, 'D')],
                 'ParC': [(80, 'S'), (84, 'E')]}

    gyra_ref = 'MSDLAREITPVNIEEELKNSYLDYAMSVIVGRALPDVRDGLKPVHRRVLYAMNVLGNDWN' \
               'KAYKKSARVVGDVIGKYHPHGDSAVYDTIVRMAQPFSLRYMLVDGQGNFGSIDGDSAAAM'
    parc_ref = 'MSDMAERLALHEFTENAYLNYSMYVIMDRALPFIGDGLKPVQRRIVYAMSELGLNASAKF' \
               'KKSARTVGDVLGKYHPHGDSACYEAMVLMAQPFSYRYPLVDGQGNWGAPDDPKSFAAMRY'

    snps = []
    hits = run_blastn(qrdr, contigs, None, min_ident)
    for hit in hits:
        _, coverage, translation = truncation_check(hit)
        if coverage > min_cov:
            if hit.gene_id == 'GyrA':
                alignments = pairwise2.align.globalds(gyra_ref, translation, blosum62, -10, -0.5)
            elif hit.gene_id == 'ParC':
                alignments = pairwise2.align.globalds(parc_ref, translation, blosum62, -10, -0.5)
            else:
                assert False
            bases_per_ref_pos = get_bases_per_ref_pos(alignments[0])
            loci = qrdr_loci[hit.gene_id]
            for pos, wt_base in loci:
                assembly_base = bases_per_ref_pos[pos]
                if pos in bases_per_ref_pos and assembly_base != wt_base and \
                        assembly_base != '-' and assembly_base != '.':
                    snps.append(hit.gene_id + '-' + str(pos) + assembly_base)
    if snps:
        hits_dict['Flq'] += snps


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


def check_for_mgrb_pmrb_gene_truncations(hits_dict, contigs, seqs, min_ident):
    best_mgrb_cov, best_pmrb_cov = 0.0, 0.0

    hits = run_blastn(seqs, contigs, None, min_ident)
    for hit in hits:
        assert hit.gene_id == 'pmrB' or hit.gene_id == 'mgrB'
        _, coverage, _ = truncation_check(hit)
        if hit.gene_id == 'mgrB' and coverage > best_mgrb_cov:
            best_mgrb_cov = coverage
        elif hit.gene_id == 'pmrB' and coverage > best_pmrb_cov:
            best_pmrb_cov = coverage

    truncations = []
    if best_mgrb_cov < 90.0:
        truncations.append('MgrB-' + ('%.0f' % best_mgrb_cov) + '%')
    if best_pmrb_cov < 90.0:
        truncations.append('PmrB-' + ('%.0f' % best_pmrb_cov) + '%')

    if truncations:
        if 'Col' not in hits_dict:
            hits_dict['Col'] = []
        hits_dict['Col'] += truncations


def check_omp_genes(hits_dict, contigs, omp, min_ident, min_cov):
    best_ompk35_cov, best_ompk36_cov = 0.0, 0.0

    hits = run_blastn(omp, contigs, None, min_ident)
    for hit in hits:
        _, coverage, translation = truncation_check(hit)
        if hit.gene_id == 'OmpK35':
            if coverage > best_ompk35_cov:
                best_ompk35_cov = coverage
        elif hit.gene_id == 'OmpK36':
            if coverage > best_ompk36_cov:
                best_ompk36_cov = coverage
            if coverage >= min_cov:
                if 'GDGDTY' in translation:
                    hits_dict['Omp'].append('OmpK36GD')
                elif 'GDTDTY' in translation:
                    hits_dict['Omp'].append('OmpK36TD')
        else:
            assert False

    truncations = []
    if best_ompk35_cov < 90.0:
        truncations.append('OmpK35-' + ('%.0f' % best_ompk35_cov) + '%')
    if best_ompk36_cov < 90.0:
        truncations.append('OmpK36-' + ('%.0f' % best_ompk36_cov) + '%')

    if truncations:
        if 'Omp' not in hits_dict:
            hits_dict['Omp'] = []
        hits_dict['Omp'] += truncations
