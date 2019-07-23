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
import os
import subprocess
import tempfile
import xml.etree.ElementTree as ElementTree

from .blastn import run_blastn
from .truncation import truncation_check


def resblast_one_assembly(contigs, gene_info, qrdr, trunc, omp, seqs, mincov, minident):
    build_blast_databases(qrdr, omp)
    hits_dict = blast_against_all(seqs, mincov, minident, contigs, gene_info)
    if qrdr:
        check_for_qrdr_mutations(hits_dict, contigs, qrdr)
    if trunc:
        check_for_mgrb_pmrb_gene_truncations(hits_dict, contigs, trunc, 95.0)
    if omp:
        check_omp_genes(hits_dict, contigs, omp)
    return hits_dict


def build_blast_databases(qrdr, omp):
    if qrdr:
        if not os.path.exists(qrdr + '.pin'):
            with open(os.devnull, 'w') as devnull:
                subprocess.check_call('makeblastdb -dbtype prot -in ' + qrdr,
                                      stdout=devnull, shell=True)
    if omp:
        if not os.path.exists(omp + '.pin'):
            with open(os.devnull, 'w') as devnull:
                subprocess.check_call('makeblastdb -dbtype prot -in ' + omp,
                                      stdout=devnull, shell=True)


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


# functions for finding snps
def get_gapped_position(seq, position):
    num_chars = 0
    i = 0
    while num_chars <= position and i < len(seq):
        if seq[i] != '-':
            num_chars += 1
        i += 1
    return i - 1


def print_header(res_classes, bla_classes):
    print('\t'.join(['strain'] + get_res_headers(res_classes, bla_classes)))


def get_res_headers(res_classes, bla_classes):
    return res_classes + bla_classes


def blast_against_all(seqs, mincov, minident, contigs, gene_info):
    hits_dict = collections.defaultdict(list)  # key = class, value = list
    hits = run_blastn(seqs, contigs, minident, ungapped=True)
    for hit in hits:
        if (hit.alignment_length / hit.ref_length * 100.0) > mincov:
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

            hits_dict[hit_class].append(hit_allele)

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


def blastx_results_as_xml_tree(database, query):
    blastx_cmd = 'blastx -db ' + database + ' -query ' + query + ' -query_gencode 11' + \
                 ' -outfmt 5 -comp_based_stats F -culling_limit 1 -max_hsps 1 -seg no'
    process = subprocess.Popen(blastx_cmd, stdout=subprocess.PIPE, stderr=None, shell=True)
    blast_output = process.communicate()[0]
    if not isinstance(blast_output, str):
        blast_output = blast_output.decode()
    return ElementTree.fromstring(blast_output)


def check_for_qrdr_mutations(hits_dict, contigs, qrdr):
    qrdr_loci = {'GyrA': [(83, 'S'), (87, 'D')],
                 'ParC': [(80, 'S'), (84, 'E')]}

    # key = (locus, pos), value = allele,
    # if found in a simple hit starting at position 1 of the protein seq
    complete_hits, incomplete_hits = collections.defaultdict(list), collections.defaultdict(list)

    root = blastx_results_as_xml_tree(qrdr, contigs)
    for query in root.find('BlastOutput_iterations'):
        for hit in query.find('Iteration_hits'):
            gene_id = hit.find('Hit_def').text
            for hsp in hit.find('Hit_hsps'):
                hsp_hit_eval = float(hsp.find('Hsp_evalue').text)
                hsp_hit_from = int(hsp.find('Hsp_hit-from').text)
                hsp_hit_to = int(hsp.find('Hsp_hit-to').text)
                hsp_gaps = int(hsp.find('Hsp_gaps').text)
                hsp_qseq = hsp.find('Hsp_qseq').text
                hsp_hseq = hsp.find('Hsp_hseq').text

                for pos, wt in qrdr_loci[gene_id]:
                    if hsp_hit_to >= pos and hsp_gaps == 0 and hsp_hit_from == 1:
                        # simple alignment
                        complete_hits[(gene_id, pos)].append(hsp_qseq[pos - 1])
                    else:
                        # not a simple alignment, need to align query and hit and extract loci
                        # manually
                        if hsp_hit_from <= pos <= hsp_hit_to and hsp_hit_eval <= 0.00001:
                            # locus is within aligned area, set evalue to filter out the junk
                            # alignments
                            pos_in_aln = get_gapped_position(hsp_hseq, pos - hsp_hit_from + 1)
                            incomplete_hits[(gene_id, pos)].append(hsp_qseq[pos_in_aln - 1])
    snps = []

    for locus in qrdr_loci:
        for pos, wt in qrdr_loci[locus]:
            if (locus, pos) in complete_hits:
                if complete_hits[(locus, pos)][0] != wt:
                    snps.append(locus + '-' + str(pos) +
                                complete_hits[(locus, pos)][0])  # record SNP at this site
            else:
                if (locus, pos) in incomplete_hits:
                    if incomplete_hits[(locus, pos)][0] != wt:
                        snps.append(locus + '-' + str(pos) +
                                    incomplete_hits[(locus, pos)][0])  # record SNP at this site
    if snps:
        hits_dict['Flq'] += snps


def check_for_mgrb_pmrb_gene_truncations(hits_dict, contigs, seqs, minident):
    best_mgrb_cov, best_pmrb_cov = 0.0, 0.0

    hits = run_blastn(seqs, contigs, minident)
    for hit in hits:
        assert hit.gene_id == 'pmrB' or hit.gene_id == 'mgrB'
        _, coverage = truncation_check(hit)
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


def check_omp_genes(hits_dict, contigs, omp):
    check_for_omp_gene_truncations(hits_dict, contigs, omp)
    check_for_ompk36_mutations(hits_dict, contigs, omp)


def check_for_omp_gene_truncations(hits_dict, contigs, omp):
    best_ompk35_cov, best_ompk36_cov = 0.0, 0.0

    root = blastx_results_as_xml_tree(omp, contigs)
    for query in root.find('BlastOutput_iterations'):
        for hit in query.find('Iteration_hits'):
            gene_id = hit.find('Hit_def').text
            gene_len = int(hit.find('Hit_len').text)
            for hsp in hit.find('Hit_hsps'):
                hsp_qseq = hsp.find('Hsp_qseq').text
                hsp_hseq = hsp.find('Hsp_hseq').text
                identity = sum([1 if a == b else 0
                                for a, b in zip(hsp_qseq, hsp_hseq)]) / len(hsp_hseq)
                hsp_hit_eval = float(hsp.find('Hsp_evalue').text)
                hit_length = max(len(x.replace('-', '')) for x in hsp_qseq.split('*'))
                coverage = 100.0 * float(hit_length) / gene_len
                if hsp_hit_eval <= 0.001 and identity >= 0.9:
                    if gene_id == 'OmpK35' and coverage > best_ompk35_cov:
                        best_ompk35_cov = coverage
                    elif (gene_id == 'OmpK36' or gene_id == 'OmpK36GD' or
                          gene_id == 'OmpK36TD') and coverage > best_ompk36_cov:
                        best_ompk36_cov = coverage

    truncations = []
    if best_ompk35_cov < 90.0:
        truncations.append('OmpK35-' + ('%.0f' % best_ompk35_cov) + '%')
    if best_ompk36_cov < 90.0:
        truncations.append('OmpK36-' + ('%.0f' % best_ompk36_cov) + '%')

    if truncations:
        if 'Omp' not in hits_dict:
            hits_dict['Omp'] = []
        hits_dict['Omp'] += truncations


def check_for_ompk36_mutations(hits_dict, contigs, omp):
    root = blastx_results_as_xml_tree(omp, contigs)
    for query in root.find('BlastOutput_iterations'):
        for hit in query.find('Iteration_hits'):
            gene_id = hit.find('Hit_def').text
            gene_len = int(hit.find('Hit_len').text)
            for hsp in hit.find('Hit_hsps'):
                hsp_qseq = hsp.find('Hsp_qseq').text
                hsp_hseq = hsp.find('Hsp_hseq').text
                identity = sum([1 if a == b else 0
                                for a, b in zip(hsp_qseq, hsp_hseq)]) / len(hsp_hseq)
                hsp_hit_eval = float(hsp.find('Hsp_evalue').text)
                hit_length = max(len(x.replace('-', '')) for x in hsp_qseq.split('*'))
                coverage = 100.0 * float(hit_length) / gene_len
                if coverage >= 90.0 and hsp_hit_eval <= 0.001 and identity >= 0.9:
                    if gene_id == 'OmpK36GD' and 'GDGDTY' in hsp_qseq:
                        hits_dict['Omp'].append('OmpK36GD')
                    if gene_id == 'OmpK36TD' and 'GDTDTY' in hsp_qseq:
                        hits_dict['Omp'].append('OmpK36TD')
