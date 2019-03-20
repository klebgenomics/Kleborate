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

import os
import subprocess
import argparse
import xml.etree.ElementTree as ElementTree


def main():
    args = parse_arguments()
    gene_info, res_classes, bla_classes = read_class_file(args.resclass)
    print_header(res_classes, bla_classes)

    for contigs in args.assemblies:
        hits_dict = resblast_one_assembly(contigs, gene_info, args.qrdr, args.trunc, args.seqs,
                                          args.mincov, args.minident)
        print_results(contigs, res_classes, bla_classes, hits_dict)


def resblast_one_assembly(contigs, gene_info, qrdr, trunc, seqs, mincov, minident):
    build_blast_databases(seqs, qrdr, trunc)
    hits_dict = blast_against_all(seqs, mincov, minident, contigs, gene_info)
    if qrdr:
        check_for_qrdr_mutations(hits_dict, contigs, qrdr)
    if trunc:
        check_for_gene_truncations(hits_dict, contigs, trunc)
    return hits_dict


def parse_arguments():
    parser = argparse.ArgumentParser(description='Klebsiella resistance screen (part of Kleborate)',
                                     add_help=False)

    parser.add_argument('assemblies', nargs='*', type=str,
                        help='FASTA file(s) for assemblies')

    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument('-s', '--seqs', type=str, required=True,
                               help='resistance gene sequences to screen for')
    required_args.add_argument('-t', '--resclass', type=str, required=True,
                               help='resistance gene classes (CSV)')

    additional_screening_args = parser.add_argument_group('Additional screening')
    additional_screening_args.add_argument('-q', '--qrdr', type=str,
                                           help='QRDR sequences for which mutations can cause '
                                                'fluoroquinolone resistance')
    additional_screening_args.add_argument('-r', '--trunc', type=str,
                                           help='MgrB and PmrB genes for which truncation can '
                                                'cause colistin resistance')

    settings_args = parser.add_argument_group('Settings')
    settings_args.add_argument('-m', '--minident', type=float, default=90.0,
                               help='minimum percent identity (default 90)')
    settings_args.add_argument('-c', '--mincov', type=float, default=80.0,
                               help='minimum percent coverage (default 80)')

    help_args = parser.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                           help='show this help message and exit')

    return parser.parse_args()


def build_blast_databases(seqs, qrdr, trunc):
    if not os.path.exists(seqs + '.nin'):
        with open(os.devnull, 'w') as devnull:
            subprocess.check_call('makeblastdb -dbtype nucl -in ' + seqs,
                                  stdout=devnull, shell=True)
    if qrdr:
        if not os.path.exists(qrdr + '.pin'):
            with open(os.devnull, 'w') as devnull:
                subprocess.check_call('makeblastdb -dbtype prot -in ' + qrdr,
                                      stdout=devnull, shell=True)
    if trunc:
        if not os.path.exists(trunc + '.pin'):
            with open(os.devnull, 'w') as devnull:
                subprocess.check_call('makeblastdb -dbtype prot -in ' + trunc,
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
    res_classes.remove('Bla')
    bla_classes.sort()
    bla_classes.remove('NA')

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
    hits_dict = {}  # key = class, value = list

    f = os.popen('blastn -task blastn -db ' + seqs + ' -query ' + contigs +
                 " -outfmt '6 sacc pident slen length score' -ungapped -dust no -evalue 1E-20 " +
                 '-word_size 32 -max_target_seqs 10000 -culling_limit 1 -perc_identity ' +
                 str(minident))
    for line in f:
        fields = line.rstrip().split('\t')
        gene_id, pcid, length, allele_length, score = \
            fields[0], float(fields[1]), float(fields[2]), float(fields[3]), float(fields[4])
        if (allele_length / length * 100) > mincov:
            (hit_allele, hit_class, hit_bla_class) = gene_info[gene_id]
            if hit_class == 'Bla':
                hit_class = hit_bla_class
            if pcid < 100.00:
                hit_allele += settings.inexact_nucleotide_match
            if allele_length < length:
                hit_allele += settings.partial_match
            if hit_class in hits_dict:
                hits_dict[hit_class].append(hit_allele)
            else:
                hits_dict[hit_class] = [hit_allele]
    f.close()
    return hits_dict


def blastx_results_as_xml_tree(database, query):
    blastx_cmd = 'blastx -db ' + database + ' -query ' + query + ' -query_gencode 11' + \
                 ' -outfmt 5 -ungapped -comp_based_stats F -culling_limit 1 -max_hsps 1 -seg no'
    process = subprocess.Popen(blastx_cmd, stdout=subprocess.PIPE, stderr=None, shell=True)
    blast_output = process.communicate()[0]
    if not isinstance(blast_output, str):
        blast_output = blast_output.decode()
    return ElementTree.fromstring(blast_output)


def check_for_qrdr_mutations(hits_dict, contigs, qrdr):
    qrdr_loci = {'GyrA': [(83, 'S'), (87, 'D')], 'ParC': [(80, 'S'), (84, 'E')]}

    # key = (locus, pos), value = allele,
    # if found in a simple hit starting at position 1 of the protein seq
    complete_hits, incomplete_hits = {}, {}

    root = blastx_results_as_xml_tree(qrdr, contigs)
    for query in root[8]:
        for hit in query[4]:
            gene_id = hit[2].text
            for hsp in hit[5]:
                hsp_hit_eval = float(hsp[5].text)
                hsp_hit_from = int(hsp[6].text)
                hsp_hit_to = int(hsp[7].text)
                hsp_gaps = int(hsp[12].text)
                hsp_qseq = hsp[14].text
                hsp_hseq = hsp[15].text

                for (pos, wt) in qrdr_loci[gene_id]:
                    if hsp_hit_to >= pos and (hsp_gaps == 0) and (hsp_hit_from == 1):
                        # simple alignment
                        if (gene_id, pos) in complete_hits:
                            complete_hits[(gene_id, pos)].append(hsp_qseq[pos - 1])
                        else:
                            complete_hits[(gene_id, pos)] = [hsp_qseq[pos - 1]]
                    else:
                        # not a simple alignment, need to align query and hit and extract loci
                        # manually
                        if (pos >= hsp_hit_from) and (pos <= hsp_hit_to) and \
                                (hsp_hit_eval <= 0.00001):
                            # locus is within aligned area, set evalue to filter out the junk
                            # alignments
                            pos_in_aln = get_gapped_position(hsp_hseq, pos - hsp_hit_from + 1)
                            if (gene_id, pos) in incomplete_hits:
                                incomplete_hits[(gene_id, pos)].append(hsp_qseq[pos_in_aln - 1])
                            else:
                                incomplete_hits[(gene_id, pos)] = [hsp_qseq[pos_in_aln - 1]]
    snps = []

    for locus in qrdr_loci:
        for (pos, wt) in qrdr_loci[locus]:
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
        if 'Flq' not in hits_dict:
            hits_dict['Flq'] = []
        hits_dict['Flq'] += snps


def check_for_gene_truncations(hits_dict, contigs, trunc):
    best_mgrb_cov, best_pmrb_cov = 0.0, 0.0

    root = blastx_results_as_xml_tree(trunc, contigs)
    for query in root[8]:
        for hit in query[4]:
            gene_id = hit[2].text
            assert gene_id == 'MgrB' or gene_id == 'PmrB'
            gene_len = int(hit[4].text)

            for hsp in hit[5]:
                hsp_qseq = hsp[14].text
                hit_length = max(len(x) for x in hsp_qseq.split('*'))
                coverage = 100.0 * float(hit_length) / gene_len

                if gene_id == 'MgrB' and coverage > best_mgrb_cov:
                    best_mgrb_cov = coverage
                elif gene_id == 'PmrB' and coverage > best_pmrb_cov:
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


def get_strain_name(full_path):
    filename = os.path.split(full_path)[1]
    if filename.endswith('_temp_decompress.fasta'):
        filename = filename[:-22]
    if filename.endswith('.gz'):
        filename = filename[:-3]
    return os.path.splitext(filename)[0]


def print_results(contigs, res_classes, bla_classes, hits_dict):
    result_string = [get_strain_name(contigs)]
    for res_class in (res_classes + bla_classes):
        if res_class in hits_dict:
            result_string.append(';'.join(hits_dict[res_class]))
        else:
            result_string.append('-')
    print('\t'.join(result_string))


if __name__ == '__main__':
    main()
