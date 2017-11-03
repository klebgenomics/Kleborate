"""
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
import sys
import gzip
import argparse
import subprocess
import distutils.spawn
from pkg_resources import resource_filename
from .contig_stats import load_fasta, get_compression_type, get_contig_stats
from .version import __version__


def main():
    args = parse_arguments()
    check_inputs_and_programs(args)

    # Find necessary resources
    data_folder = resource_filename(__name__, 'data')
    mlstblast = resource_filename(__name__, 'mlstBLAST.py')
    resblast = resource_filename(__name__, 'resBLAST.py')
    clusterblast = resource_filename(__name__, 'clusterBLAST.py')

    # Output in two places: stdout (less verbose) and file (more verbose)
    stdout_header, full_header, res_headers = build_output_headers(args, resblast, data_folder)
    print '\t'.join(stdout_header)
    o = file(args.outfile, 'w')
    o.write('\t'.join(full_header))
    o.write('\n')

    for contigs in args.assemblies:
        (_, filename) = os.path.split(contigs)
        (name, ext) = os.path.splitext(filename)

        # If the contigs are in a gz file, make a temporary decompressed FASTA file.
        if get_compression_type(contigs) == 'gz':
            new_contigs = contigs + '_temp_decompress.fasta'
            decompress_file(contigs, new_contigs)
            contigs = new_contigs
            temp_decompress = True
        else:
            temp_decompress = False

        contig_count, n50, longest_contig = get_contig_stats(contigs)

        # Check the species with Mash
        if args.species:
            species, species_hit_strength = get_klebsiella_species(contigs, data_folder)

        # run chromosome MLST
        f = os.popen('python ' + mlstblast +
                     ' -s ' + data_folder + '/Klebsiella_pneumoniae.fasta' +
                     ' -d ' + data_folder + '/kpneumoniae.txt' +
                     ' -i no' +
                     ' --maxmissing 3' +
                     ' ' + contigs)
        chr_st = ''
        chr_st_detail = []
        for line in f:
            fields = line.rstrip().split('\t')
            if fields[1] != 'ST':  # skip header
                (strain, chr_st) = (fields[0], fields[1])
                chr_st_detail = fields[2:]
                if chr_st != '0':
                    chr_st = 'ST' + chr_st
        f.close()

        # run ybt MLST
        f = os.popen('python ' + mlstblast +
                     ' -s ' + data_folder + '/ybt_alleles.fasta' +
                     ' -d ' + data_folder + '/YbST_profiles.txt' +
                     ' -i yes' +
                     ' --maxmissing 3' +
                     ' ' + contigs)
        yb_st = ''
        yb_group = ''
        yb_st_detail = []
        for line in f:
            fields = line.rstrip().split('\t')
            if fields[2] != 'ST':  # skip header
                strain, yb_st, yb_group = fields[0], fields[2], fields[1]
                yb_st_detail = fields[3:]

                # If no ybt group was found but at least 8 out of the 11 ybt genes are present,
                # then this strain is labelled as 'ybt unknown'.
                if yb_group == '':
                    if sum(0 if x == '-' else 1 for x in yb_st_detail) >= 8:
                        yb_group = 'ybt unknown'
                        yb_st = 'unknown'
                    else:
                        yb_group = '-'
        f.close()

        # run colibactin MLST
        f = os.popen('python ' + mlstblast +
                     ' -s ' + data_folder + '/colibactin_alleles.fasta' +
                     ' -d ' + data_folder + '/CbST_profiles.txt' +
                     ' -i yes' +
                     ' --maxmissing 3' +
                     ' ' + contigs)
        cb_st = ''
        cb_group = ''
        cb_st_detail = []
        for line in f:
            fields = line.rstrip().split('\t')
            if fields[2] != 'ST':  # skip header
                strain, cb_st, cb_group = fields[0], fields[2], fields[1]
                cb_st_detail = fields[3:]

                # If no clb group was found but at least 12 out of the 15 clb genes are present,
                # then this strain is labelled as 'clb unknown'.
                if cb_group == '':
                    if sum(0 if x == '-' else 1 for x in cb_st_detail) >= 12:
                        cb_group = 'clb unknown'
                        cb_st = 'unknown'
                    else:
                        cb_group = '-'
        f.close()

        # screen for other virulence genes (binary calls)
        aerobactin, salmochelin, hypermucoidy = '-', '-', '-'
        f = os.popen('python ' + clusterblast +
                     ' -s ' + data_folder + '/other_vir_clusters.fasta' +
                     ' ' + contigs)
        for line in f:
            fields = line.rstrip().split('\t')
            if fields[1] != 'aerobactin':  # skip header
                aerobactin = fields[1]
                salmochelin = fields[2]
                hypermucoidy = fields[3]
        f.close()

        # screen for wzi allele
        f = os.popen('python ' + mlstblast +
                     ' -s ' + data_folder + '/wzi.fasta' +
                     ' -d ' + data_folder + '/wzi.txt' +
                     ' -i yes' +
                     ' --maxmissing 0' +
                     ' -m 99' +
                     ' ' + contigs)
        for line in f:
            fields = line.rstrip().split('\t')
            if fields[0] != 'ST':  # skip header
                (strain, wzi_st, k_type) = (fields[0], 'wzi' + fields[2], fields[1])
                if fields[2] == '0':
                    wzi_st = '0'
                if k_type == '':
                    k_type = '-'

        # screen for resistance genes
        res_hits = []
        if args.resistance:
            f = os.popen('python ' + resblast +
                         ' -s ' + data_folder + '/ARGannot_r2.fasta' +
                         ' -t ' + data_folder + '/ARGannot_clustered80_r2.csv' +
                         ' -q ' + data_folder + '/QRDR_120.aa' +
                         ' ' + contigs)
            for line in f:
                fields = line.rstrip().split('\t')
                if fields[0] != 'strain':  # skip header
                    res_hits = fields[1:]
            f.close()

        # # TO DO: run Kaptive
        # if args.kaptive:
        #     pass

        # Summarise virulence and resistance.
        virulence_score = get_virulence_score(yb_group, cb_group, aerobactin, salmochelin,
                                              hypermucoidy)
        resistance_score = get_resistance_score(res_headers, res_hits)
        resistance_class_count = get_resistance_class_count(res_headers, res_hits)
        resistance_gene_count = get_resistance_gene_count(res_headers, res_hits)

        # Print results to screen.
        stdout_results = [name]
        if args.species:
            stdout_results += [species]
        stdout_results += [chr_st, str(virulence_score)]
        if args.resistance:
            stdout_results.append(str(resistance_score))
        stdout_results += [yb_group, yb_st, cb_group, cb_st, aerobactin, salmochelin, hypermucoidy,
                           wzi_st, k_type]
        if args.resistance:
            stdout_results += res_hits
        print '\t'.join(stdout_results)

        # Save results to file.
        full_results = [name]
        if args.species:
            full_results += [species, species_hit_strength]
        full_results += [str(contig_count), str(n50), str(longest_contig), chr_st,
                         str(virulence_score)]
        if args.resistance:
            full_results.append(str(resistance_score))
            full_results.append(str(resistance_class_count))
            full_results.append(str(resistance_gene_count))
        full_results += [yb_group, yb_st, cb_group, cb_st, aerobactin, salmochelin, hypermucoidy,
                         wzi_st, k_type, chr_st] + chr_st_detail + yb_st_detail + cb_st_detail
        if args.resistance:
            full_results += res_hits
        o.write('\t'.join(full_results))
        o.write('\n')

        # If we've been working on a temporary decompressed file, delete it now.
        if temp_decompress:
            os.remove(contigs)

    o.close()


def parse_arguments():
    parser = argparse.ArgumentParser(description='Kleborate: a tool for characterising '
                                                 'virulence and resistance in Klebsiella',
                                     add_help=False)

    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument('-a', '--assemblies', nargs='+', type=str, required=True,
                               help='FASTA file(s) for assemblies')

    optional_args = parser.add_argument_group('Optional arguments')
    optional_args.add_argument('-o', '--outfile', type=str, default='Kleborate_results.txt',
                        help='File for detailed output (default: Kleborate_results.txt)')
    optional_args.add_argument('-r', '--resistance', action='store_true',
                        help='Turn on resistance genes screening (default: no resistance gene '
                             'screening)')
    optional_args.add_argument('-s', '--species', action='store_true',
                        help='Turn on Klebsiella species identification (requires Mash, default: '
                             'no species identification)')
    # TO DO
    # optional_args.add_argument('-k', '--kaptive', action='store_true',
    #                            help='Turn on capsule typing with Kaptive (default: no capsule typing')

    help_args = parser.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                           help='Show this help message and exit')
    help_args.add_argument('--version', action='version', version=__version__)

    # If no arguments were used, print the entire help (argparse default is to just give an error
    # like '-a is required').
    if len(sys.argv) == 1:
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    return parser.parse_args()


def check_inputs_and_programs(args):
    for assembly in args.assemblies:
        if not os.path.isfile(assembly):
            sys.exit('Error: could not find ' + assembly)
        fasta = load_fasta(assembly)
        if len(fasta) < 1:
            sys.exit('Error: invalid FASTA file: ' + assembly)
        for record in fasta:
            header, seq = record
            if len(seq) == 0:
                sys.exit('Error: invalid FASTA file (contains a zero-length sequence): ' + assembly)
    if not distutils.spawn.find_executable('blastn'):
        sys.exit('Error: could not find makeblastdb')
    if not distutils.spawn.find_executable('makeblastdb'):
        sys.exit('Error: could not find blastn')
    if args.resistance:
        if not distutils.spawn.find_executable('blastx'):
            sys.exit('Error: could not find blastx')
    if args.species:
        if not distutils.spawn.find_executable('mash'):
            sys.exit('Error: could not find mash')

    # # TO DO
    # if args.kaptive:
    #     try:
    #         imp.find_module('Bio')
    #     except ImportError:
    #         sys.exit('Error: could not find BioPython (required for Kaptive)')


def build_output_headers(args, resblast, data_folder):
    """
    There are two levels of output:
      * stdout is simpler and displayed to the console
      * full contains more and is saved to file
    This function returns headers for both. It also returns the resistance headers in a separate
    list, as they are used to total up some resistance summaries.
    """
    stdout_header = ['strain']
    full_header = ['strain']
    if args.species:
        stdout_header += ['species']
        full_header += ['species', 'species_match']
    stdout_header += ['ST', 'virulence_score']
    full_header += ['contig_count', 'N50', 'largest_contig', 'ST', 'virulence_score']
    if args.resistance:
        stdout_header.append('resistance_score')
        full_header.append('resistance_score')
        full_header.append('num_resistance_classes')
        full_header.append('num_resistance_genes')

    other_columns = ['Yersiniabactin', 'YbST', 'Colibactin', 'CbST', 'aerobactin', 'salmochelin',
                     'hypermucoidy', 'wzi', 'KL']
    stdout_header += other_columns
    full_header += other_columns

    mlst_header = ['Chr_ST', 'gapA', 'infB', 'mdh', 'pgi', 'phoE', 'rpoB', 'tonB', 'ybtS', 'ybtX',
                   'ybtQ', 'ybtP', 'ybtA', 'irp2', 'irp1', 'ybtU', 'ybtT', 'ybtE', 'fyuA', 'clbA',
                   'clbB', 'clbC', 'clbD', 'clbE', 'clbF', 'clbG', 'clbH', 'clbI', 'clbL', 'clbM',
                   'clbN', 'clbO', 'clbP', 'clbQ']
    full_header += mlst_header

    # If resistance genes are on, run the resBLAST.py script to get its headers.
    if args.resistance:
        res_out = subprocess.check_output('python ' + resblast +
                                          ' -s ' + data_folder + '/ARGannot_r2.fasta' +
                                          ' -t ' + data_folder + '/ARGannot_clustered80_r2.csv',
                                          shell=True)
        if res_out is bytes:
            res_out = res_out.decode()
        fields = res_out.split('\n')[0].rstrip().split('\t')
        res_headers = fields[1:]
        stdout_header += res_headers
        full_header += res_headers
    else:
        res_headers = []

    # # TO DO
    # if args.kaptive:
    #     pass

    return stdout_header, full_header, res_headers


def get_virulence_score(yb_group, cb_group, aerobactin, salmochelin, hypermucoidy):
    """
    The yersiniabactin locus counts as 1 and having any of the others counts as 2.
    This gives 4 possible scores:
      * 0 = no virulence
      * 1 = just yersiniabactin
      * 2 = one or more non-yersiniabactin
      * 3 = yersiniabactin and one or more non-yersiniabactin
    """
    score = 0
    if yb_group != '-':
        score += 1
    if cb_group != '-' or aerobactin != '-' or salmochelin != '-' or hypermucoidy != '-':
        score += 2
    return score


def get_resistance_score(res_headers, res_hits):
    """
    Three possible resistance scores:
      * 0 = no ESBL, no Carbepenemase
      * 1 = ESBL, no Carbepenemase
      * 2 = Carbepenemase (whether or not ESBL is present)
    """
    if not res_headers:
        return '-'

    # Look for a hit in any 'ESBL' column (e.g. 'Bla_ESBL' or 'Bla_ESBL_inhR')
    esbl_header_indices = [i for i, h in enumerate(res_headers) if 'esbl' in h.lower()]
    has_esbl = any(res_hits[i] != '-' for i in esbl_header_indices)

    # Look for a hit in any 'Carb' column (e.g. 'Bla_Carb')
    carb_header_indices = [i for i, h in enumerate(res_headers) if 'carb' in h.lower()]
    has_carb = any(res_hits[i] != '-' for i in carb_header_indices)

    if has_carb:
        return 2
    elif has_esbl:
        return 1
    else:
        return 0


def get_resistance_class_count(res_headers, res_hits):
    """
    Counts up all resistance gene classes, excluding the 'Bla' class which is intrinsic.
    """
    if not res_headers:
        return '-'
    res_indices = [i for i, h in enumerate(res_headers) if h.lower() != 'bla']
    return sum(0 if res_hits[i] == '-' else 1 for i in res_indices)


def get_resistance_gene_count(res_headers, res_hits):
    """
    Counts up all resistance genes, excluding the 'Bla' class which is intrinsic.
    """
    if not res_headers:
        return '-'
    res_indices = [i for i, h in enumerate(res_headers) if h.lower() != 'bla']
    return sum(0 if res_hits[i] == '-' else len(res_hits[i].split(';')) for i in res_indices)


def decompress_file(in_file, out_file):
    with gzip.GzipFile(in_file, 'rb') as i, open(out_file, 'wb') as o:
        s = i.read()
        o.write(s)


def get_klebsiella_species(contigs, data_folder):
    f = os.popen('mash dist ' + data_folder + '/species_mash_sketches.msh ' + contigs)

    best_species = None
    best_distance = 1.0

    for line in f:
        line_parts = line.split('\t')
        reference = line_parts[0]
        if len(line_parts) < 4:
            continue
        species = reference.split('/')[0]
        distance = float(line_parts[2])

        # Fix up the species name formatting a bit.
        species = species.replace('Escherichia_coli', 'Escherichia coli / Shigella')
        species = species.replace('_', ' ')
        species = species.replace(' subsp ', ' subsp. ')

        if distance < best_distance:
            best_distance = distance
            best_species = species

    if best_distance <= 0.01:
        return best_species, 'strong'
    elif best_distance <= 0.03:
        return best_species, 'weak'
    else:
        return 'unknown', ''


if __name__ == '__main__':
    main()
