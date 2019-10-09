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

import argparse
import distutils.spawn
import gzip
import os
import subprocess
import sys
import tempfile

from pkg_resources import resource_filename
from .contig_stats import get_contig_stat_results
from .kaptive import get_kaptive_paths, get_kaptive_results
from .species import get_species_results, is_kp_complex
from .mlstBLAST import mlst_blast
from .resBLAST import read_class_file, get_res_headers, resblast_one_assembly
from .rmpA import rmpa_blast
from .misc import get_compression_type, load_fasta
from .version import __version__


def main():
    args = parse_arguments()
    check_inputs_and_programs(args)
    data_folder = get_data_path()
    kaptive_py, kaptive_k_db, kaptive_o_db = get_kaptive_paths()

    stdout_header, full_header, res_headers = get_output_headers(args, data_folder)
    output_headers(stdout_header, full_header, args.outfile)

    with tempfile.TemporaryDirectory() as tmp_dir:
        for contigs in args.assemblies:
            contigs = gunzip_contigs_if_necessary(contigs, tmp_dir)

            # All results are stored in a dictionary where the key is the column name and the value
            # is the result. The results are outputted in order of the header rows. This means that
            # the column orders can be easily changed by modifying the get_output_headers function.

            results = {'strain': get_strain_name(contigs)}
            results.update(get_contig_stat_results(contigs))
            results.update(get_species_results(contigs, data_folder))

            kp_complex = is_kp_complex(results)

            results.update(get_chromosome_mlst_results(data_folder, contigs, kp_complex))
            results.update(get_ybt_mlst_results(data_folder, contigs))
            results.update(get_clb_mlst_results(data_folder, contigs))
            results.update(get_iuc_mlst_results(data_folder, contigs))
            results.update(get_iro_mlst_results(data_folder, contigs))
            results.update(get_hypermucoidy_results(data_folder, contigs))
            results.update(get_wzi_and_k_locus_results(data_folder, contigs))
            results.update(get_resistance_results(data_folder, contigs, args, res_headers, kp_complex))
            results.update(get_summary_results(results, res_headers))
            results.update(get_kaptive_results('K', kaptive_py, kaptive_k_db, contigs, args))
            results.update(get_kaptive_results('O', kaptive_py, kaptive_o_db, contigs, args))

            output_results(stdout_header, full_header, args.outfile, results)


def parse_arguments():
    parser = argparse.ArgumentParser(description='Kleborate: a tool for characterising '
                                                 'virulence and resistance in Klebsiella',
                                     add_help=False)

    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument('-a', '--assemblies', nargs='+', type=str, required=True,
                               help='FASTA file(s) for assemblies')

    screening_args = parser.add_argument_group('Screening options')
    screening_args.add_argument('-r', '--resistance', action='store_true',
                                help='Turn on resistance genes screening (default: no resistance '
                                     'gene screening)')
    screening_args.add_argument('--kaptive_k', action='store_true',
                                help='Turn on Kaptive screening of K loci (default: do not run '
                                     'Kaptive for K loci)')
    screening_args.add_argument('--kaptive_o', action='store_true',
                                help='Turn on Kaptive screening of O loci (default: do not run '
                                     'Kaptive for O loci)')
    screening_args.add_argument('-k', '--kaptive', action='store_true',
                                help='Equivalent to --kaptive_k --kaptive_o')
    screening_args.add_argument('--all', action='store_true',
                                help='Equivalent to --resistance --kaptive')

    output_args = parser.add_argument_group('Output options')
    output_args.add_argument('-o', '--outfile', type=str, default='Kleborate_results.txt',
                             help='File for detailed output (default: Kleborate_results.txt)')
    output_args.add_argument('--kaptive_k_outfile', type=str,
                             help='File for full Kaptive K locus output (default: do not '
                                  'save Kaptive K locus results to separate file)')
    output_args.add_argument('--kaptive_o_outfile', type=str,
                             help='File for full Kaptive O locus output (default: do not '
                                  'save Kaptive O locus results to separate file)')

    help_args = parser.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                           help='Show this help message and exit')
    help_args.add_argument('--version', action='version', version=__version__,
                           help="Show program's version number and exit")

    # If no arguments were used, print the entire help (argparse default is to just give an error
    # like '-a is required').
    if len(sys.argv) == 1:
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    if args.kaptive:
        args.kaptive_k = True
        args.kaptive_o = True
    if args.all:
        args.resistance = True
        args.kaptive_k = True
        args.kaptive_o = True

    if args.kaptive_k_outfile and not args.kaptive_k:
        sys.exit('Error: you must use --kaptive_k (or --kaptive) to use --kaptive_k_outfile')
    if args.kaptive_o_outfile and not args.kaptive_o:
        sys.exit('Error: you must use --kaptive_o (or --kaptive) to use --kaptive_o_outfile')

    return args


def check_inputs_and_programs(args):
    for assembly in args.assemblies:
        if os.path.isdir(assembly):
            sys.exit('Error: ' + assembly + ' is a directory (please specify assembly files)')
        if not os.path.isfile(assembly):
            sys.exit('Error: could not find ' + assembly)
        fasta = load_fasta(assembly)
        if len(fasta) < 1:
            sys.exit('Error: invalid FASTA file: ' + assembly)
        for record in fasta:
            header, seq = record
            if len(seq) == 0:
                sys.exit('Error: invalid FASTA file (contains a zero-length sequence): ' + assembly)
    if not distutils.spawn.find_executable('makeblastdb'):
        sys.exit('Error: could not find makeblastdb')
    if not distutils.spawn.find_executable('blastn'):
        sys.exit('Error: could not find blastn')
    if args.resistance:
        if not distutils.spawn.find_executable('blastx'):
            sys.exit('Error: could not find blastx')
    if not distutils.spawn.find_executable('mash'):
        sys.exit('Error: could not find mash')
    major, minor, patch = get_blast_version()
    if major < 2 or (major == 2 and minor < 7):
        sys.exit('Error: you have BLAST v{}.{}.{} installed, but Kleborate requires v2.7.1 or '
                 'later'.format(major, minor, patch))


def get_blast_version():
    command = ['blastn', '-version']
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, _ = process.communicate()
    out = out.decode()
    try:
        version = out.split(': ')[1].split()[0].split('+')[0]
        major, minor, patch = version.split('.')
        return int(major), int(minor), patch
    except (IndexError, ValueError):
        sys.exit('Error: could not determine BLAST version')


def get_output_headers(args, data_folder):
    """
    There are two levels of output:
      * stdout is simpler and displayed to the console
      * full contains more and is saved to file
    This function returns headers for both. It also returns the resistance headers in a separate
    list, as they are used to total up some resistance summaries.
    """
    stdout_header = ['strain', 'species']
    full_header = ['strain', 'species', 'species_match']
    stdout_header += ['ST', 'virulence_score']
    full_header += ['contig_count', 'N50', 'largest_contig', 'ambiguous_bases', 'ST',
                    'virulence_score']

    if args.resistance:
        stdout_header.append('resistance_score')
        full_header.append('resistance_score')
        full_header.append('num_resistance_classes')
        full_header.append('num_resistance_genes')

    other_columns = ['Yersiniabactin', 'YbST',
                     'Colibactin', 'CbST',
                     'Aerobactin', 'AbST',
                     'Salmochelin', 'SmST',
                     'rmpA', 'rmpA2',
                     'wzi', 'K_locus']
    stdout_header += other_columns
    full_header += other_columns

    if args.kaptive_k:
        # K_locus is already in the header.
        stdout_header.append('K_locus_confidence')
        full_header.append('K_locus_problems')
        full_header.append('K_locus_confidence')
        full_header.append('K_locus_identity')
        full_header.append('K_locus_missing_genes')

    if args.kaptive_o:
        stdout_header.append('O_locus')
        stdout_header.append('O_locus_confidence')
        full_header.append('O_locus')
        full_header.append('O_locus_problems')
        full_header.append('O_locus_confidence')
        full_header.append('O_locus_identity')
        full_header.append('O_locus_missing_genes')

    full_header.append('Chr_ST')
    full_header += get_chromosome_mlst_header()
    full_header += get_ybt_mlst_header()
    full_header += get_clb_mlst_header()
    full_header += get_iuc_mlst_header()
    full_header += get_iro_mlst_header()

    # If resistance genes are on, run the resBLAST.py script to get its headers.
    if args.resistance:
        gene_info, res_classes, bla_classes = \
            read_class_file(data_folder + '/ARGannot_clustered80_r3.csv')
        res_headers = get_res_headers(res_classes, bla_classes)
        stdout_header += res_headers
        full_header += res_headers
    else:
        res_headers = []

    return stdout_header, full_header, res_headers


def get_virulence_score(yersiniabactin, colibactin, aerobactin):
    """
    Six possible virulence scores:
      * 0 = no virulence (no yersiniabactin, colibactin or aerobactin)
      * 1 = just yersiniabactin (no colibactin or aerobactin)
      * 2 = colibactin but no aerobactin (regardless of yersiniabactin, which is probably present)
      * 3 = just aerobactin (no yersiniabactin or colibactin)
      * 4 = aerobactin and yersiniabactin (but not colibactin)
      * 5 = colibactin and aerobactin (regardless of yersiniabactin, which is probably present)
    """
    has_ybt = (yersiniabactin != '-')
    has_aero = (aerobactin != '-')
    has_coli = (colibactin != '-')

    if has_coli and has_aero:
        return 5
    elif has_aero and has_ybt:
        return 4
    elif has_aero:
        return 3
    elif has_coli:
        return 2
    elif has_ybt:
        return 1
    else:
        return 0


def get_resistance_score(res_headers, res_hits):
    """
    Four possible resistance scores:
      * 0 = no ESBL, no carbapenemase (regardless of colistin resistance)
      * 1 = ESBL, no carbapenemase (regardless of colistin resistance)
      * 2 = Carbapenemase without colistin resistance
      * 3 = Carbapenemase and colistin resistance
    """
    if not res_headers:
        return '-'

    # Look for a hit in any 'ESBL' column (e.g. 'Bla_ESBL' or 'Bla_ESBL_inhR').
    esbl_header_indices = [i for i, h in enumerate(res_headers) if 'esbl' in h.lower()]
    has_esbl = any(res_hits[i] != '-' for i in esbl_header_indices)

    # Look for a hit in any 'Carb' column (e.g. 'Bla_Carb').
    carb_header_indices = [i for i, h in enumerate(res_headers) if 'carb' in h.lower()]
    has_carb = any(res_hits[i] != '-' for i in carb_header_indices)

    # Look for a hit in the 'Col' column.
    col_header_indices = [i for i, h in enumerate(res_headers) if h.lower() == 'col']
    has_col = any(res_hits[i] != '-' for i in col_header_indices)

    if has_carb and has_col:
        return 3
    elif has_carb:
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
    res_indices = [i for i, h in enumerate(res_headers)
                   if h.lower() != 'bla' and h.lower() != 'omp']
    return sum(0 if res_hits[i] == '-' else 1 for i in res_indices)


def get_resistance_gene_count(res_headers, res_hits):
    """
    Counts up all resistance genes, excluding the 'Bla' class which is intrinsic.
    """
    if not res_headers:
        return '-'
    res_indices = [i for i, h in enumerate(res_headers)
                   if h.lower() != 'omp']

    gene_list = []
    for i in res_indices:
        genes = res_hits[i].split(';')

        # Exclude mutation-based flq resistance.
        genes = [g for g in genes if 'gyra-' not in g.lower()]
        genes = [g for g in genes if 'parc-' not in g.lower()]

        # Exclude truncation-based resistance.
        genes = [g for g in genes if '%' not in g]

        # Exclude intrinsic bla genes.
        genes = [g for g in genes if 'shv-' not in g.lower()]
        genes = [g for g in genes if 'okp-' not in g.lower()]
        genes = [g for g in genes if 'len-' not in g.lower()]

        genes = [g for g in genes if g != '-']
        gene_list += genes

    return len(gene_list)


def get_data_path():
    return resource_filename(__name__, 'data')


def get_chromosome_mlst_header():
    return ['gapA', 'infB', 'mdh', 'pgi', 'phoE', 'rpoB', 'tonB']


def get_ybt_mlst_header():
    return ['ybtS', 'ybtX', 'ybtQ', 'ybtP', 'ybtA', 'irp2', 'irp1', 'ybtU', 'ybtT', 'ybtE', 'fyuA']


def get_clb_mlst_header():
    return ['clbA', 'clbB', 'clbC', 'clbD', 'clbE', 'clbF', 'clbG', 'clbH', 'clbI', 'clbL', 'clbM',
            'clbN', 'clbO', 'clbP', 'clbQ']


def get_iuc_mlst_header():
    return ['iucA', 'iucB', 'iucC', 'iucD', 'iutA']


def get_iro_mlst_header():
    return ['iroB', 'iroC', 'iroD', 'iroN']


def gunzip_contigs_if_necessary(contigs, temp_dir):
    if get_compression_type(contigs) == 'gz':
        name = get_strain_name(contigs)
        new_contigs = temp_dir + '/' + name + '.fasta'
        decompress_file(contigs, new_contigs)
        return new_contigs
    else:
        return contigs


def decompress_file(in_file, out_file):
    with gzip.GzipFile(in_file, 'rb') as i, open(out_file, 'wb') as o:
        s = i.read()
        o.write(s)


def get_chromosome_mlst_results(data_folder, contigs, kp_complex):
    chromosome_mlst_header = get_chromosome_mlst_header()

    if kp_complex:
        seqs = data_folder + '/Klebsiella_pneumoniae.fasta'
        database = data_folder + '/kpneumoniae.txt'
        results = mlst_blast(seqs, database, 'no', [contigs], minident=95, maxmissing=3,
                             print_header=False)
        chr_st, chr_st_detail = results[1], results[2:]
        if chr_st != '0':
            chr_st = 'ST' + chr_st

        # ST67 and ST90 get special 'subspecies' names.
        chr_st_with_subsp = chr_st
        if chr_st_with_subsp == 'ST90':
            chr_st_with_subsp = 'ST90 (subsp. ozanae)'
        if chr_st_with_subsp == 'ST67':
            chr_st_with_subsp = 'ST67 (subsp. rhinoscleromatis)'
            
        assert len(chromosome_mlst_header) == len(chr_st_detail)

        results = {'ST': chr_st_with_subsp,
                   'Chr_ST': chr_st}

    else:
        results = {'ST': "NA",
                   'Chr_ST': "NA"}
        chr_st_detail = ['-'] * len(chromosome_mlst_header)

    results.update(dict(zip(get_chromosome_mlst_header(), chr_st_detail)))
    return results


def get_virulence_cluster_results(data_folder, contigs, alleles_fasta, profiles_txt,
                                  vir_name, vir_st_name, unknown_group_name, min_gene_count,
                                  header_function):
    seqs = data_folder + '/' + alleles_fasta
    database = data_folder + '/' + profiles_txt
    results = mlst_blast(seqs, database, 'yes', [contigs], minident=95, maxmissing=3,
                         print_header=False)
    group, st, st_detail = results[1], results[2], results[3:]
    if group == '':
        if sum(0 if x == '-' else 1 for x in st_detail) >= min_gene_count:
            group = unknown_group_name
            st = '0'
        else:
            group = '-'

    mlst_header = header_function()
    assert len(mlst_header) == len(st_detail)

    results = {vir_name: group,
               vir_st_name: st}
    results.update(dict(zip(mlst_header, st_detail)))
    return results


def get_ybt_mlst_results(data_folder, contigs):
    return get_virulence_cluster_results(data_folder, contigs, 'ybt_alleles.fasta',
                                         'YbST_profiles.txt', 'Yersiniabactin', 'YbST',
                                         'ybt unknown', 8, get_ybt_mlst_header)


def get_clb_mlst_results(data_folder, contigs):
    return get_virulence_cluster_results(data_folder, contigs, 'clb_alleles.fasta',
                                         'CbST_profiles.txt', 'Colibactin', 'CbST',
                                         'clb unknown', 12, get_clb_mlst_header)


def get_iuc_mlst_results(data_folder, contigs):
    return get_virulence_cluster_results(data_folder, contigs, 'iuc_alleles.fasta',
                                         'AbST_profiles.txt', 'Aerobactin', 'AbST',
                                         'iuc unknown', 3, get_iuc_mlst_header)


def get_iro_mlst_results(data_folder, contigs):
    return get_virulence_cluster_results(data_folder, contigs, 'iro_alleles.fasta',
                                         'SmST_profiles.txt', 'Salmochelin', 'SmST',
                                         'iro unknown', 3, get_iro_mlst_header)


def get_hypermucoidy_results(data_folder, contigs):
    seqs = data_folder + '/hypermucoidy.fasta'
    database = data_folder + '/hypermucoidy_rmpA.txt'
    results = rmpa_blast(seqs, database, [contigs], 95.0)
    rmpa_allele, rmpa2_allele = results[1], results[2]

    return {'rmpA': rmpa_allele,
            'rmpA2': rmpa2_allele}


def get_wzi_and_k_locus_results(data_folder, contigs):
    seqs = data_folder + '/wzi.fasta'
    database = data_folder + '/wzi.txt'
    results = mlst_blast(seqs, database, 'yes', [contigs], minident=95, maxmissing=0,
                         print_header=False)
    k_type = results[1]
    if results[2] == '0':
        wzi_st = '-'
    else:
        wzi_st = 'wzi' + results[2]
    if k_type == '':
        k_type = '-'

    return {'wzi': wzi_st,
            'K_locus': k_type}


def get_resistance_results(data_folder, contigs, args, res_headers, kp_complex):
    if args.resistance:
        gene_info, _, _ = read_class_file(data_folder + '/ARGannot_clustered80_r3.csv')

        # Only do mutation/truncation tests for Kp complex species.
        if kp_complex:
            qrdr = data_folder + '/QRDR_120.fasta'
            trunc = data_folder + '/MgrB_and_PmrB.fasta'
            omp = data_folder + '/OmpK.fasta'
        else:
            qrdr, trunc, omp = None, None, None

        seqs = data_folder + '/ARGannot_r3.fasta'
        res_hits = resblast_one_assembly(contigs, gene_info, qrdr, trunc, omp, seqs, 80.0, 90.0)
        return {r: ';'.join(sorted(res_hits[r])) if r in res_hits else '-'
                for r in res_headers}
    else:
        return {}


def get_summary_results(results, res_headers):
    res_hits = [results[x] for x in res_headers]
    return {'virulence_score': str(get_virulence_score(results['Yersiniabactin'],
                                                       results['Colibactin'],
                                                       results['Aerobactin'])),
            'resistance_score': str(get_resistance_score(res_headers, res_hits)),
            'num_resistance_classes': str(get_resistance_class_count(res_headers, res_hits)),
            'num_resistance_genes': str(get_resistance_gene_count(res_headers, res_hits))}


def output_headers(stdout_header, full_header, outfile):
    print('\t'.join(stdout_header))
    with open(outfile, 'wt') as o:
        o.write('\t'.join(full_header))
        o.write('\n')


def output_results(stdout_header, full_header, outfile, results):
    print('\t'.join([results[x] for x in stdout_header]))
    with open(outfile, 'at') as o:
        o.write('\t'.join([results[x] for x in full_header]))
        o.write('\n')


def get_strain_name(full_path):
    filename = os.path.split(full_path)[1]
    if filename.endswith('_temp_decompress.fasta'):
        filename = filename[:-22]
    if filename.endswith('.gz'):
        filename = filename[:-3]
    return os.path.splitext(filename)[0]


if __name__ == '__main__':
    main()
