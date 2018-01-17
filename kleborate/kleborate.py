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
    data_folder, mlstblast, resblast, clusterblast = get_resource_paths()
    kaptive_py, kaptive_k_db, kaptive_o_db = get_kaptive_paths()

    stdout_header, full_header, res_headers = get_output_headers(args, resblast, data_folder)
    output_headers(stdout_header, full_header, args.outfile)

    for contigs in args.assemblies:
        temp_decompressed_file, contigs = gunzip_contigs_if_necessary(contigs)

        # All results are stored in a dictionary where the key is the column name and the value
        # is the result. The results are outputted in order of the header rows. This means that the
        # column orders can be easily changed by modifying the get_output_headers function.

        results = {'strain': get_strain_name(contigs)}
        results.update(get_contig_stat_results(contigs))
        results.update(get_species_results(contigs, data_folder, args))
        results.update(get_chromosome_mlst_results(mlstblast, data_folder, contigs))
        results.update(get_ybt_mlst_results(mlstblast, data_folder, contigs))
        results.update(get_clb_mlst_results(mlstblast, data_folder, contigs))
        results.update(get_iuc_mlst_results(mlstblast, data_folder, contigs))
        results.update(get_iro_mlst_results(mlstblast, data_folder, contigs))
        results.update(get_hypermucoidy_results(clusterblast, data_folder, contigs))
        results.update(get_wzi_and_k_locus_results(mlstblast, data_folder, contigs))
        results.update(get_resistance_results(resblast, data_folder, contigs, args, res_headers))
        results.update(get_summary_results(results, res_headers))
        results.update(get_kaptive_results('K', kaptive_py, kaptive_k_db, contigs, args))
        results.update(get_kaptive_results('O', kaptive_py, kaptive_o_db, contigs, args))

        output_results(stdout_header, full_header, args.outfile, results)
        if temp_decompressed_file:
            os.remove(contigs)


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
    screening_args.add_argument('-s', '--species', action='store_true',
                                help='Turn on Klebsiella species identification (requires Mash, '
                                     'default: no species identification)')
    screening_args.add_argument('--kaptive_k', action='store_true',
                                help='Turn on Kaptive screening of K loci (default: do not run '
                                     'Kaptive for K loci)')
    screening_args.add_argument('--kaptive_o', action='store_true',
                                help='Turn on Kaptive screening of O loci (default: do not run '
                                     'Kaptive for O loci)')
    screening_args.add_argument('-k', '--kaptive', action='store_true',
                                help='Equivalent to --kaptive_k --kaptive_o')
    screening_args.add_argument('--all', action='store_true',
                                help='Equivalent to --resistance --species --kaptive')

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
        args.species = True
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
    if args.species:
        if not distutils.spawn.find_executable('mash'):
            sys.exit('Error: could not find mash')


def get_output_headers(args, resblast, data_folder):
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

    other_columns = ['Yersiniabactin', 'YbST',
                     'Colibactin', 'CbST',
                     'Aerobactin', 'AbST',
                     'Salmochelin', 'SmST',
                     'hypermucoidy', 'wzi', 'K_locus']
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

    # If resistance genes are on, run the resBLAST.py script to get its headers.
    if args.resistance:
        res_out = subprocess.check_output('python ' + resblast +
                                          ' -s ' + data_folder + '/ARGannot_r2.fasta' +
                                          ' -t ' + data_folder + '/ARGannot_clustered80_r2.csv',
                                          shell=True)
        if not isinstance(res_out, str):
            res_out = res_out.decode()
        fields = res_out.split('\n')[0].rstrip().split('\t')
        res_headers = fields[1:]
        stdout_header += res_headers
        full_header += res_headers
    else:
        res_headers = []

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
      * 0 = no ESBL, no carbapenemase
      * 1 = ESBL, no carbapenemase
      * 2 = Carbapenemase (whether or not ESBL is present)
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


def get_resource_paths():
    data_folder = resource_filename(__name__, 'data')
    mlstblast = resource_filename(__name__, 'mlstBLAST.py')
    resblast = resource_filename(__name__, 'resBLAST.py')
    clusterblast = resource_filename(__name__, 'clusterBLAST.py')
    return data_folder, mlstblast, resblast, clusterblast


def get_kaptive_paths():
    this_file = os.path.realpath(__file__)
    kaptive_dir = os.path.join(os.path.dirname(os.path.dirname(this_file)), 'kaptive')
    if not os.path.isdir(kaptive_dir):
        sys.exit('Error: could not find Kaptive directory. Did you git clone with --recursive?')
    kaptive_py = os.path.join(kaptive_dir, 'kaptive.py')
    if not os.path.isfile(kaptive_py):
        sys.exit('Error: could not find kaptive.py')
    db_dir = os.path.join(kaptive_dir, 'reference_database')
    kaptive_k_db = os.path.join(db_dir, 'Klebsiella_k_locus_primary_reference.gbk')
    if not os.path.isfile(kaptive_k_db):
        sys.exit('Error: could not find Klebsiella_k_locus_primary_reference.gbk')
    kaptive_o_db = os.path.join(db_dir, 'Klebsiella_o_locus_primary_reference.gbk')
    if not os.path.isfile(kaptive_o_db):
        sys.exit('Error: could not find Klebsiella_o_locus_primary_reference.gbk')
    return kaptive_py, kaptive_k_db, kaptive_o_db


def run_kaptive(kaptive_py, kaptive_db, contigs, output_file):
    kaptive_prefix = 'temp_kaptive_results_' + str(os.getpid())
    kaptive_table = kaptive_prefix + '_table.txt'

    p = subprocess.Popen('python ' + kaptive_py +
                         ' -a ' + contigs +
                         ' -k ' + kaptive_db +
                         ' -o ' + kaptive_prefix +
                         ' --verbose' +
                         ' --no_seq_out --no_json',
                         shell=True,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()

    # Make sure the output is a string, whether we are in Python 2 or 3.
    if not isinstance(stdout, str):
        stdout = stdout.decode()
    if not isinstance(stderr, str):
        stderr = stderr.decode()

    if p.returncode != 0:
        if stderr:
            sys.exit('Error: Kaptive failed to run with the following error:\n' + stderr.strip())
        else:
            sys.exit('Error: Kaptive failed to run')

    locus, confidence, problems, identity = None, None, None, None
    missing = []

    # Parse the required information from the Kaptive verbose output.
    output_lines = stdout.splitlines()
    missing_gene_lines = False
    for line in output_lines:
        if 'Best match locus:' in line:
            locus = line.split('Best match locus:')[1].strip()
        if 'Match confidence:' in line:
            confidence = line.split('Match confidence:')[1].strip()
        if 'Problems:' in line:
            problems = line.split('Problems:')[1].strip()
            if problems == 'None':
                problems = problems.lower()
        if 'Identity:' in line:
            identity = line.split('Identity:')[1].strip()
        if 'Other genes in locus:' in line:
            missing_gene_lines = False
        if missing_gene_lines:
            missing_gene = line.strip()
            if missing_gene:
                missing.append(missing_gene)
        if 'Missing expected genes:' in line:
            missing_gene_lines = True

    if not output_file:
        try:
            os.remove(kaptive_table)
        except OSError:
            pass

    else:  # if we are saving Kaptive results to file...
        with open(kaptive_table, 'rt') as f:
            kaptive_table_lines = f.readlines()
        assert len(kaptive_table_lines) == 2
        if not os.path.isfile(output_file):
            with open(output_file, 'wt') as f:
                f.write(kaptive_table_lines[0])  # write header line
        with open(output_file, 'at') as f:
            f.write(kaptive_table_lines[1])      # write data line

    if locus is None or confidence is None or problems is None or identity is None:
        sys.exit('Error: Kaptive failed to produce the expected output')

    return [locus, confidence, problems, identity, ','.join(missing)]


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


def gunzip_contigs_if_necessary(contigs):
    if get_compression_type(contigs) == 'gz':
        new_contigs = contigs + '_temp_decompress.fasta'
        decompress_file(contigs, new_contigs)
        contigs = new_contigs
        temp_decompress = True
    else:
        temp_decompress = False
    return temp_decompress, contigs


def get_contig_stat_results(contigs):
    contig_count, n50, longest_contig = get_contig_stats(contigs)
    return {'contig_count': str(contig_count),
            'N50': str(n50),
            'largest_contig': str(longest_contig)}


def get_species_results(contigs, data_folder, args):
    if args.species:
        species, species_hit_strength = get_klebsiella_species(contigs, data_folder)
    else:
        species, species_hit_strength = '', ''
    return {'species': species,
            'species_match': species_hit_strength}


def get_chromosome_mlst_results(mlstblast, data_folder, contigs):
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

    chromosome_mlst_header = get_chromosome_mlst_header()
    assert len(chromosome_mlst_header) == len(chr_st_detail)

    results = {'ST': chr_st,
               'Chr_ST': chr_st}
    results.update(dict(zip(get_chromosome_mlst_header(), chr_st_detail)))
    return results


def get_virulence_cluster_results(mlstblast, data_folder, contigs, alleles_fasta, profiles_txt,
                                  vir_name, vir_st_name, unknown_group_name, min_gene_count,
                                  header_function):
    f = os.popen('python ' + mlstblast +
                 ' -s ' + data_folder + '/' + alleles_fasta +
                 ' -d ' + data_folder + '/' + profiles_txt +
                 ' -i yes' +
                 ' --maxmissing 3' +
                 ' ' + contigs)
    st, group = '', ''
    st_detail = []
    for line in f:
        fields = line.rstrip().split('\t')
        if fields[2] != 'ST':  # skip header
            strain, st, group = fields[0], fields[2], fields[1]
            st_detail = fields[3:]

            # If no group was found but enough of the genes are present, then this strain is
            # labelled as an unknown lineage.
            if group == '':
                if sum(0 if x == '-' else 1 for x in st_detail) >= min_gene_count:
                    group = unknown_group_name
                    st = '0'
                else:
                    group = '-'
    f.close()

    mlst_header = header_function()
    assert len(mlst_header) == len(st_detail)

    results = {vir_name: group,
               vir_st_name: st}
    results.update(dict(zip(mlst_header, st_detail)))
    return results


def get_ybt_mlst_results(mlstblast, data_folder, contigs):
    return get_virulence_cluster_results(mlstblast, data_folder, contigs, 'ybt_alleles.fasta',
                                         'YbST_profiles.txt', 'Yersiniabactin', 'YbST',
                                         'ybt unknown', 8, get_ybt_mlst_header)


def get_clb_mlst_results(mlstblast, data_folder, contigs):
    return get_virulence_cluster_results(mlstblast, data_folder, contigs, 'clb_alleles.fasta',
                                         'CbST_profiles.txt', 'Colibactin', 'CbST',
                                         'clb unknown', 12, get_clb_mlst_header)


def get_iuc_mlst_results(mlstblast, data_folder, contigs):
    return get_virulence_cluster_results(mlstblast, data_folder, contigs, 'iuc_alleles.fasta',
                                         'AbST_profiles.txt', 'Aerobactin', 'AbST',
                                         'iuc unknown', 3, get_iuc_mlst_header)


def get_iro_mlst_results(mlstblast, data_folder, contigs):
    return get_virulence_cluster_results(mlstblast, data_folder, contigs, 'iro_alleles.fasta',
                                         'SmST_profiles.txt', 'Salmochelin', 'SmST',
                                         'iro unknown', 3, get_iro_mlst_header)


def get_hypermucoidy_results(clusterblast, data_folder, contigs):
    hypermucoidy = '-'
    f = os.popen('python ' + clusterblast +
                 ' -s ' + data_folder + '/hypermucoidy.fasta' +
                 ' ' + contigs)
    for line in f:
        fields = line.rstrip().split('\t')
        if fields[1] != 'hypermucoidy':  # skip header
            hypermucoidy = fields[1]
    f.close()

    return {'hypermucoidy': hypermucoidy}


def get_wzi_and_k_locus_results(mlstblast, data_folder, contigs):
    wzi_st, k_type = '-', '-'
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
    f.close()

    return {'wzi': wzi_st,
            'K_locus': k_type}


def get_resistance_results(resblast, data_folder, contigs, args, res_headers):
    res_hits = []
    if args.resistance:
        f = os.popen('python ' + resblast +
                     ' -s ' + data_folder + '/ARGannot_r2.fasta' +
                     ' -t ' + data_folder + '/ARGannot_clustered80_r2.csv' +
                     ' -q ' + data_folder + '/QRDR_120.aa' +
                     ' -r ' + data_folder + '/MgrB_and_PmrB.aa' +
                     ' ' + contigs)
        for line in f:
            fields = line.rstrip().split('\t')
            if fields[0] != 'strain':  # skip header
                res_hits = fields[1:]
        f.close()
    else:
        res_hits = [''] * len(res_headers)

    assert len(res_headers) == len(res_hits)
    return dict(zip(res_headers, res_hits))


def get_kaptive_results(locus_type, kaptive_py, kaptive_db, contigs, args):
    assert locus_type == 'K' or locus_type == 'O'

    headers = ['K_locus', 'K_locus_confidence', 'K_locus_problems', 'K_locus_identity',
               'K_locus_missing_genes']
    if locus_type == 'O':
        headers = [x.replace('K_locus', 'O_locus') for x in headers]

    kaptive_results = ['', '', '', '', '']
    if (args.kaptive_k and locus_type == 'K') or (args.kaptive_o and locus_type == 'O'):
        kaptive_results = run_kaptive(kaptive_py, kaptive_db, contigs, args.kaptive_k_outfile)
        assert len(headers) == len(kaptive_results)
        return dict(zip(headers, kaptive_results))
    else:
        return {}


def get_summary_results(results, res_headers):
    res_hits = [results[x] for x in res_headers]
    return {'virulence_score': str(get_virulence_score(results['Yersiniabactin'],
                                                       results['Colibactin'],
                                                       results['Aerobactin'],
                                                       results['Salmochelin'],
                                                       results['hypermucoidy'])),
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
