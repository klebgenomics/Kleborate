"""
Copyright 2024 Mary Maranga, Kat Holt, Ryan Wick
https://github.com/klebgenomics/KleborateModular/

This file is part of Kleborate. Kleborate is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Kleborate is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Kleborate. If
not, see <https://www.gnu.org/licenses/>.
"""

import argparse
import graphlib
import gzip
import importlib
import importlib.metadata
import os
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile
import textwrap
import uuid
from glob import glob

from .shared.help_formatter import MyParser, MyHelpFormatter
from .shared.misc import get_compression_type, load_fasta,reverse_complement
from .shared.species_defs import is_kp_complex, is_ko_complex, is_escherichia


def parse_arguments(args, all_module_names, modules):
    """
    This function does the CLI argument parsing for Kleborate. Module-specific arguments are added
    by each module's add_cli_options function.
    """
    parser = MyParser(description='Kleborate: a tool for characterising virulence and resistance '
                                  'in pathogen assemblies',
                      formatter_class=MyHelpFormatter, add_help=False, epilog=paper_refs())

    if '--helpall' in args or '--allhelp' in args or '--all_help' in args:
        args.append('--help_all')

    io_args = parser.add_argument_group('Input/output')
    io_args.add_argument('-a', '--assemblies', nargs='+', type=str,
                         help='FASTA file(s) for assemblies')

    io_args.add_argument('-o', '--outdir', type=str,
                         help='Directory for storing output files')

    io_args.add_argument('-r', '--resume', action='store_true',
                         help='append the output files')

    io_args.add_argument('--trim_headers', action='store_true',
                         help='Trim headers in the output files')

    module_args = parser.add_argument_group('Modules')
    module_args.add_argument('--list_modules', action='store_true',
                             help='Print a list of all available modules and then quit')
    module_args.add_argument('-p', '--preset', type=str,
                             help=f'Module presets, choose from: ' + ', '.join(get_presets()))
    module_args.add_argument('-m', '--modules', type=str,
                             help='Comma-delimited list of Kleborate modules to use')

    add_module_cli_arguments(parser, args, all_module_names, modules)

    help_args = parser.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                           help='Show this help message and exit')
    help_args.add_argument('--help_all', action='help',
                           help='Show a help message with all module options')
    help_args.add_argument('--version', action='version', version=f'Kleborate v{get_version()}',
                           help="Show program's version number and exit")

<<<<<<< HEAD
    if not args:
        parser.print_help(file=sys.stderr)
        sys.exit(1)
=======
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
    full_header += ['contig_count', 'N50', 'largest_contig', 'total_size', 'ambiguous_bases',
                    'QC_warnings', 'ST', 'virulence_score']

    if args.resistance:
        stdout_header.append('resistance_score')
        full_header.append('resistance_score')
        full_header.append('num_resistance_classes')
        full_header.append('num_resistance_genes')

    other_columns = ['Yersiniabactin', 'YbST',
                     'Colibactin', 'CbST',
                     'Aerobactin', 'AbST',
                     'Salmochelin', 'SmST',
                     'RmpADC', 'RmST',
                     'rmpA2',
                     'wzi', 'K_locus']
    stdout_header += other_columns
    full_header += other_columns

    if args.kaptive_k:
        stdout_header.append('K_locus_confidence')
        full_header.append('K_type')
        full_header.append('K_locus_problems')
        full_header.append('K_locus_confidence')
        full_header.append('K_locus_identity')
        full_header.append('K_locus_missing_genes')

    if args.kaptive_o:
        stdout_header.append('O_locus')
        stdout_header.append('O_locus_confidence')
        full_header.append('O_locus')
        full_header.append('O_type')
        full_header.append('O_locus_problems')
        full_header.append('O_locus_confidence')
        full_header.append('O_locus_identity')
        full_header.append('O_locus_missing_genes')

    if args.resistance:
        gene_info, res_classes, bla_classes = \
            read_class_file(data_folder + '/CARD_AMR_clustered.csv')
        res_headers = get_res_headers(res_classes, bla_classes)
        res_headers += ['truncated_resistance_hits', 'spurious_resistance_hits']
        res_headers += ['Cipro_prediction','Cipro_prediction_support']
        stdout_header += res_headers
        full_header += res_headers
    else:
        res_headers = []
>>>>>>> 212fb361b4dc2dc6eda10c33582d38d9a6b132f5

    return parser.parse_args(args)


def main():
    all_module_names, modules = import_modules()
    args = parse_arguments(sys.argv[1:], all_module_names, modules)
    print_modules(args, all_module_names, modules)

    module_names, check_module_list, pass_modules = get_used_module_names(args, all_module_names, get_presets())

    # Define preset_check_modules
    preset_check_modules = []
    if args.preset:
        presets = get_presets()
        preset_check_modules = [module for module, _ in presets[args.preset]['check']]

    module_names, module_run_order, external_programs = check_modules(args, modules, module_names, check_module_list, pass_modules)

    full_headers, stdout_headers = get_headers(module_names, modules)
    print('\t'.join([h.split('__')[-1] for h in stdout_headers]))

    # Ensure the output directory exists
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    # If the resume flag is not set, remove existing output files
    if args.modules:
        module_name = args.modules.split(',')[0]  
        out_files_suffixes = [f'{module_name}_output.txt']
    else:
        out_files_suffixes = ['klebsiella_pneumo_complex_output.txt',
                              'klebsiella_oxytoca_complex_output.txt',
                              'escherichia_output.txt']
    if not args.resume:
        for suffix in out_files_suffixes:
            for file in glob(f'{args.outdir}/*{suffix}'):
                os.remove(file)


    for assembly in args.assemblies:
        check_assembly(assembly)  # Check assembly before processing

        with tempfile.TemporaryDirectory() as temp_dir:
            unzipped_assembly = gunzip_assembly_if_necessary(assembly, temp_dir)
            minimap2_index = build_minimap2_index(assembly, unzipped_assembly, external_programs, temp_dir)
            results = {'strain': get_strain_name(assembly)}

            pass_check = True  # default, assume no check and run all modules

            # if we have 'check' modules in the preset, run these
            if args.preset and len(check_module_list) > 0:
                for module, check in presets[args.preset]['check']:
                    try:
                        module_results = modules[module].get_results(unzipped_assembly, minimap2_index, args, results)

                        results.update({f'{module}__{header}': result for header, result in module_results.items()})
                        check_function = globals()[check]

                        if not check_function(module_results):
                            pass_check = False
                            print(f"Assembly {assembly} failed in check {check}.")
                            break  # Exit the for loop since this assembly failed the check

                    except Exception as e:
                        print(f"Error encountered while processing {assembly} with {module}: {e}.")
                        pass_check = False
                        break  # Exit the for loop since an error occurred

            # proceed through all other modules
            if pass_check:
                for module in module_run_order:
                    if module not in preset_check_modules:
                        module_results = modules[module].get_results(unzipped_assembly, minimap2_index, args, results)
                        results.update({f'{module}__{header}': result for header, result in module_results.items()})
            else:
                # Populate results with "Not Tested" for modules that did not run
                for module in module_run_order:
                    if module not in preset_check_modules:
                        module_headers = [header for header in full_headers if header.startswith(module)]
                        for header in module_headers:
                            results[header] = 'Not Tested'

            # Split the results based on species
            if args.modules:
                module_name = args.modules.split(',')[0] 
                outfile_suffix = f'{module_name}_output.txt'
            else:
                # Determine the appropriate output file suffix based on species
                species = results.get('enterobacterales__species__species', None)
                if species and is_kp_complex({'species': species}):
                    outfile_suffix = 'klebsiella_pneumo_complex_output.txt'
                elif species and is_ko_complex({'species': species}):
                    outfile_suffix = 'klebsiella_oxytoca_complex_output.txt'
                elif species and is_escherichia({'species': species}):
                    outfile_suffix = 'escherichia_output.txt'
                else:
                    print(f"Assembly {assembly} does not match any specified species. Skipping to next assembly.")
                    continue

            # write results
            output_file = os.path.join(args.outdir, outfile_suffix)
            output_results(full_headers, stdout_headers, output_file, results, args.trim_headers)


# def main(): 
#     all_module_names, modules = import_modules()
#     args = parse_arguments(sys.argv[1:], all_module_names, modules)
#     print_modules(args, all_module_names, modules)
    
#     module_names, check_module_list, pass_modules = get_used_module_names(args, all_module_names, get_presets())  
    
#     # Define preset_check_modules 
#     preset_check_modules = []
#     if args.preset:
#         presets = get_presets()
#         preset_check_modules = [module for module, _ in presets[args.preset]['check']]

#     module_names, module_run_order, external_programs = check_modules(args, modules, module_names, check_module_list, pass_modules)
#     check_assemblies(args)

#     full_headers, stdout_headers = get_headers(module_names, modules)  
#     output_headers(full_headers, stdout_headers, args.outfile) 

#     for assembly in args.assemblies:
#         with tempfile.TemporaryDirectory() as temp_dir:
#             unzipped_assembly = gunzip_assembly_if_necessary(assembly, temp_dir)
#             minimap2_index = build_minimap2_index(assembly, unzipped_assembly, external_programs, temp_dir)
#             results = {'strain': get_strain_name(assembly)}

#             pass_check = True  # default, assume no check and run all modules

#             # if we have 'check' modules in the preset, run these
#             if args.preset and len(check_module_list) > 0:
#                 for module, check in presets[args.preset]['check']:
#                     try:
#                         module_results = modules[module].get_results(unzipped_assembly, minimap2_index, args, results)

#                         results.update({f'{module}__{header}': result for header, result in module_results.items()})
#                         check_function = globals()[check]

#                         if not check_function(module_results):
#                             pass_check = False
#                             print(f"Assembly {assembly} failed in check {check}. Continuing with next assembly.")
#                             break  # Exit the for loop since this assembly failed the check

#                     except Exception as e:
#                         print(f"Error encountered while processing {assembly} with {module}: {e}. Continuing with next assembly.")
#                         pass_check = False
#                         break  # Exit the for loop since an error occurred

#             # proceed through all other modules
#             if pass_check:
#                 for module in module_run_order:
#                     if module not in preset_check_modules:
#                         module_results = modules[module].get_results(unzipped_assembly, minimap2_index, args, results)

#                         results.update({f'{module}__{header}': result for header, result in module_results.items()})

#             # write results
#             output_results(full_headers, stdout_headers, args.outfile, results)



def print_modules(args, all_module_names, modules):
    if args.list_modules:
        print()
        print('Available modules for Kleborate')
        print('-------------------------------')
        terminal_width = shutil.get_terminal_size().columns
        end_formatting, bold = '\033[0m', '\033[1m'
        for m in all_module_names:
            description = modules[m].description()
            text = f'{bold}{m}{end_formatting}: {description}'
            print('\n'.join(textwrap.wrap(text, width=terminal_width - 1)))
            print()
        sys.exit(0)
    elif not args.assemblies:
        sys.exit('Error: you must provide one or more assembly files using --assemblies')



def get_presets():
    kpsc_modules = {
        'check': [('enterobacterales__species', 'is_kp_complex')],
        'pass': [
            'general__contig_stats','klebsiella_pneumo_complex__mlst',
            'klebsiella__ybst', 'klebsiella__cbst', 'klebsiella__abst', 'klebsiella__smst', 'klebsiella__rmst', 'klebsiella_pneumo_complex__virulence_score',
            'klebsiella__rmpa2','klebsiella_pneumo_complex__amr', 'klebsiella_pneumo_complex__resistance_score', 'klebsiella_pneumo_complex__resistance_class_count',
            'klebsiella_pneumo_complex__resistance_gene_count', 'klebsiella_pneumo_complex__wzi','klebsiella_pneumo_complex__kaptive', 'klebsiella_pneumo_complex__cipro_prediction'
        ]
    }

    kosc_modules = {
        'check': [('enterobacterales__species', 'is_ko_complex')],
        'pass': [
            'general__contig_stats',
            'klebsiella_oxytoca_complex__mlst', 'klebsiella__ybst', 'klebsiella__cbst', 'klebsiella__abst', 'klebsiella__smst', 'klebsiella__rmst','klebsiella__rmpa2'
        ]
    }

    escherichia_modules = {
        'check': [('enterobacterales__species', 'is_escherichia')],
        'pass': [
            'general__contig_stats',
            'escherichia__mlst_achtman', 'escherichia__mlst_pasteur'
        ]
    }

    return {
        'kpsc': kpsc_modules,
        'kosc': kosc_modules,
        'escherichia': escherichia_modules
    }


def add_module_cli_arguments(parser, args, all_module_names, modules):
    """
    This function add CLI argument for modules. Each modules that has options gets its own argument
    group. These are only displayed in the help text if the user used --help_all, otherwise they
    are hidden.
    """
    for m in all_module_names:
        group = modules[m].add_cli_options(parser)
        if '--help_all' not in args and group is not None:
            for a in group._group_actions:
                a.help = argparse.SUPPRESS


def get_used_module_names(args, all_module_names, presets): 
    if args.preset is None and args.modules is None:
        sys.exit('Error: either --preset or --modules is required')

    # Initialize empty lists for module names, check modules, and pass modules
    module_names = []
    check_modules = []
    pass_modules = []

    if args.preset:
        if args.preset not in presets:
            sys.exit(f'Error: {args.preset} is not a valid preset')

        # Assuming presets[args.preset] is a dictionary with 'check' and 'pass' keys
        check_modules = [module[0] for module in presets[args.preset].get('check', [])]  # Extract module names from check modules
        pass_modules = presets[args.preset].get('pass', [])  # Directly assign pass modules

        module_names += check_modules + pass_modules  # Combine check and pass modules for the overall list

    if args.modules:
        for m in args.modules.split(','):
            if m not in all_module_names:
                sys.exit(f'Error: {m} is not a valid module name')
            if m not in module_names:
                module_names.append(m)

    return module_names, check_modules, pass_modules


def get_all_module_names():
    """
    Looks for all Kleborate modules and returns their names. To qualify as a module, it must be in
    the 'modules' directory, in a subdirectory that matches the filename. For example:
    * modules/contig_stats/contig_stats.py       <- is a module
    * modules/kpsc_mlst/kpsc_mlst.py             <- is a module
    * modules/contig_stats/test_contig_stats.py  <- not a module

    """
    module_dir = pathlib.Path(__file__).parents[0] / 'modules'
    module_names = []
    for module_file in module_dir.glob('*/*.py'):
        dir_name = module_file.parts[-2]
        if module_file.parts[-1][:-3] == dir_name:
            module_names.append(dir_name)
    if 'template' in module_names:
        module_names.remove('template')
    return sorted(module_names)


def get_strain_name(full_path):
    filename = os.path.split(full_path)[1]
    if filename.endswith('_temp_decompress.fasta'):
        filename = filename[:-22]
    if filename.endswith('.gz'):
        filename = filename[:-3]
    return os.path.splitext(filename)[0]

def import_modules():
    """
    This function imports all Kleborate modules (whether or not they are used in this run).
    """
    all_module_names = get_all_module_names()
    modules = {}
    for m in all_module_names:
        modules[m] = importlib.import_module(f'..modules.{m}.{m}', __name__)
    return all_module_names, modules


def check_modules(args, modules, module_names, preset_check_modules, preset_pass_modules):
    """
    This function checks the options, prerequisites external requirements of the used modules. If
    any fail, the program will quit with an error. It returns:
    * A list of all module names. Probably the same list as given, but if any prerequisite modules
      were missing, they've been added to the end.
    * A topologically sorted list of module names. This is the run order so each module will have
      its prerequisites run first.
    * A list of all external programs used.
    """
    all_external_programs = set()

    for m in module_names:
        modules[m].check_cli_options(args)
        all_external_programs.update(modules[m].check_external_programs())

    new_module_names = module_names.copy()
    for m in module_names:
        for prereq in modules[m].prerequisite_modules():
            if prereq not in new_module_names:
                new_module_names.append(prereq)
    module_names = new_module_names

    dependency_graph = {m: modules[m].prerequisite_modules() for m in module_names}
    
    return module_names, get_run_order(dependency_graph), sorted(all_external_programs)


def get_run_order(dependency_graph):
    try:
        ts = graphlib.TopologicalSorter(dependency_graph)
        return list(ts.static_order())
    except graphlib.CycleError:
        sys.exit('Error: module dependency graph contains a cycle')


def check_assembly(assembly):
    """
    This function does a quick check to make sure that the input assembly looks good.
    """
    # for assembly in args.assemblies:
    if os.path.isdir(assembly):
        sys.exit('Error: ' + assembly + ' is a directory (please specify assembly files)')
    if not os.path.isfile(assembly):
        sys.exit('Error: could not find ' + assembly)
    fasta = load_fasta(assembly)
    if len(fasta) < 1:
        sys.exit('Error: invalid FASTA file: ' + assembly)
    for _, seq in fasta:
        if len(seq) == 0:
            sys.exit('Error: invalid FASTA file (contains a zero-length sequence): ' + assembly)


def get_headers(module_names, modules):
    """
    This function returns three lists of headers:
    * top_headers: for Kleborate's output file, shows module names
    * full_headers: for Kleborate's output file, show column names
    * stdout_headers: for Kleborate's stdout

    Each used module contributes headers to these lists. To ensure that there aren't any duplicate
    headers, the module name is added before each header in full_headers and stdout_header,
    separated by a double-underscore.
    """
    _,full_headers, stdout_headers = [''], ['strain'], ['strain']
    for module_name in module_names:
        module_full, module_stdout = modules[module_name].get_headers()
        #top_headers.append(module_name)
        #top_headers += [''] * (len(module_full) - 1)
        full_headers += [f'{module_name}__{h}' for h in module_full]
        stdout_headers += [f'{module_name}__{h}' for h in module_stdout]
    return full_headers, stdout_headers


def gunzip_assembly_if_necessary(assembly, temp_dir):
    if get_compression_type(assembly) == 'gz':
        unzipped_assembly = pathlib.Path(temp_dir) / (uuid.uuid4().hex + '.fasta')
        decompress_file(assembly, unzipped_assembly)
        return unzipped_assembly
    else:
        return assembly


def build_minimap2_index(assembly, unzipped_assembly, external_programs, temp_dir):
    """
    A lot of the modules use minimap2 alignment, so pre-building the index for this assembly once
    can save a bit of time.
    """
    if 'minimap2' not in external_programs:
        return None
    minimap2_index = (pathlib.Path(temp_dir) / (uuid.uuid4().hex + '.mmi')).resolve()
    command = ['minimap2', '-d', minimap2_index, unzipped_assembly]
    p = subprocess.run(command, capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f'\nError: minimap2 failed to index sample {assembly}:\n{p.stderr}')
    return minimap2_index


<<<<<<< HEAD
def decompress_file(in_file, out_file):
    with gzip.GzipFile(in_file, 'rb') as i, open(out_file, 'wb') as o:
        s = i.read()
        o.write(s)
=======
def get_wzi_and_k_locus_results(data_folder, contigs, args):
    seqs = data_folder + '/wzi.fasta'
    database = data_folder + '/wzi.txt'
    bst, _, k_type, _ = mlst_blast(seqs, database, 'yes', [contigs], min_cov=args.min_coverage,
                                   min_ident=args.min_identity, max_missing=0, allow_multiple=False)
    if bst == '0':
        wzi_st = '-'
    else:
        wzi_st = 'wzi' + bst
    if k_type == '':
        k_type = '-'
    return {'wzi': wzi_st,
            'K_locus': k_type}


def get_resistance_results(data_folder, contigs, args, res_headers, kp_complex):
    if not pathlib.Path(contigs).is_file():
        raise OSError
    if args.resistance:
        gene_info, _, _ = read_class_file(data_folder + '/CARD_AMR_clustered.csv')

        # Only do mutation/truncation tests for Kp complex species.
        if kp_complex:
            qrdr = data_folder + '/QRDR_120.fasta'
            trunc = data_folder + '/MgrB_and_PmrB.fasta'
            omp = data_folder + '/OmpK.fasta'
        else:
            qrdr, trunc, omp = None, None, None

        seqs = data_folder + '/CARD_v3.1.13.fasta'
        res_hits = resblast_one_assembly(contigs, gene_info, qrdr, trunc, omp, seqs,
                                         args.min_coverage,  args.min_identity,
                                         args.min_spurious_coverage, args.min_spurious_identity)

        # Double check that there weren't any results without a corresponding output header.
        for h in res_hits.keys():
            if h not in res_headers:
                sys.exit( f'Error: results contained a value ({h}) that is not covered by the '
                          f'output headers')

        res_hits_dict=({r: ';'.join(sorted(res_hits[r])) if r in res_hits else '-' for r in res_headers})

        #Ciprofloxacin resistance prediction

        #Adding aac(6')-Ib-cr to Flq_acquired column
        if "aac(6')-Ib-cr.v1" in res_hits_dict["AGly_acquired"]:
            if ("Flq_acquired", "-") in res_hits_dict.items():
                res_hits_dict["Flq_acquired"]="aac(6')-Ib-cr.v1"
            elif len([res_hits_dict["Flq_acquired"]]) == 1:
                res_hits_dict["Flq_acquired"]=res_hits_dict["Flq_acquired"]+";aac(6')-Ib-cr.v1"
            else:
                res_hits_dict.setdefault("Flq_acquired",[]).append("aac(6')-Ib-cr.v1")
        elif "aac(6')-Ib-cr.v2" in res_hits_dict["AGly_acquired"]:
            if ("Flq_acquired", "-") in res_hits_dict.items():
                res_hits_dict["Flq_acquired"]="aac(6')-Ib-cr.v2"
            elif len([res_hits_dict["Flq_acquired"]]) == 1:
                res_hits_dict["Flq_acquired"]=res_hits_dict["Flq_acquired"]+";aac(6')-Ib-cr.v2"
            else:
                res_hits_dict.setdefault("Flq_acquired",[]).append("aac(6')-Ib-cr.v2")

        # 0 Flq_mutations, 0 Flq_acquired, 0 aac6
            # test: KP_NORM_URN_105939
        if ("Flq_mutations", "-") in res_hits_dict.items() and ("Flq_acquired", "-") in res_hits_dict.items() and "aac(6')-Ib-cr" not in res_hits_dict["AGly_acquired"]:
            res_hits_dict["Cipro_prediction"]="S"
            res_hits_dict["Cipro_prediction_support"]="95"

        # 0 Flq_mutations, 0 Flq_acquired, 1 aac6
            # test: ERR4635459
        elif ("Flq_mutations", "-") in res_hits_dict.items() and ("Flq_acquired", "-") in res_hits_dict.items() and "aac(6')-Ib-cr" in res_hits_dict["AGly_acquired"]:
            res_hits_dict["Cipro_prediction"]="S"
            res_hits_dict["Cipro_prediction_support"]="81"

        # 1 Flq_mutations, 0 Flq_acquired, >=0 aac6
            # test: ERR4046108_NHP1975 (with aac), NK_H5_026 (no aac)
        elif len(set(res_hits_dict["Flq_mutations"].split(";")))==1 and ("Flq_acquired", "-") in res_hits_dict.items():
            res_hits_dict["Cipro_prediction"]="R"
            res_hits_dict["Cipro_prediction_support"]="72"

        # >=2 Flq_mutations, 0 Flq_acquired, >=0 aac6
            # test: G20250619 (with aac), G20250926 (no aac)
        elif len(set(res_hits_dict["Flq_mutations"].split(";")))>=2 and ("Flq_acquired", "-") in res_hits_dict.items():
            res_hits_dict["Cipro_prediction"]="R"
            res_hits_dict["Cipro_prediction_support"]="99"

        # 0 Flq_mutations, 1 Flq_acquired, 0 aac6
            # test: ERR486365
        elif ("Flq_mutations", "-") in res_hits_dict.items() and len(set(res_hits_dict["Flq_acquired"].split(";")))==1  and "aac(6')-Ib-cr" not in res_hits_dict["AGly_acquired"]:
            res_hits_dict["Cipro_prediction"]="R"
            res_hits_dict["Cipro_prediction_support"]="70"

        # 0 Flq_mutations, 1 Flq_acquired, 1 aac6
            # test: HE205
        elif ("Flq_mutations", "-") in res_hits_dict.items() and len(set(res_hits_dict["Flq_acquired"].split(";")))==1  and "aac(6')-Ib-cr" in res_hits_dict["AGly_acquired"]:
            res_hits_dict["Cipro_prediction"]="R"
            res_hits_dict["Cipro_prediction_support"]="93"

        # 0 Flq_mutations, >=2 Flq_acquired, >=0 aac6
            # test: ERR3480591 (with aac), ERR4635120 (no aac)
        elif ("Flq_mutations", "-") in res_hits_dict.items() and len(set(res_hits_dict["Flq_acquired"].split(";")))>=2:
            res_hits_dict["Cipro_prediction"]="R"
            res_hits_dict["Cipro_prediction_support"]="97"

        # >=1 Flq_mutations, >=1 Flq_acquired, >=0 aac6
            # test: ERR3567346_NHP139 (with aac), ERR3585150_BMP790 (no aac)
        elif len(set(res_hits_dict["Flq_mutations"].split(";")))>=1 and len(set(res_hits_dict["Flq_acquired"].split(";")))>=1:
            res_hits_dict["Cipro_prediction"]="R"
            res_hits_dict["Cipro_prediction_support"]="99"

        #return {r: ';'.join(sorted(res_hits[r])) if r in res_hits else '-' for r in res_headers}
        return res_hits_dict
        
    else:
        return {}

>>>>>>> 212fb361b4dc2dc6eda10c33582d38d9a6b132f5

def output_headers(full_headers, stdout_headers, outfile):
    """
    This function prints headers to stdout and writes headers to the output file. Module names are
    trimmed off to make the headers shorter and easier to read.
    """
    trimmed_stdout_headers = [h.split('__')[-1] for h in stdout_headers]
    trimmed_full_headers = [h.split('__')[-1] for h in full_headers]
    print('\t'.join(trimmed_stdout_headers))
    with open(outfile, 'wt') as o:
        o.write('\t'.join(trimmed_full_headers))


def output_results(full_headers, stdout_headers, outfile, results, trim_headers=False):
    """
    This function writes the results to stdout and the output file.
    Always prints stdout headers and writes full headers to the file if the file is new (empty).
    """
    # Print results to the terminal using stdout_headers
    print('\t'.join([str(results.get(x, "-")).strip("[] ") for x in stdout_headers]))

    # Determine headers based on trim_headers option
    headers_to_write = full_headers
    if trim_headers:
        headers_to_write = [h.split('__')[-1] for h in full_headers]

    # Write results to the output file
    with open(outfile, 'at') as o:
        if o.tell() == 0:  # Write headers if file is empty
            o.write('\t'.join(headers_to_write) + '\n')
        o.write('\t'.join([str(results.get(x, "-")).strip("[] ") for x in full_headers]) + '\n')

    # Check for any headers in results that are not in full_headers
    for h in results.keys():
        if h not in full_headers:
            sys.exit(f'Error: results contained a value ({h}) that is not covered by the output headers')


# def output_results(full_headers, stdout_headers, outfile, results):
#     """
#     This function writes the results to stdout and the output file.
#     """
#     print('\t'.join([str(results.get(x, "-")).strip("[] ").replace("assembly", "strain") for x in stdout_headers]))
#     with open(outfile, 'at') as o:
#         if o.tell() > 0:  # Check if the file is not empty
#             o.write('\n')
#         o.write('\t'.join([str(results.get(x, "-")).strip("[] ").replace("assembly", "strain") for x in full_headers]))

#     for h in results.keys():
#         if h not in full_headers:
#             sys.exit(f'Error: results contained a value ({h}) that is not covered by the output headers')


def paper_refs():
    """
    This function prepares the references for the program's help text, with clean line wrapping.
    """
    terminal_width = shutil.get_terminal_size().columns
    text = 'If you use Kleborate, please cite the paper:\n' \
           'Lam MMC, et al. A genomic surveillance framework and genotyping tool for Klebsiella ' \
           'pneumoniae and its related species complex. Nature Communications. 2021. ' \
           'doi:10.1038/s41467-021-24448-3.\n\n' \
           'If you turn on the Kaptive option for full K and O typing, please also cite:\n' \
           'Wyres KL, et al. Identification of Klebsiella capsule synthesis loci from whole ' \
           'genome data. Microbial Genomics. 2016. doi:10.1099/mgen.0.000102.'
    wrapped_text = ''
    for line in text.split('\n'):
        wrapped_text += '\n'.join(textwrap.wrap(line, width=terminal_width - 1))
        wrapped_text += '\n'
    return 'R|' + wrapped_text


def get_version():
    """
    This function returns the version of Kleborate as a string, without the leading 'v'.
    """
    # First try to get the version from the pyproject.toml file, in case this is being run directly
    # from the repo using the kleborate-runner.py file.
    pyproject = pathlib.Path(__file__).parents[1] / 'pyproject.toml'
    if pyproject.is_file():
        with open(pyproject, 'rt') as f:
            pyproject_text = f.read()
        if 'name = "kleborate"' in pyproject_text:  # make sure it's the right pyproject.toml
            match = re.search(r'version = "(\d+\.\d+\.\d+)"', pyproject_text)
            if match:
                return match.group(1)

    # If there wasn't a 'pyproject.toml' file, then Kleborate is probably installed, so use
    # importlib to get the version of the installed package.
    try:
        from importlib import metadata
    except ImportError:
        import importlib_metadata as metadata
    return metadata.version(__package__ or __name__)


if __name__ == '__main__':
    main()
