"""
Copyright 2025 Mary Maranga
https://github.com/klebgenomics/Kleborate/

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

    if not args:
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    parsed_args = parser.parse_known_args(args)[0]

    help_requested = any(
        flag in args for flag in ('-h', '--help', '--help_all', '--version', '--list_modules')
    )
    if not help_requested and '-o' not in args and '--outdir' not in args:
        sys.exit('Error: The --outdir (-o) option is required. Please specify an output directory.')

    return parser.parse_args(args)


def main():
    all_module_names, modules = import_modules()
    args = parse_arguments(sys.argv[1:], all_module_names, modules)
    print_modules(args, all_module_names, modules)

    module_names, check_module_list, pass_modules = get_used_module_names(args, all_module_names, get_presets())

    preset_check_modules = []
    if args.preset:
        presets = get_presets()
        preset_check_modules = [module for module, _ in presets[args.preset]['check']]

    module_names, module_run_order, external_programs = check_modules(
        args, modules, module_names, check_module_list, pass_modules
    )

    full_headers, stdout_headers = get_headers(module_names, modules)
    print('\t'.join([h.split('__')[-1] for h in stdout_headers]))

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    if args.modules:
        module_name = args.modules.split(',')[0]
        out_files_suffixes = [f'{module_name}_output.txt']
    else:
        out_files_suffixes = [
            'klebsiella_pneumo_complex_output.txt',
            'klebsiella_pneumo_complex_hAMRonization_output.txt',
            'klebsiella_oxytoca_complex_output.txt',
            'escherichia_output.txt'
        ]

    if not args.resume:
        for suffix in out_files_suffixes:
            for file in glob(f'{args.outdir}/*{suffix}'):
                os.remove(file)

    for assembly in args.assemblies:
        check_assembly(assembly)

        with tempfile.TemporaryDirectory() as temp_dir:
            unzipped_assembly = gunzip_assembly_if_necessary(assembly, temp_dir)
            minimap2_index = build_minimap2_index(assembly, unzipped_assembly, external_programs, temp_dir)
            results = {'strain': get_strain_name(assembly)}

            pass_check = True  # default, assume no check and run all modules
            if args.preset and len(check_module_list) > 0:
                for module, check in presets[args.preset]['check']:
                    try:
                        module_results = modules[module].get_results(
                            unzipped_assembly, minimap2_index, args, results
                        )
                        results.update({f'{module}__{header}': result for header, result in module_results.items()})
                        check_function = globals()[check]

                        if not check_function(module_results):
                            pass_check = False
                            print(f"Assembly {assembly} failed in check {check}.")
                            break

                    except Exception as e:
                        print(f"Error encountered while processing {assembly} with {module}: {e}.")
                        pass_check = False
                        break

            if pass_check:
                for module in module_run_order:
                    if module not in preset_check_modules:
                        module_results = modules[module].get_results(
                            unzipped_assembly, minimap2_index, args, results
                        )
                        results.update({f'{module}__{header}': result for header, result in module_results.items()})
            else:
                for module in module_run_order:
                    if module not in preset_check_modules:
                        module_headers = [header for header in full_headers if header.startswith(module)]
                        for header in module_headers:
                            results[header] = 'Not Tested'

            if args.modules:
                module_name = args.modules.split(',')[0]
                outfile_suffix = f'{module_name}_output.txt'
                output_file = os.path.join(args.outdir, outfile_suffix)
                filtered_headers = [h for h in full_headers if h.split('__')[-1] not in annotation_fields]
                filtered_results = {k: v for k, v in results.items() if k in filtered_headers}
                output_results(filtered_headers, stdout_headers, output_file, filtered_results, args.trim_headers)
            else:
                species = results.get('enterobacterales__species__species', None)
                if species and is_kp_complex({'species': species}):
                    harmonization_file = os.path.join(args.outdir, 'klebsiella_pneumo_complex_hAMRonization_output.txt')
                    klebsiella_pneumo_file = os.path.join(args.outdir, 'klebsiella_pneumo_complex_output.txt')

                    output_results_klebsiella_pneumo_complex_hAMRonization(
                        full_headers, stdout_headers, harmonization_file, results, args.trim_headers
                    )

                    selective_headers = [
                        header for header in full_headers
                        if not header.startswith('klebsiella_pneumo_complex__amr') or
                        header.split('__')[-1] in res_headers
                    ]
                    filtered_results = {header: results.get(header, "-") for header in selective_headers}
                    output_results(selective_headers, stdout_headers, klebsiella_pneumo_file, filtered_results, args.trim_headers)

                else:
                    if species and is_ko_complex({'species': species}):
                        outfile_suffix = 'klebsiella_oxytoca_complex_output.txt'
                    elif species and is_escherichia({'species': species}):
                        outfile_suffix = 'escherichia_output.txt'
                    else:
                        print(f"Assembly {assembly} does not match any specified species. Skipping to next assembly.")
                        continue

                    output_file = os.path.join(args.outdir, outfile_suffix)
                    output_results(full_headers, stdout_headers, output_file, results, args.trim_headers)


# def main():
#     all_module_names, modules = import_modules()
#     args = parse_arguments(sys.argv[1:], all_module_names, modules)
#     print_modules(args, all_module_names, modules)

#     module_names, check_module_list, pass_modules = get_used_module_names(args, all_module_names, get_presets())

#     preset_check_modules = []
#     if args.preset:
#         presets = get_presets()
#         preset_check_modules = [module for module, _ in presets[args.preset]['check']]

#     module_names, module_run_order, external_programs = check_modules(args, modules, module_names, check_module_list, pass_modules)

#     full_headers, stdout_headers = get_headers(module_names, modules)
#     print('\t'.join([h.split('__')[-1] for h in stdout_headers]))

#     if not os.path.exists(args.outdir):
#         os.makedirs(args.outdir)

#     # If the resume flag is not set, remove existing output files
#     if args.modules:
#         module_name = args.modules.split(',')[0]
#         out_files_suffixes = [f'{module_name}_output.txt']
#     else:
#         out_files_suffixes = [
#             'klebsiella_pneumo_complex_output.txt',
#             'klebsiella_pneumo_complex_hAMRonization_output.txt',
#             'klebsiella_oxytoca_complex_output.txt',
#             'escherichia_output.txt'
#         ]

#     if not args.resume:
#         for suffix in out_files_suffixes:
#             for file in glob(f'{args.outdir}/*{suffix}'):
#                 os.remove(file)

#     for assembly in args.assemblies:
#         check_assembly(assembly)

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
#                             print(f"Assembly {assembly} failed in check {check}.")
#                             break # Exit the for loop since this assembly failed the check

#                     except Exception as e:
#                         print(f"Error encountered while processing {assembly} with {module}: {e}.")
#                         pass_check = False
#                         break

#             # proceed through all other modules
#             if pass_check:
#                 for module in module_run_order:
#                     if module not in preset_check_modules:
#                         module_results = modules[module].get_results(unzipped_assembly, minimap2_index, args, results)
#                         results.update({f'{module}__{header}': result for header, result in module_results.items()})
#             else:
#                 # Populate results with "Not Tested" for modules that did not run
#                 for module in module_run_order:
#                     if module not in preset_check_modules:
#                         module_headers = [header for header in full_headers if header.startswith(module)]
#                         for header in module_headers:
#                             results[header] = 'Not Tested'

#             species = results.get('enterobacterales__species__species', None)
#             if species and is_kp_complex({'species': species}):
#                 harmonization_file = os.path.join(args.outdir, 'klebsiella_pneumo_complex_hAMRonization_output.txt')
#                 klebsiella_pneumo_file = os.path.join(args.outdir, 'klebsiella_pneumo_complex_output.txt')

#                 # Write all results for harmonization file with required processing
#                 output_results_klebsiella_pneumo_complex_hAMRonization(full_headers, stdout_headers, harmonization_file, results, args.trim_headers)


#                 # Select headers for klebsiella_pneumo_complex_output.txt
#                 selective_headers = [
#                     header for header in full_headers 
#                     if not header.startswith('klebsiella_pneumo_complex__amr') or 
#                     header.split('__')[-1] in res_headers
#                 ]

#                 # Filter the results directly based on selective headers
#                 filtered_results = {header: results.get(header, "-") for header in selective_headers}

#                 # Output filtered results
#                 output_results(selective_headers, stdout_headers, klebsiella_pneumo_file, filtered_results, args.trim_headers)

#             else:
#                 if species and is_ko_complex({'species': species}):
#                     outfile_suffix = 'klebsiella_oxytoca_complex_output.txt'
#                 elif species and is_escherichia({'species': species}):
#                     outfile_suffix = 'escherichia_output.txt'
#                 else:
#                     print(f"Assembly {assembly} does not match any specified species. Skipping to next assembly.")
#                     continue

#                 output_file = os.path.join(args.outdir, outfile_suffix)
#                 output_results(full_headers, stdout_headers, output_file, results, args.trim_headers)






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
            'klebsiella_pneumo_complex__resistance_gene_count', 'klebsiella_pneumo_complex__cipro_prediction', 'klebsiella_pneumo_complex__wzi','klebsiella_pneumo_complex__kaptive'
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
            'escherichia__mlst_achtman', 'escherichia__mlst_pasteur', 'escherichia__pathovar', 'escherichia__mlst_lee', 'escherichia__ezclermont','escherichia__stxtyper', 'escherichia__ectyper', 'escherichia__amr'
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

        if args.preset == 'escherichia' and args.modules is None:
            pass_modules = [module for module in pass_modules if module != 'escherichia__mlst_pasteur']

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


def decompress_file(in_file, out_file):
    with gzip.GzipFile(in_file, 'rb') as i, open(out_file, 'wb') as o:
        s = i.read()
        o.write(s)


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


res_headers = [
    'AGly_acquired', 'Col_acquired', 'Fcyn_acquired', 'Flq_acquired', 'Gly_acquired', 
    'MLS_acquired', 'Phe_acquired', 'Rif_acquired', 'Sul_acquired', 'Tet_acquired', 
    'Tgc_acquired', 'Tmt_acquired', 'Bla_acquired', 'Bla_inhR_acquired', 'Bla_ESBL_acquired', 
    'Bla_ESBL_inhR_acquired', 'Bla_Carb_acquired', 'Bla_chr', 'SHV_mutations', 
    'Omp_mutations', 'Col_mutations', 'Flq_mutations', 'truncated_resistance_hits', 
    'spurious_resistance_hits'
]


annotation_fields = [
        'Genetic_variation_type','Drug_class','Input_sequence_ID','Input_gene_length', 'Input_gene_start', 'Input_gene_stop', 'Reference_gene_length',
        'Reference_gene_start', 'Reference_gene_stop', 'Sequence_identity', 'Coverage','Reference_accession','Strand_orientation',
        'Software_name', 'Software_version', 'Reference_database_name',
        'Reference_database_version','Input_protein_length','Reference_protein_length','Input_protein_start', 'Input_protein_stop','Antimicrobial_agent', 
        'Coverage_depth', 'Coverage_ratio','Predicted_phenotype','predicted_phenotype_confidence_level', 
        'Reference_protein_start', 'Reference_protein_stop','Resistance_mechanism'
]


def output_results_klebsiella_pneumo_complex_hAMRonization(full_headers, stdout_headers, outfile, results, trim_headers=False):
    # Rename strain to Input_file_name
    input_file_name = results['strain']

    res_headers = [
        'AGly_acquired', 'Col_acquired', 'Fcyn_acquired', 'Flq_acquired', 'Gly_acquired',
        'MLS_acquired', 'Phe_acquired', 'Rif_acquired', 'Sul_acquired', 'Tet_acquired',
        'Tgc_acquired', 'Tmt_acquired', 'Bla_acquired', 'Bla_inhR_acquired',
        'Bla_ESBL_acquired', 'Bla_ESBL_inhR_acquired', 'Bla_Carb_acquired', 'Bla_chr', 'truncated_resistance_hits',
        'spurious_resistance_hits'
    ]

    mutation_variant_headers = [
        'SHV_mutations', 'Omp_mutations', 'Col_mutations', 'Flq_mutations'
    ]

    prefix = 'klebsiella_pneumo_complex__amr__'

    def clean(val):
    # strip any “-digits%” suffix
        val = re.sub(r'-\d+%$', '', val)
        return val.translate(str.maketrans('', '', '^*?')).strip()

    # Parse: build a dict of {field: [list of entries]}
    parsed = {}
    for k, v in results.items():
        if isinstance(v, list) and v:
            val = v[0]
            if isinstance(val, str) and val != '-':
                parsed[k] = [clean(x) for x in val.split(';')]
            else:
                parsed[k] = []
        else:
            parsed[k] = []

    rows = []

    software_name = results.get('Software_name', ['Kleborate'])[0]
    software_version = results.get('Software_version', ['3.1.3'])[0]
    db_name = results.get('Reference_database_name', ['CARD'])[0]
    db_version = results.get('Reference_database_version', ['3.2.9'])[0]

    # ----- ACQUIRED GENES -----
    for header in res_headers:
        key = prefix + header
        variants = parsed.get(key, [])
        for variant in variants:
            if not variant or variant == '-':
                continue

            row = {
                'Input_file_name': input_file_name,
                'Gene_symbol': variant,
                'Mutation': '-',
                'Software_name': software_name,
                'Software_version': software_version,
                'Reference_database_name': db_name,
                'Reference_database_version': db_version
            }

            # set annottaion field for each variant
            for field in annotation_fields:
                ann_key = prefix + field
                ann_values = parsed.get(ann_key, [])
                matched_val = '-'
                for item in ann_values:
                    # Custom parsing for Reference_accession to handle nested colons
                    if field == 'Reference_accession' and 'ARO:' in item:
                        idx = item.find('ARO:')
                        g = item[:idx-1]
                        v = item[idx:]
                    else:
                        if ':' in item:
                            idx = item.rfind(':')
                            g = item[:idx]
                            v = item[idx+1:]
                        else:
                            g = item
                            v = item
                    if clean(g) == variant:
                        matched_val = v
                        break
                row[field] = matched_val

            if header == 'Flq_acquired':
                row['Drug_class'] = 'Fluoroquinolone antibiotic'

            rows.append(row)

    # ----- MUTATIONS-----
    for header in mutation_variant_headers:
        key = prefix + header
        variants = parsed.get(key, [])
        for variant in variants:
            if not variant or variant == '-':
                continue

            # Extract gene symbol and mutation after first colon
            if ':' in variant:
                Gene_symbol, mutation = variant.split(':', 1)
            elif '-' in variant:
                Gene_symbol = variant.split('-')[0]
                mutation = variant
            else:
                Gene_symbol = variant
                mutation = variant

            row = {
                'Input_file_name': input_file_name,
                'Gene_symbol': Gene_symbol,
                'Mutation': mutation,
                'Software_name': software_name,
                'Software_version': software_version,
                'Reference_database_name': db_name,
                'Reference_database_version': db_version
            }

            # add annotation field for each variant
            for field in annotation_fields:
                ann_key = prefix + field
                ann_values = parsed.get(ann_key, [])
                matched_val = '-'
                for item in ann_values:
                    if field == 'Reference_accession' and 'ARO:' in item:
                        idx = item.find('ARO:')
                        g = item[:idx-1]
                        v = item[idx:]
                    else:
                        if ':' in item:
                            idx = item.rfind(':')
                            g = item[:idx]
                            v = item[idx+1:]
                        else:
                            g = item
                            v = item
                    if clean(g) == variant:
                        matched_val = v
                        break
                row[field] = matched_val
            rows.append(row)

    # set software and database information
    for row in rows:
        if row['Software_name'] in ['-', ''] or row['Software_version'] in ['-', ''] \
           or row['Reference_database_name'] in ['-', ''] or row['Reference_database_version'] in ['-', '']:
            row['Software_name'] = software_name
            row['Software_version'] = software_version
            row['Reference_database_name'] = db_name
            row['Reference_database_version'] = db_version

    headers = ['Input_file_name', 'Gene_symbol', 'Mutation'] + annotation_fields

    with open(outfile, 'a') as f:
        # Only write headers if the file is empty
        if f.tell() == 0:
            f.write('\t'.join(headers) + '\n')
        for row in rows:
            f.write('\t'.join([str(row.get(h, '-')) for h in headers]) + '\n')


def output_results(full_headers, stdout_headers, outfile, results, trim_headers=False):
    """
    This function writes the results to stdout and the output file.
    Always prints stdout headers and writes full headers to the file if the file is new (empty).
    """
    # Print results to the terminal using stdout_headers
    print('\t'.join([
        str(results.get(x, "-")).strip("[]'\"") if not isinstance(results.get(x, "-"), list) 
        else ";".join(map(str, results.get(x, "-")))
        for x in stdout_headers
    ]))

    # Determine headers based on trim_headers option
    headers_to_write = full_headers
    if trim_headers:
        headers_to_write = [h.split('__')[-1] for h in full_headers]

    # Write results to the output file
    with open(outfile, 'at') as o:
        if o.tell() == 0:  # Write headers if file is empty
            o.write('\t'.join(headers_to_write) + '\n')
        o.write('\t'.join([
            str(results.get(x, "-")).strip("[]'\"") if not isinstance(results.get(x, "-"), list) 
            else ";".join(map(str, results.get(x, "-")))
            for x in full_headers
        ]) + '\n')

    # Check for any headers in results that are not in full_headers
    for h in results.keys():
        if h not in full_headers:
            sys.exit(f'Error: results contained a value ({h}) that is not covered by the output headers')

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
