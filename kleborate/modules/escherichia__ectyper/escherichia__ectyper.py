"""
Copyright 2025 Mary Maranga (gathonimaranga@gmail.com)
https://github.com/klebgenomics/Kleborate

This file is part of Kleborate. Kleborate is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Kleborate is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Kleborate. If
not, see <https://www.gnu.org/licenses/>.
"""

import os
import sys
import shutil
from pathlib import Path
import pandas as pd
import tempfile
from collections import defaultdict
import subprocess
import csv
import re


def description():
    return 'E. coli serotyping of O and H antigens'


def prerequisite_modules():
    return []


def get_headers():

    full_headers = [
        'O-type','H-type','Serotype','QC','Evidence','GeneScores',
        'AlleleKeys','GeneIdentities(%)','GeneCoverages(%)','GeneLengths','Warnings'
    ]
    stdout_headers = []
    
    return full_headers, stdout_headers


def add_cli_options(parser):
    module_name = os.path.basename(__file__)[:-3]
    group = parser.add_argument_group(f'{module_name} module')
    group.add_argument('--escherichia__serotyping_cores', type=int, default=4,
                       help='Number of CPU cores to use for ectyper (default: 4)')
    return


    
def check_cli_options(args):
    """
    Validate the command-line arguments.
    """

    if args.escherichia__serotyping_cores < 4:
        raise ValueError("The number of threads must be at least 4")

    if not shutil.which('ectyper'):
        sys.exit('Error: ectyper is not installed or not in PATH.')


def check_external_programs():
    """
    Ensure the required external programs are available.
    """
    if not shutil.which('ectyper'):
        sys.exit('Error: could not find ectyper executable.')
    return ['ectyper']


def data_dir():
    return pathlib.Path(__file__).parents[0] / 'data'


def run_ectyper(input_fasta, output_dir, quiet,cores) :
    """
    Run ectyper to serotype E. coli from a FASTA file.

    Parameters:
        input_fasta (str): Path to the input FASTA file.
        output_dir (str): Directory where ectyper will write its output.
        quiet (bool): If True, suppress ectyper stdout/stderr.

    Raises:
        subprocess.CalledProcessError: if ectyper exits with non-zero.
    """


    os.makedirs(output_dir, exist_ok=True)
    command = [
        "ectyper",
        "-i", input_fasta,
        "-o", output_dir,
        "-c", str(cores)
    ]
    # Choose output redirection
    if quiet:
        stdout_dest = subprocess.DEVNULL
        stderr_dest = subprocess.DEVNULL
    else:
        stdout_dest = subprocess.PIPE
        stderr_dest = subprocess.PIPE

    subprocess.run(
        command,
        stdout=stdout_dest,
        stderr=stderr_dest,
        check=True,
        text=not quiet
    )

def get_results(assembly, minimap2_index, args, previous_results):
    results = {}
    quiet = getattr(args, 'quiet', True)

    with tempfile.TemporaryDirectory() as tmpdir:
        run_ectyper(assembly, tmpdir, quiet=quiet, cores=args.escherichia__serotyping_cores)

        tsv_path = os.path.join(tmpdir, "output.tsv")
        full_headers, _ = get_headers()

        if not os.path.exists(tsv_path):
            if not quiet:
                print(f"[get_results] output.tsv not found in temporary directory")
            return {}

        with open(tsv_path, newline='') as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            try:
                row = next(reader)
            except StopIteration:
                return {}

            for col in full_headers:
                raw = row.get(col, "") or ""

                # General cleaning:
                # 1. Remove "-:" at the beginning
                # 2. Remove trailing semicolons
                # 3. Strip leading/trailing whitespace
                raw = re.sub(r"^-:", "", raw)
                clean = raw.rstrip(";").strip()

                results[col] = clean

    return results




# from .ectyper import load_json_file, alignment, minimap2_output_to_df, ectyper_dict_to_df, setAlleleMeta, get_prediction, predict_serotype, process_sample_data

# LOG = logging.getLogger(__name__)
# LOG.setLevel(logging.DEBUG)


# def description():
#     return 'E. coli serotyping of O and H antigens'


# def prerequisite_modules():
#     return []


# def get_headers():
#     full_headers = ['Otype', 'Htype', 'Serotype', 'Evidence', 'AlleleKeys',
#     'GeneIdentities'] 
#     stdout_headers = []
#     return full_headers, stdout_headers


# def add_cli_options(parser):
#     module_name = os.path.basename(__file__)[:-3]
#     group = parser.add_argument_group(f'{module_name} module')
#     group.add_argument('--escherichia__serotyping_min_identity', type=float, default=90.0,
#                        help='Minimum alignment percent identity for detecting serotypes')
#     group.add_argument('--escherichia__serotyping_min_coverage', type=float, default=80.0,
#                        help='Minimum alignment percent coverage for detecting serotypes')
#     # New arguments
#     group.add_argument('--escherichia__serotyping_percentIdentityOtype', type=float, default=90.0,
#                        help='Minimum percent identity for O-type detection')
#     group.add_argument('--escherichia__serotyping_percentIdentityHtype', type=float, default=95.0,
#                        help='Minimum percent identity for H-type detection')
#     group.add_argument('--escherichia__serotyping_percentCoverageOtype', type=float, default=90.0,
#                        help='Minimum percent coverage for O-type detection')
#     group.add_argument('--escherichia__serotyping_percentCoverageHtype', type=float, default=50.0,
#                        help='Minimum percent coverage for H-type detection')
#     group.add_argument('--escherichia__serotyping_HIGH_SIMILARITY_THRESHOLD_O', type=float, default=0.00771,
#                        help='High similarity threshold for O-type')
#     group.add_argument('--escherichia__serotyping_MIN_O_IDENTITY_LS', type=float, default=95.0,
#                        help='Minimum percent identity for O-type (LS threshold)')
#     group.add_argument('--escherichia__serotyping_MIN_O_COVERAGE_LS', type=float, default=48.0,
#                        help='Minimum percent coverage for O-type (LS threshold)')
#     return


# def check_cli_options(args):
#     if not (50.0 <= args.escherichia__serotyping_min_identity <= 100.0):
#         sys.exit('Error: --escherichia__serotyping_min_identity must be between 50.0 and 100.0')
#     if not (50.0 <= args.escherichia__serotyping_min_coverage <= 100.0):
#         sys.exit('Error: --escherichia__serotyping_min_coverage must be between 50.0 and 100.0')

#     # Additional checks for new arguments
#     if not (50.0 <= args.escherichia__serotyping_percentIdentityOtype <= 100.0):
#         sys.exit('Error: --percentIdentityOtype must be between 50.0 and 100.0')
#     if not (50.0 <= args.escherichia__serotyping_percentIdentityHtype <= 100.0):
#         sys.exit('Error: --percentIdentityHtype must be between 50.0 and 100.0')
#     if not (50.0 <= args.escherichia__serotyping_percentCoverageOtype <= 100.0):
#         sys.exit('Error: --percentCoverageOtype must be between 50.0 and 100.0')
#     if not (50.0 <= args.escherichia__serotyping_percentCoverageHtype <= 100.0):
#         sys.exit('Error: --percentCoverageHtype must be between 50.0 and 100.0')
#     if not (50.0 <= args.escherichia__serotyping_MIN_O_IDENTITY_LS <= 100.0):
#         sys.exit('Error: --MIN_O_IDENTITY_LS must be between 50.0 and 100.0')
#     if not (10.0 <= args.escherichia__serotyping_MIN_O_COVERAGE_LS <= 100.0):
#         sys.exit('Error: --MIN_O_COVERAGE_LS must be between 10.0 and 100.0')

#     return

    

# def check_external_programs():
#     if not shutil.which('minimap2'):
#         sys.exit('Error: could not find minimap2')
#     return ['minimap2']



# def data_dir():
#     return pathlib.Path(__file__).parents[0] / 'data'



# def get_results(assembly, minimap2_index, args, previous_results):
#     # Step 1: Load the JSON database
#     ref_db = data_dir() / 'ectyper_alleles_db.json'
#     ref_file = data_dir() / 'ectyper_alleles_db.fasta'
#     ectyper_dict = load_json_file(ref_db)
#     if ectyper_dict is None:
#         sys.exit(f"Error: Could not load the JSON reference database from {ref_db}")

#     # Step 2: Perform the alignment
#     alignment_hits = alignment(
#         ref_file,
#         minimap2_index,
#         assembly,
#         args.escherichia__serotyping_min_identity,
#         args.escherichia__serotyping_min_coverage
#     )

#     # Step 3: Predict the serotype
#     predictions_dict, output_df = predict_serotype(alignment_hits, ectyper_dict, args, previous_results)

#     # Step 4: Process each data for final serotype prediction
#     final_results_dict = {}
#     for genome_name in predictions_dict:
#         final_results_dict[genome_name] = process_sample_data(genome_name, predictions_dict, ectyper_dict)

#     # Step 5: Generate the output dictionary including the full headers
#     combined_results = {}  # Initialize combined_results dictionary
#     output_results = []
#     for genome_name, result in final_results_dict.items():
#         # Extract the necessary values based on headers
#         otype = result.get('Otype', '-')
#         htype = result.get('Htype', '-')
#         serotype = result.get('Serotype', '-')
#         # gene_scores = result.get('gene_scores', '-')
#         evidence = result.get('O_allele_No', ['-'])[0]
#         allele_keys = result.get('alleles', '-')
#         # qc = result.get('QC', '-')

#         # Gene identities, contig names, coordinates, and lengths
#         gene_identities = []
#         gene_contig_names = []
#         gene_ranges = []
#         gene_lengths = []

#         # Iterate over all genes
#         for gene, gene_data in result.items():
#             if isinstance(gene_data, dict):
#                 gene_identities.append(str(gene_data.get('identity', '-')))
#                 gene_contig_names.append(gene_data.get('contigname', '-'))
#                 gene_ranges.append(gene_data.get('coordinates', '-'))
#                 gene_lengths.append(str(gene_data.get('length', '-')))

#         # If no gene data is found, append '-'
#         if not gene_identities:
#             gene_identities.append('-')
#         if not gene_contig_names:
#             gene_contig_names.append('-')
#         if not gene_ranges:
#             gene_ranges.append('-')
#         if not gene_lengths:
#             gene_lengths.append('-')

#         # Join the lists with semicolons
#         gene_identities_str = ';'.join(gene_identities)
#         gene_contig_names_str = ';'.join(gene_contig_names)
#         gene_ranges_str = ';'.join(gene_ranges)
#         gene_lengths_str = ';'.join(gene_lengths)

#         results = {
#             'Otype': otype,
#             'Htype': htype,
#             'Serotype': serotype,
#             # 'gene_scores': gene_scores,
#             'Evidence': evidence,
#             'AlleleKeys': allele_keys,
#             'GeneIdentities': gene_identities_str,
#             # 'GeneContigNames': gene_contig_names_str,
#             # 'Gene ranges': gene_ranges_str,
#             # 'GeneLengths': gene_lengths_str
#             # 'QC':qc
#         }

#         full_headers, _ = get_headers()
#         for h in results.keys():
#             if h not in full_headers:
#                 sys.exit(f'Error: results contained a value ({h}) that is not covered by the full headers')

#         combined_results.update({
#             r: ';'.join(sorted(results[r].split(';'))) if isinstance(results[r], str) and r in [
#             'GeneIdentities' 
#             ] else results[r] for r in full_headers
#         })


#         # combined_results.update({
#         #     r: ';'.join(sorted(results[r].split(';'))) if isinstance(results[r], str) and r in [
#         #         'GeneIdentities', 'GeneContigNames', 'Gene ranges', 'GeneLengths'
#         #     ] else results[r] for r in full_headers
#         # })

#     return combined_results



# def get_results(assembly, minimap2_index, args, previous_results):
#     """
#     Run ectyper on the given assembly in a temporary directory, parse output.tsv,
#     then clean up the temporary files.

#     Parameters:
#         assembly (str): Path to the FASTA assembly.

#     Returns:
#         dict: Mapping each header from get_headers() to its value in output.tsv.
#     """
#     results = {}
#     # Determine whether to suppress ectyper output
#     quiet = getattr(args, 'quiet', True)

#     # Create a temporary directory for ectyper output
#     with tempfile.TemporaryDirectory() as tmpdir:
#         # Run ectyper; it will produce output.tsv in tmpdir
#         run_ectyper(assembly, tmpdir, quiet=quiet)

#         # Path to the results TSV
#         tsv_path = os.path.join(tmpdir, "output.tsv")
#         full_headers, _ = get_headers()

#         if not os.path.exists(tsv_path):
#             if not quiet:
#                 print(f"[get_results] output.tsv not found in temporary directory")
#             return {}

#         with open(tsv_path, newline='') as fh:
#             reader = csv.DictReader(fh, delimiter='\t')
#             try:
#                 row = next(reader)
#             except StopIteration:
#                 return {}

#             for col in full_headers:
#                 results[col] = row.get(col)

#     # Temporary directory and its contents are automatically removed here
#     return results




