"""
Copyright 2020 Kat Holt
Copyright 2020 Ryan Wick (rrwick@gmail.com)
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
import sys
import tempfile


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


def get_kaptive_results(locus_type, kaptive_py, kaptive_db, contigs, args):
    assert locus_type == 'K' or locus_type == 'O'

    headers = ['K_locus', 'K_type', 'K_locus_confidence', 'K_locus_problems', 'K_locus_identity',
               'K_locus_missing_genes']
    if locus_type == 'O':
        headers = [x.replace('K_', 'O_') for x in headers]

    if (args.kaptive_k and locus_type == 'K') or (args.kaptive_o and locus_type == 'O'):
        if locus_type == 'K':
            outfile = args.kaptive_k_outfile
        else:  # locus_type == 'O':
            outfile = args.kaptive_o_outfile

        # Try with multiple threads, and if that fails (possible due to a multithreading BLAST
        # error), try again with one thread.
        kaptive_results = run_kaptive(kaptive_py, kaptive_db, contigs, outfile,
                                      args.min_kaptive_confidence, one_thread=False)
        if kaptive_results is None:
            kaptive_results = run_kaptive(kaptive_py, kaptive_db, contigs, outfile,
                                          args.min_kaptive_confidence, one_thread=True)

        assert len(headers) == len(kaptive_results)
        return dict(zip(headers, kaptive_results))
    else:
        return {}


def run_kaptive(kaptive_py, kaptive_db, contigs, output_file, min_confidence, one_thread):
    thread_option = ' --threads 1' if one_thread else ''
    with tempfile.TemporaryDirectory() as tmp_dir:
        kaptive_prefix = tmp_dir + '/kaptive'
        kaptive_table = kaptive_prefix + '_table.txt'

        p = subprocess.Popen(sys.executable + ' ' + kaptive_py +
                             ' -a ' + contigs +
                             ' -k ' + kaptive_db +
                             ' -o ' + kaptive_prefix +
                             ' --verbose --no_seq_out --no_json' +
                             thread_option,
                             shell=True,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()

        # Make sure the output is a string, whether we are in Python 2 or 3.
        if not isinstance(stdout, str):
            stdout = stdout.decode()
        if not isinstance(stderr, str):
            stderr = stderr.decode()

        # If we hit the BLAST threading problem, return None and we'll try again with one thread.
        if 'tblastn crashed!' in stderr and not one_thread:
            return None

        if p.returncode != 0:
            if stderr:
                sys.exit('Error: Kaptive failed to run with the following error:\n' +
                         stderr.strip())
            else:
                sys.exit('Error: Kaptive failed to run')

        locus, serotype, confidence, problems, identity = None, None, None, None, None
        missing = []

        # Parse the required information from the Kaptive verbose output.
        output_lines = stdout.splitlines()
        missing_gene_lines = False
        for line in output_lines:
            if 'Best match locus:' in line:
                locus = line.split('Best match locus:')[1].strip()
            if 'Best match type:' in line:
                serotype = line.split('Best match type:')[1].strip()
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

        if output_file:	 # if we are saving Kaptive results to file...
            with open(kaptive_table, 'rt') as f:
                kaptive_table_lines = f.readlines()
            assert len(kaptive_table_lines) == 2
            if not os.path.isfile(output_file):
                with open(output_file, 'wt') as f:
                    f.write(kaptive_table_lines[0])	 # write header line
            with open(output_file, 'at') as f:
                f.write(kaptive_table_lines[1])		 # write data line

    if locus is None or serotype is None or confidence is None or problems is None or identity is None:
        sys.exit('Error: Kaptive failed to produce the expected output')

    if not confidence_meets_threshold(confidence, min_confidence):
        locus = f'unknown ({locus})'
        if not serotype.startswith('unknown'):
            serotype = f'unknown ({serotype})'

    return [locus, serotype, confidence, problems, identity, ','.join(missing)]


def confidence_meets_threshold(confidence, min_confidence):
    """
    Returns True if the confidence level meets or exceeds the minimum confidence level.
    """
    min_confidence = min_confidence.replace('_', ' ')
    scores = {'None': 0, 'Low': 1, 'Good': 2, 'High': 3, 'Very high': 4, 'Perfect': 5}
    return scores[confidence] >= scores[min_confidence]
