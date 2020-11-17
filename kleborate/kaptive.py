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

    headers = ['K_locus', 'K_locus_confidence', 'K_locus_problems', 'K_locus_identity',
               'K_locus_missing_genes']
    if locus_type == 'O':
        headers = [x.replace('K_locus', 'O_locus') for x in headers]
        headers.append('O_type')

    if (args.kaptive_k and locus_type == 'K') or (args.kaptive_o and locus_type == 'O'):
        if locus_type == 'K':
            outfile = args.kaptive_k_outfile
        else:  # locus_type == 'O':
            outfile = args.kaptive_o_outfile

        kaptive_results = run_kaptive(kaptive_py, kaptive_db, contigs, outfile, one_thread=False)
        if kaptive_results is None:
            kaptive_results = run_kaptive(kaptive_py, kaptive_db, contigs, outfile, one_thread=True)
        if locus_type == 'O':
            o_locus = kaptive_results[0]
            kaptive_results.append(get_o_type(o_locus))

        assert len(headers) == len(kaptive_results)
        return dict(zip(headers, kaptive_results))
    else:
        return {}


def run_kaptive(kaptive_py, kaptive_db, contigs, output_file, one_thread):
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

        if output_file:	 # if we are saving Kaptive results to file...
            with open(kaptive_table, 'rt') as f:
                kaptive_table_lines = f.readlines()
            assert len(kaptive_table_lines) == 2
            if not os.path.isfile(output_file):
                with open(output_file, 'wt') as f:
                    f.write(kaptive_table_lines[0])	 # write header line
            with open(output_file, 'at') as f:
                f.write(kaptive_table_lines[1])		 # write data line

    if locus is None or confidence is None or problems is None or identity is None:
        sys.exit('Error: Kaptive failed to produce the expected output')

    return [locus, confidence, problems, identity, ','.join(missing)]


def get_o_type(o_locus):
    """
    This function returns an O type using the O locus. In many cases, they are the same, except for:
    * loci O1v1 and O1v2 = type O1
    * loci O2v1 and O2v2 = type O2
    * loci O1/O2v1 and O1/O2v2 = either O1/O2 or 'unknown' (thoughts?)
    """
    if 'O1/O2' in o_locus:
        return 'unknown'
    if 'v1' in o_locus:
        return o_locus.replace('v1', '')
    if 'v2' in o_locus:
        return o_locus.replace('v2', '')
    return o_locus
