#!/usr/bin/env python3

import os
import time
import subprocess
import argparse
import copy
from ftplib import FTP
from io import BytesIO

def main():
    args = get_arguments()
    try:
        os.makedirs(args.assembly_dir)
    except OSError:
        pass
    os.chdir(args.assembly_dir)

    ncbi_header, assembly_data = download_klebs_metadata()

    assemblies_not_yet_downloaded = [x for x in assembly_data if not file_exists(x[7])]
    print(str(len(assemblies_not_yet_downloaded)) + ' out of ' + str(len(assembly_data)) + ' will be downloaded')
    for sample in assemblies_not_yet_downloaded:
        full_header = copy.deepcopy(ncbi_header)
        sample_data = [sample[i] for i in range(len(ncbi_header))]
        temp_file_name, final_file_name = download_genome(sample[7])
        print()
        print(final_file_name, flush=True)

        print('  Running Kleborate...', flush=True)
        run_kleborate(temp_file_name, sample_data, args.kleborate_dir, full_header)

        print('  Running Kaptive...', flush=True)
        run_kaptive(temp_file_name, sample_data, args.kaptive_dir, full_header)

        # Save data to file.
        if not os.path.isfile(args.table):
            with open(args.table, 'wt') as table:
                table.write('\t'.join(full_header))
                table.write('\n')
        with open(args.table, 'at') as table:
            table.write('\t'.join(sample_data))
            table.write('\n')

        # Now that we're done, give the file its final name. We don't do this sooner in case the
        # program was interrupted - the final assembly file should only exist if Kleborate and
        # Kaptive finished and the data was saved to the table.
        os.rename(temp_file_name, final_file_name)
        print('  Done!', flush=True)
    
    ftp.quit()

def connect_ftp():
    # print('Connecting to ftp.ncbi.nlm.nih.gov...', flush=True)  # TEMP
    ftp = FTP('ftp.ncbi.nlm.nih.gov', timeout=1)
    # print('Logging in...', flush=True)  # TEMP
    ftp.login()
    return ftp

def get_arguments():
    """
    Parse the command line arguments.
    """
    parser = argparse.ArgumentParser(description='NCBI Pathogen Detection Klebsiella downloader')
    parser.add_argument('--assembly_dir', required=True)
    parser.add_argument('--table', required=True)
    parser.add_argument('--kleborate_dir', required=True)
    parser.add_argument('--kaptive_dir', required=True)
    args = parser.parse_args()
    args.assembly_dir = os.path.abspath(args.assembly_dir)
    args.table = os.path.abspath(args.table)
    args.kleborate_dir = os.path.abspath(args.kleborate_dir)
    args.kaptive_dir = os.path.abspath(args.kaptive_dir)
    return args


def download_klebs_metadata():
    ftp = connect_ftp()
    path = 'pathogen/Results/Klebsiella/latest_kmer/Metadata/'
    metadata_file = [x for x in ftp.nlst(path) if x.endswith('.tsv')][0]
    r = BytesIO()
    print('Downloading metadata...', flush=True)
    ftp.retrbinary('RETR ' + metadata_file, r.write)
    ftp.quit()
    all_lines = r.getvalue().decode().splitlines()
    header = all_lines[0].strip().split('\t')
    data = [x.strip().split('\t') for x in all_lines[1:]]
    print('Found ' + str(len(data)) + ' total records')
    assembly_data = [x for x in data if x[7] != 'NULL']
    print(str(len(assembly_data)) + ' of which are assemblies', flush=True)
    return header, assembly_data


def download_genome(assembly_accession):
    ftp = connect_ftp()
    accession_part_1 = assembly_accession.split('_')[0]
    accession_end = assembly_accession.split('_')[1]
    accession_part_2 = accession_end[0:3]
    accession_part_3 = accession_end[3:6]
    accession_part_4 = accession_end[6:9]
    path = 'genomes/all/' + accession_part_1 + '/' + accession_part_2 + '/' + accession_part_3 + '/' + accession_part_4 + '/'
    path = ftp.nlst(path)[0]
    assembly_files = [x for x in ftp.nlst(path) if x.endswith('genomic.fna.gz') and not x.endswith('_from_genomic.fna.gz')]
    assert len(assembly_files) == 1
    assembly_file = assembly_files[0]
    temp_file_name = 'temp_' + str(os.getpid()) + '.fna.gz'
    local_file_name = assembly_file.split('/')[-1]
    with open(temp_file_name, 'wb') as f:
        ftp.retrbinary('RETR ' + assembly_file, f.write)
    ftp.quit()
    return temp_file_name, local_file_name


def file_exists(assembly_accession):
    files_in_dir = [f for f in os.listdir() if os.path.isfile(f)]
    return any(assembly_accession in f for f in files_in_dir)


def run_kleborate(assembly_filename, sample_data, kleborate_dir, full_header):
    temp_kleborate_out = 'temp_kleborate_' + str(os.getpid()) + '.tsv'
    kleborate_command = ['python', os.path.join(kleborate_dir, 'Kleborate.py'),
                         '-p', kleborate_dir, '-o', temp_kleborate_out, '-r', 'on',
                         assembly_filename]
    process = subprocess.Popen(kleborate_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    _, _ = process.communicate()
    with open(temp_kleborate_out, 'rt') as f:
        kleborate_results = f.readlines()
    header = kleborate_results[0].strip().split('\t')
    sample = kleborate_results[1].strip().split('\t')
    for i, h in enumerate(header):
        if h == 'strain':
            continue
        full_header.append(h)
        sample_data.append(sample[i])
    os.remove(temp_kleborate_out)


def run_kaptive(assembly_filename, sample_data, kaptive_dir, full_header):
    temp_kaptive_out = 'temp_kaptive_' + str(os.getpid())
    kleborate_command = [os.path.join(kaptive_dir, 'kaptive.py'),
                         '-a', assembly_filename, '-o', temp_kaptive_out,
                         '-k', os.path.join(kaptive_dir, 'reference_database', 'Klebsiella_k_locus_primary_reference.gbk'),
                         '--no_seq_out', '--no_json']
    process = subprocess.Popen(kleborate_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if stderr:
        print('Kaptive error:\n' + stderr.decode())
    kaptive_table = temp_kaptive_out + '_table.txt'
    with open(kaptive_table, 'rt') as f:
        kaptive_results = f.readlines()
    header = kaptive_results[0].strip().split('\t')
    sample = kaptive_results[1].strip().split('\t')
    for i, h in enumerate(header):
        if h == 'Assembly' or 'details' in h:
            continue
        full_header.append(h)
        sample_data.append(sample[i])

    os.remove(kaptive_table)


if __name__ == '__main__':
    main()
