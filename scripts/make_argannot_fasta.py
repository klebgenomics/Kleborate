#!/usr/bin/env python3

"""
This script modifies the ARGannot_r2.fasta file (comes with SRST2) for use with Kleborate.
Specifically, it removes alleles which are not in the ARGannot_clustered80_r2.csv file.

It takes two arguments:
  1) SRST2's copy of ARGannot_r2.fasta
  2) the ARGannot_clustered80_r2.csv file made by make_argannot_csv.py

Usage:
  ./make_argannot_fasta.py path/to/srst2/data/ARGannot_r2.fasta ARGannot_clustered80_r2.csv > ARGannot_r2.fasta
"""

import sys


def main():
    argannot_fasta_filename = sys.argv[1]
    argannot_csv_filename = sys.argv[2]

    allele_names = set()
    with open(argannot_csv_filename, 'rt') as argannot_csv_file:
        for line in argannot_csv_file:
            allele_names.add(line.split(',')[4])

    fasta_seqs = load_fasta(argannot_fasta_filename)
    fasta_seqs = sorted(fasta_seqs, key=fasta_sorting_key)

    for name, seq in fasta_seqs:
        name = name.replace('NimB_Nitroimidazole_Gene', 'NimB')  # fix verbose gene name
        allele_name = name.split('__')[2]
        if allele_name in allele_names:
            print('>' + name)
            print(seq, end='')


def load_fasta(filename):
    """Returns the names and sequences for the given fasta file."""
    fasta_seqs = []
    with open(filename, 'rt') as fasta_file:
        name = ''
        sequence = ''
        for line in fasta_file:
            stripped_line = line.strip()
            if not stripped_line:
                continue
            if stripped_line[0] == '>':  # Header line = start of new contig
                if name:
                    fasta_seqs.append((name, sequence))
                    sequence = ''
                name = stripped_line[1:]
            else:
                sequence += line
        if name:
            fasta_seqs.append((name, sequence))
    return fasta_seqs


def fasta_sorting_key(fasta_record):
    """
    This function defines a sorting key for fasta so 'TEM-2' will come before 'TEM-10', etc.
    """
    name_parts = fasta_record[0].split()[0].split('__')
    return int(name_parts[0]), int(name_parts[3])


if __name__ == '__main__':
    main()
