"""
Copyright 2025 Kat Holt, Ryan Wick
https://github.com/katholt/Kleborate/

This file is part of Kleborate. Kleborate is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Kleborate is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Kleborate. If
not, see <https://www.gnu.org/licenses/>.
"""

import gzip
import sys
import subprocess
import re
import os
import pathlib
import datetime

def load_fasta(filename):
    """
    Returns the names and sequences for the given fasta file as a list of tuples (name, seq).
    """
    fasta_seqs = []
    with get_open_func(filename)(filename, 'rt') as fasta_file:
        name = ''
        sequence = ''
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    fasta_seqs.append((name.split()[0], sequence.upper()))
                    sequence = ''
                name = line[1:]
            else:
                sequence += line
        if name:
            fasta_seqs.append((name.split()[0], sequence.upper()))
    return fasta_seqs


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)
    unknown_file = open(filename, 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == 'bz2':
        sys.exit('Error: cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('Error: cannot use zip format - use gzip instead')
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == 'gz':
        return gzip.open
    else:  # plain text
        return open


REV_COMP_DICT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
                 'D': 'H', 'H': 'D', 'N': 'N', 'r': 'y', 'y': 'r', 's': 's', 'w': 'w', 'k': 'm',
                 'm': 'k', 'b': 'v', 'v': 'b', 'd': 'h', 'h': 'd', 'n': 'n', '.': '.', '-': '-',
                 '?': '?'}


def complement_base(base):
    try:
        return REV_COMP_DICT[base]
    except KeyError:
        return 'N'


def reverse_complement(seq):
    return ''.join([complement_base(x) for x in seq][::-1])



# get tool versions for genome spec file
def get_tool_version(command):
    try:
        proc = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        output = proc.stdout + proc.stderr
        match = re.search(r'(\d+\.\d+(\.\d+)?)', output)
        return match.group(1) if match else "Unknown"
    except FileNotFoundError:
        return "Not Installed"



def get_mlst_db_version():
    """
    Returns the download date of the kpsc__mlst MLST database.
    Verifies all .fasta files and profiles.tsv exist before returning version.
    """
    import kleborate

    data_dir = pathlib.Path(kleborate.__file__).parent / 'modules' / 'kpsc__mlst' / 'data'
    if not data_dir.exists():
        return None

    db_files = sorted(data_dir.glob('*.fasta')) + [data_dir / 'profiles.tsv']

    if not all(f.exists() for f in db_files):
        return None

    return max(datetime.datetime.fromtimestamp(f.stat().st_mtime).strftime('%Y-%m-%d') for f in db_files)



def get_cgmlst_db_version():
    """
    Returns the download date of the kpsc__cgmlst database,
    """
    import kleborate

    data_dir = pathlib.Path(kleborate.__file__).parent / 'modules' / 'kpsc__cgmlst' / 'data'
    if not data_dir.exists():
        return None

    db_files = sorted(data_dir.glob('kleb_scgmlst_s-index*'))

    if not db_files:
        return None

    return max(datetime.datetime.fromtimestamp(f.stat().st_mtime).strftime('%Y-%m-%d') for f in db_files)



KLEBSIELLA_TYPING_SPEC = {
    "species": {
        "genotyping_method": "In silico species detection",
        "genotyping_schema_taxon": "Enterobacterales [NCBITaxon:91347]",
        "genotyping_database_name": "kleborate_enterobacterales__species",
        "genotyping_database_version": get_tool_version(['kleborate', '--version']),
        "genotyping_schema_name": "enterobacterales__species",
        "genotyping_software_name": "Kleborate",
        "genotyping_software_version": get_tool_version(['kleborate', '--version'])
    },  

    "ST": {
        "genotyping_method": "MLST",
        "genotyping_schema_taxon": "Klebsiella pneumoniae species complex [NCBITaxon:3390273]",
        "genotyping_database_name": "pubmlst_klebsiella_seqdef",
        "genotyping_database_version": get_mlst_db_version(), 
        "genotyping_database_version": "2024-12-31",
        "genotyping_schema_name": "MLST",
        "genotyping_software_name": "Kleborate",
        "genotyping_software_version": get_tool_version(['kleborate', '--version'])
    },

    "subspecies": {
        "genotyping_method": "In silico subspecies detection",
        "genotyping_schema_taxon": "Klebsiella pneumoniae species complex [NCBITaxon:3390273]",
        "genotyping_database_name": "pubmlst_klebsiella_seqdef",
        "genotyping_database_version": get_mlst_db_version(), 
        "genotyping_database_version": "2024-12-31",
        "genotyping_schema_name": "MLST",
        "genotyping_software_name": "Kleborate",
        "genotyping_software_version": get_tool_version(['kleborate', '--version']),
        "genotyping_details": "subspecies call based on MLST"
    },

    "K_locus": {
        "genotyping_method": "In silico serotyping",
        "genotyping_schema_taxon": "Klebsiella pneumoniae species complex [NCBITaxon:3390273]",
        "genotyping_database_name": "kaptive_Klebsiella_k_locus_primary_reference", 
        "genotyping_database_version": get_tool_version(['kaptive', '--version']),
        "genotyping_schema_name": "Klebsiella_k_locus_primary_reference",
        "genotyping_software_name": "Kaptive",
        "genotyping_software_version": get_tool_version(['kaptive', '--version']),
        "genotype_predicted_phenotype": "K_type"
    },

    "O_locus": {
        "genotyping_method": "In silico serotyping",
        "genotyping_schema_taxon": "Klebsiella pneumoniae species complex [NCBITaxon:3390273]",
        "genotyping_database_name": "kaptive_Klebsiella_o_locus_primary_reference", 
        "genotyping_database_version": get_tool_version(['kaptive', '--version']),
        "genotyping_schema_name": "Klebsiella_o_locus_primary_reference",
        "genotyping_software_name": "Kaptive",
        "genotyping_software_version": get_tool_version(['kaptive', '--version']),
        "genotype_predicted_phenotype": "O_type"
    },

    "cgST": {
        "genotyping_method": "cgMLST",
        "genotyping_schema_taxon": "Klebsiella pneumoniae species complex [NCBITaxon:3390273]",
        "genotyping_database_name": "pubmlst_klebsiella_seqdef",
        "genotyping_database_version": get_cgmlst_db_version(),
        "genotyping_database_version": "2024-12-31",
        "genotyping_schema_name": "scgMLST629_S",
        "genotyping_software_name": "MiST",
        "genotyping_software_version": get_tool_version(['mist', '--version'])
    },

    "LIN code": {
        "genotyping_method": "LIN code",
        "genotyping_schema_taxon": "Klebsiella pneumoniae species complex [NCBITaxon:3390273]",
        "genotyping_database_name": "pubmlst_klebsiella_seqdef",
        "genotyping_database_version": get_cgmlst_db_version(),
        "genotyping_database_version": "2024-12-31",
        "genotyping_schema_name": "scgMLST629_S",
        "genotyping_software_name": "MiST",
        "genotyping_software_version": get_tool_version(['mist', '--version'])
    },

    "Sublineage": {
        "genotyping_method": "Sublineage",
        "genotyping_schema_taxon": "Klebsiella pneumoniae species complex [NCBITaxon:3390273]",
        "genotyping_database_name": "pubmlst_klebsiella_seqdef",
        "genotyping_database_version": get_cgmlst_db_version(),
        "genotyping_database_version": "2024-12-31",
        "genotyping_schema_name": "scgMLST629_S",
        "genotyping_software_name": "MiST",
        "genotyping_software_version": get_tool_version(['mist', '--version'])
    },

    "Clonal group": {
        "genotyping_method": "Clonal group",
        "genotyping_schema_taxon": "Klebsiella pneumoniae species complex [NCBITaxon:3390273]",
        "genotyping_database_name": "pubmlst_klebsiella_seqdef",
        "genotyping_database_version": get_cgmlst_db_version(),
        "genotyping_database_version": "2024-12-31",
        "genotyping_schema_name": "scgMLST629_S",
        "genotyping_software_name": "MiST",
        "genotyping_software_version": get_tool_version(['mist', '--version'])
    }
}



# HaRmronization headers
res_headers = [
    'AGly_acquired', 'Col_acquired', 'Fcyn_acquired', 'Flq_acquired', 'Gly_acquired', 
    'MLS_acquired', 'Phe_acquired', 'Rif_acquired', 'Sul_acquired', 'Tet_acquired', 
    'Tgc_acquired', 'Tmt_acquired', 'Bla_acquired', 'Bla_inhR_acquired', 'Bla_ESBL_acquired', 
    'Bla_ESBL_inhR_acquired', 'Bla_Carb_acquired', 'Bla_chr', 'SHV_mutations', 
    'Omp_mutations', 'Col_mutations', 'Flq_mutations', 'truncated_resistance_hits', 
    'spurious_resistance_hits'
]


annotation_fields = [
        'Genetic Variation Type','Drug Class','Input Sequence ID','Input Gene Length', 'Input Gene Start', 'Input Gene Stop', 'Reference Gene Length',
        'Reference Gene Start', 'Reference Gene Stop', 'Sequence Identity', 'Coverage','Reference Accession','Strand Orientation',
        'Analysis Software Name', 'Analysis Software Version', 'Reference Database Name',
        'Reference Database Version','Input Protein Length','Reference Protein Length','Input Protein Start', 'Input Protein Stop',
        'Reference Protein Start', 'Reference Protein Stop'
]
