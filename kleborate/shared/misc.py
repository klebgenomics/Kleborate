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


KLEBSIELLA_TYPING_SPEC = {
    "species": {
        "genotyping_method": "In silico species detection",
        "genotyping_schema_taxon": "Enterobacterales",
        "genotyping_database_name": "kleborate_enterobacterales__species",
        "genotyping_database_version": get_tool_version(['kleborate', '--version']),
        "genotyping_schema_name": "enterobacterales__species",
        "genotyping_software_name": "Kleborate",
        "genotyping_software_version": get_tool_version(['kleborate', '--version'])
    },  

    "ST": {
        "genotyping_method": "MLST",
        "genotyping_schema_taxon": "Klebsiella pneumoniae species complex[NCBITaxon:3390273]",
        "genotyping_database_name": "pubmlst_klebsiella_seqdef",
        "genotyping_database_version": "2024-12-31", 
        "genotyping_schema_name": "MLST",
        "genotyping_software_name": "Kleborate",
        "genotyping_software_version": get_tool_version(['kleborate', '--version'])
    },

    "K_locus": {
        "genotyping_method": "in silico serotyping",
        "genotyping_schema_taxon": "Klebsiella pneumoniae species complex[NCBITaxon:3390273]",
        "genotyping_database_name": "kaptive_Klebsiella_k_locus_primary_reference", 
        "genotyping_database_version": get_tool_version(['kaptive', '--version']),
        "genotyping_schema_name": "Klebsiella_k_locus_primary_reference",
        "genotyping_software_name": "Kaptive",
        "genotyping_software_version": get_tool_version(['kaptive', '--version']),
        "genotype_predicted_phenotype": "K_type"
    },

    "O_locus": {
        "genotyping_method": "in silico serotyping",
        "genotyping_schema_taxon": "Klebsiella pneumoniae species complex[NCBITaxon:3390273]",
        "genotyping_database_name": "kaptive_Klebsiella_o_locus_primary_reference", 
        "genotyping_database_version": get_tool_version(['kaptive', '--version']),
        "genotyping_schema_name": "Klebsiella_o_locus_primary_reference",
        "genotyping_software_name": "Kaptive",
        "genotyping_software_version": get_tool_version(['kaptive', '--version']),
        "genotype_predicted_phenotype": "O_type"
    },

    "cgST": {
        "genotyping_method": "cgMLST",
        "genotyping_schema_taxon": "Klebsiella pneumoniae species complex[NCBITaxon:3390273]",
        "genotyping_database_name": "pubmlst_klebsiella_seqdef",
        "genotyping_database_version": get_tool_version(['mist', '--version']),
        "genotyping_schema_name": "scgMLST629_S",
        "genotyping_software_name": "MiST",
        "genotyping_software_version": get_tool_version(['mist', '--version'])
    },

    "LINcodes": {
        "genotyping_method": "LINcode",
        "genotyping_schema_taxon": "Klebsiella pneumoniae species complex[NCBITaxon:3390273]",
        "genotyping_database_name": "pubmlst_klebsiella_seqdef",
        "genotyping_database_version": get_tool_version(['mist', '--version']),
        "genotyping_schema_name": "scgMLST629_S",
        "genotyping_software_name": "MiST",
        "genotyping_software_version": get_tool_version(['mist', '--version'])
    },

    "Sublineage": {
        "genotyping_method": "Sublineage",
        "genotyping_schema_taxon": "Klebsiella pneumoniae species complex[NCBITaxon:3390273]",
        "genotyping_database_name": "pubmlst_klebsiella_seqdef",
        "genotyping_database_version": get_tool_version(['mist', '--version']),
        "genotyping_schema_name": "scgMLST629_S",
        "genotyping_software_name": "MiST",
        "genotyping_software_version": get_tool_version(['mist', '--version'])
    },

    "Clonal group": {
        "genotyping_method": "Clonal group",
        "genotyping_schema_taxon": "Klebsiella pneumoniae species complex[NCBITaxon:3390273]",
        "genotyping_database_name": "pubmlst_klebsiella_seqdef",
        "genotyping_database_version": get_tool_version(['mist', '--version']),
        "genotyping_schema_name": "scgMLST629_S",
        "genotyping_software_name": "MiST",
        "genotyping_software_version": get_tool_version(['mist', '--version'])
    }
}
# KLEBSIELLA_TYPING_SPEC = {
#     "species": {
#         "genotyping_method": "In silico species detection",
#         "genotyping_schema_taxon": "Enterobacterales",
#         "genotyping_database_name": "kleborate_enterobacterales__species",
#         "genotyping_database_version": get_tool_version(['kleborate', '--version']),
#         "genotyping_schema_name": "enterobacterales__species",
#         "genotyping_software_name": "Kleborate",
#         "genotyping_software_version": get_tool_version(['kleborate', '--version'])
#     },  

#     "ST": {
#         "genotyping_method": "MLST",
#         "genotyping_schema_taxon": "Klebsiella pneumoniae species complex[NCBITaxon:3390273]",
#         "genotyping_database_name": "pubmlst_klebsiella_seqdef",
#         "genotyping_database_version": "2024-12-31", 
#         "genotyping_schema_name": "MLST",
#         "genotyping_software_name": "Kleborate",
#         "genotyping_software_version": get_tool_version(['kleborate', '--version'])
#     },

#     "K_locus": {
#         "genotyping_method": "in silico serotyping",
#         "genotyping_schema_taxon": "Klebsiella pneumoniae species complex[NCBITaxon:3390273]",
#         "genotyping_database_name": "kaptive_Klebsiella_k_locus_primary_reference", 
#         "genotyping_database_version": get_tool_version(['kaptive', '--version']),
#         "genotyping_schema_name": "Klebsiella_k_locus_primary_reference",
#         "genotyping_software_name": "Kaptive",
#         "genotyping_software_version": get_tool_version(['kaptive', '--version']),
#         "genotype_predicted_phenotype": "K_type"
#     },

#     "O_locus": {
#         "genotyping_method": "in silico serotyping",
#         "genotyping_schema_taxon": "Klebsiella pneumoniae species complex[NCBITaxon:3390273]",
#         "genotyping_database_name": "kaptive_Klebsiella_o_locus_primary_reference", 
#         "genotyping_database_version": get_tool_version(['kaptive', '--version']),
#         "genotyping_schema_name": "Klebsiella_o_locus_primary_reference",
#         "genotyping_software_name": "Kaptive",
#         "genotyping_software_version": get_tool_version(['kaptive', '--version']),
#         "genotype_predicted_phenotype": "O_type"
#     },

#     "cgST": {
#         "genotyping_method": "cgMLST",
#         "genotyping_schema_taxon": "Klebsiella pneumoniae species complex[NCBITaxon:3390273]",
#         "genotyping_database_name": "pubmlst_klebsiella_seqdef",
#         "genotyping_database_version": get_tool_version(['mist', '--version']),
#         "genotyping_schema_name": "scgMLST629_S",
#         "genotyping_software_name": "MiST",
#         "genotyping_software_version": get_tool_version(['mist', '--version'])
#     },

#     "LINcodes": {
#         "genotyping_method": "LINcode",
#         "genotyping_schema_taxon": "Klebsiella pneumoniae species complex[NCBITaxon:3390273]",
#         "genotyping_database_name": "pubmlst_klebsiella_seqdef",
#         "genotyping_database_version": get_tool_version(['mist', '--version']),
#         "genotyping_schema_name": "scgMLST629_S",
#         "genotyping_software_name": "MiST",
#         "genotyping_software_version": get_tool_version(['mist', '--version'])
#     }
# }



def get_db_version(db_path=None, module_name=None):
    """
    Returns the version of a database. 
    """
    # pre-installed versions of the databases (the "defaults")
    DEFAULTS = {
        'kpsc__mlst': '2024-12-31',
        'kpsc__cgmlst': '2024-12-31'
    }

    # if the user has provided a custom DB path, get the download date
    if db_path and os.path.exists(db_path):
        timestamp = os.path.getmtime(db_path)
        return datetime.datetime.fromtimestamp(timestamp).strftime('%Y-%m-%d')

    return DEFAULTS.get(module_name)

def get_all_db_versions(custom_paths=None):
    """
    Iterates through all modules defined in get_presets and assigns a version.
    custom_paths: a dict mapping module name to file path, e.g. {'kpsc__mlst': '/path/to/db'}
    """
    presets = get_presets()
    version_tracking = {}
    paths = custom_paths or {}

    for species, groups in presets.items():
        all_modules = groups['check'] + groups['pass']
        
        for item in all_modules:
            module_id = item[0] if isinstance(item, tuple) else item
            
            # Record the version
            version_tracking[module_id] = get_db_version(
                db_path=paths.get(module_id), 
                module_name=module_id
            )
            
    return version_tracking


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
        'Genetic_variation_type','Drug_class','Input_sequence_ID','Input_gene_length', 'Input_gene_start', 'Input_gene_stop', 'Reference_gene_length',
        'Reference_gene_start', 'Reference_gene_stop', 'Sequence_identity', 'Coverage','Reference_accession','Strand_orientation',
        'Software_name', 'Software_version', 'Reference_database_name',
        'Reference_database_version','Input_protein_length','Reference_protein_length','Input_protein_start', 'Input_protein_stop','Antimicrobial_agent', 
        'Coverage_depth', 'Coverage_ratio','Predicted_phenotype','predicted_phenotype_confidence_level', 
        'Reference_protein_start', 'Reference_protein_stop','Resistance_mechanism'
]
