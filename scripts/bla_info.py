#!/usr/bin/env python3

"""
This script generates a table of beta-lactamase class information.

It takes two arguments:
  1) The Lahey info file at ftp://ftp.ncbi.nlm.nih.gov/pathogen/betalactamases/Lahey.tab
  2) A GenBank file downloaded from these records:
     https://www.ncbi.nlm.nih.gov/nuccore?term=313047%5BBioProject%5D
     (use 'Send to' then 'File' and set the format to GenBank)

Usage:
  ./bla_info.py Lahey.tab sequence.gb > bla_info_table

The output has three columns:
  * allele name
  * allele description
  * bla class

The bla class column can have these six possible values:
  * Bla_broad (broad spectrum beta-lactamases)
  * Bla_broad_inhR (broad spectrum beta-lactamases with resistance to beta-lactamase inhibitors)
  * Bla_ESBL (extended spectrum beta-lactamases)
  * Bla_ESBL_inhR (extended spectrum beta-lactamases with resistance to beta-lactamase inhibitors)
  * Bla_Carb (carbapenemase)
  * Bla (beta-lactamase that is none of the above)
"""

import sys
import re
from Bio import SeqIO


def main():
    lahey_filename = sys.argv[1]
    genbank_filename = sys.argv[2]

    lahey_descriptions = get_lahey_descriptions(lahey_filename)
    genbank_descriptions = get_genbank_descriptions(genbank_filename)

    all_alleles = sorted(set(list(lahey_descriptions.keys()) + list(genbank_descriptions.keys())),
                         key=allele_name_sorting_key)

    for allele in all_alleles:

        # If the allele is in only one of the two sources, then we just get the description/class
        # from the single source.
        if allele in lahey_descriptions and allele not in genbank_descriptions:
            description = lahey_descriptions[allele]
            bla_class = bla_class_from_description(description)
        elif allele in genbank_descriptions and allele not in lahey_descriptions:
            description = genbank_descriptions[allele]
            bla_class = bla_class_from_description(description)

        # If the allele is in both Lahey and the GenBank...
        else:
            lahey_description = lahey_descriptions[allele]
            genbank_description = genbank_descriptions[allele]
            lahey_bla_class = bla_class_from_description(lahey_description)
            genbank_bla_class = bla_class_from_description(genbank_description)

            # If they agree on the class, then we use the longer description (assuming that the
            # longer one is more descriptive).
            if lahey_bla_class == genbank_bla_class:
                bla_class = lahey_bla_class
                description = lahey_description
                if len(genbank_description) > len(lahey_description):
                    description = genbank_description

            # If they disagree on the class, then we choose whichever class is 'worse'.
            else:
                class_values  = {'Bla': 0,
                                 'Bla_broad': 1, 'Bla_broad_inhR': 2,
                                 'Bla_ESBL': 3, 'Bla_ESBL_inhR': 4,
                                 'Bla_Carb': 5}
                if class_values[lahey_bla_class] > class_values[genbank_bla_class]:
                    bla_class = lahey_bla_class
                    description = lahey_description
                else:
                    bla_class = genbank_bla_class
                    description = genbank_description

        print('\t'.join([allele, description, bla_class]))


def get_lahey_descriptions(lahey_filename):
    descriptions = {}
    with open(lahey_filename, 'rt') as lahey_file:
        for line in lahey_file:
            if line.startswith('#'):
                continue
            line_parts = line.strip().split('\t')
            allele = line_parts[0]
            descriptions[allele] = line_parts[7].rsplit(' ', 1)[0]
    return descriptions


def get_genbank_descriptions(genbank_filename):
    descriptions = {}
    for seq_record in SeqIO.parse(sys.argv[2], 'genbank'):
        for feature in seq_record.features:
            if feature.type == 'CDS':
                allele = feature.qualifiers.get('allele', ['none'])[0]
                gene = feature.qualifiers.get('gene', ['none'])[0]
                product = feature.qualifiers.get('product', [''])[0]
                note = feature.qualifiers.get('note', [''])[0]
                function = feature.qualifiers.get('function', [''])[0]

                if allele.lower().startswith('bla') and gene.lower().startswith('bla'):
                    allele = allele[3:]  # remove 'bla'
                    if allele not in descriptions:  # Don't overwrite the Lahey description.

                        # The product likely ends with the allele name. If so, remove it.
                        allele_name_len = len(allele)
                        if allele != 'none' and allele_name_len and \
                                product.lower().endswith(allele.lower()):
                            product = product[:-allele_name_len]
                        product = product.rstrip()
                        if product.lower().endswith(' bla'):
                            product = product[:-4]

                        # The 'product' feature probably has all of the description we need.
                        # However, if there is additional information in the 'note' or 'function'
                        # features, then add those too.
                        description = product
                        if ('carbapenem' in note and 'carbapenem' not in description) or \
                                ('extended' in note and 'extended' not in description) or \
                                ('broad' in note and 'broad' not in description) or \
                                ('inhibitor' in note and 'inhibitor' not in description):
                            description += ' ' + note
                        if ('carbapenem' in function and 'carbapenem' not in description) or \
                                ('extended' in function and 'extended' not in description) or \
                                ('broad' in function and 'broad' not in description) or \
                                ('inhibitor' in function and 'inhibitor' not in description):
                            description += ' ' + function
                        description = description.replace(',', ' ')
                        description = description.replace(';', ' ')
                        description = re.sub('\s+', ' ', description).strip()

                        # Special case to fix a particular long description.
                        if description == 'subclass B1 metallo-beta-lactamase New Delhi ' + \
                                          'metallo-beta-lactamase+ New Delhi MBL carbapenemase':
                            description = 'subclass B1 metallo-beta-lactamase carbapenemase'

                        descriptions[allele] = description
    return descriptions


def bla_class_from_description(description):
    if 'carbapenem-hydrolyzing' in description or \
            'carbapenemase' in description or \
            'confers carbapenem resistance' in description or \
            'hydrolyses beta-lactams including carbapenems' in description:
        return 'Bla_Carb'
    elif 'extended-spectrum' in description or \
            'extended spectrum' in description or \
            'ESBL' in description:
        if 'inhibitor-resistant' in description:
            return 'Bla_ESBL_inhR'
        else:
            return 'Bla_ESBL'
    elif 'broad-spectrum' in description or \
            'broad spectrum' in description:
        if 'inhibitor-resistant' in description:
            return 'Bla_broad_inhR'
        else:
            return 'Bla_broad'
    return 'Bla'


def allele_name_sorting_key(allele_name):
    """
    This function defines a sorting key for allele names. It puts padding zeros in front of
    numbers so 'TEM-2' will come before 'TEM-10', etc.
    """
    new_name_parts = []
    for name_part in allele_name.split('-'):
        try:
            name_part_int = int(name_part)
            new_name_parts.append('%05d' % name_part_int)
        except ValueError:
            new_name_parts.append(name_part)
    return '-'.join(new_name_parts)


if __name__ == '__main__':
    main()
