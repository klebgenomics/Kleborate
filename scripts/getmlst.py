#!/usr/bin/env python

'''
Download MLST datasets from this site: http://pubmlst.org/data/ by
parsing an xml file (http://pubmlst.org/data/dbases.xml).

Data is downloaded for a species determined by the user:
- profiles (maps STs to allele numbers)
- numbered sequences for each locus in the scheme

In addition, the alleles are concatenated together for use with Kleborate or SRST2.

If the species name input by the user matches multiple <species> in the
xml file, the script simply reports the possible matches so the user can
try again.
'''

from argparse import ArgumentParser
import xml.dom.minidom as xml
#import urllib2 as url
from urllib.request import urlopen
import re, os, glob
#from urlparse import urlparse
from urllib.parse import urlparse

def parse_args():
	parser = ArgumentParser(description='Download MLST datasets by species'
										'from pubmlst.org.')

	parser.add_argument('--repository_url',
						metavar = 'URL',
						default = 'http://pubmlst.org/data/dbases.xml',
						help = 'URL for MLST repository XML index')

	parser.add_argument('--species',
						metavar = 'NAME',
						required = True,
						help = 'The name of the species that you want to download (e.g. "Klebsiella pneumoniae")')
						
	parser.add_argument('--force_scheme_name',
						action="store_true",
						default = False,
						help = 'Flage to force downloading of specific scheme name (e.g. "Clostridium difficile")')
						
	return parser.parse_args()

# test if a node is an Element and that it has a specific tag name
def testElementTag(node, name):
		return node.nodeType == node.ELEMENT_NODE and node.localName == name

# Get the text from an element node with a text node child
def getText(element):
	result = ''
	for node in element.childNodes:
		if node.nodeType == node.TEXT_NODE:
			result += node.data
	return normaliseText(result)

# remove unwanted whitespace including linebreaks etc.
def normaliseText(s):
	return ' '.join(s.split())

# A collection of interesting information about a taxa
class SpeciesInfo(object):
	def __init__(self):
		self.name = None # String name of species
		self.database_url = None # URL as string
		self.retrieved = None # date as string 
		self.profiles_url = None # URL as string 
		self.profiles_count = None # positive integer
		self.loci = [] # list of loci 

class LocusInfo(object):
	def __init__(self):
		self.url = None
		self.name = None

# retrieve the interesting information for a given sample element
def getSpeciesInfo(species_node, species, exact):
	this_name = getText(species_node)
	store = False
	if exact:
		if this_name == species:
			store = True
	else:
		if this_name.startswith(species):
			store = True
	if store:
		info = SpeciesInfo()
		info.name = this_name
		for mlst_node in species_node.getElementsByTagName('mlst'):
			for database_node in mlst_node.getElementsByTagName('database'):
				for database_child_node in database_node.childNodes:
					if testElementTag(database_child_node, 'url'):
						info.database_url = getText(database_child_node)
					elif testElementTag(database_child_node, 'retrieved'):
						info.retrieved = getText(database_child_node)
					elif testElementTag(database_child_node, 'profiles'):
						for profile_count in database_child_node.getElementsByTagName('count'):
							info.profiles_count = getText(profile_count)
						for profile_url in database_child_node.getElementsByTagName('url'):
							info.profiles_url = getText(profile_url)
					elif testElementTag(database_child_node, 'loci'):
						for locus_node in database_child_node.getElementsByTagName('locus'):
							locus_info = LocusInfo()
							locus_info.name = getText(locus_node)
							for locus_url in locus_node.getElementsByTagName('url'):
								locus_info.url = getText(locus_url)
							info.loci.append(locus_info)
		return info
	else:
		return None


def main():
	args = parse_args()
	docFile = urlopen(args.repository_url)
	doc = xml.parse(docFile)
	root = doc.childNodes[0]
	found_species = []
	for species_node in root.getElementsByTagName('species'):
		info = getSpeciesInfo(species_node, args.species, args.force_scheme_name)
		if info != None:
			found_species.append(info)
	if len(found_species) == 0:
		print ("No species matched your query.")
		exit(1)
	if len(found_species) > 1:
		print ("The following {} species match your query, please be more specific:".format(len(found_species)))
		for info in found_species:
			print (info.name)
		exit(2)

	assert len(found_species) == 1
	species_info = found_species[0]
	species_name_underscores = species_info.name.replace(' ', '_')
	species_name_underscores = species_name_underscores.replace('/', '_')

	# Remove any Bowtie/Samtools index files that already exist.
	for filename in glob.glob(species_name_underscores + '*.bt2'):
		os.remove(filename)
	for filename in glob.glob(species_name_underscores + '*.fai'):
		os.remove(filename)

	# output information for the single matching species
	species_all_fasta_filename = species_name_underscores + '.fasta'
	species_all_fasta_file = open(species_all_fasta_filename, 'w')
	profile_path = urlparse(species_info.profiles_url).path
	profile_filename = profile_path.split('/')[-1]
	profile_doc = urlopen(species_info.profiles_url)
	profile_doc_string = profile_doc.read().decode('utf-8')
	profile_file = open(profile_filename, 'w')
	profile_file.write(profile_doc_string)
	profile_file.close()
	profile_doc.close()
	for locus in species_info.loci:
		locus_path = urlparse(locus.url).path
		locus_filename = locus_path.split('/')[-1]
		locus_doc = urlopen(locus.url)
		locus_fasta_content = locus_doc.read().decode('utf-8')
		species_all_fasta_file.write(locus_fasta_content)
		locus_doc.close()
	species_all_fasta_file.close()


if __name__ == '__main__':
	main()
