"""
Copyright 2017 Kat Holt
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/katholt/Kleborate/

This file is part of Kleborate. Kleborate is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Kleborate is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Kleborate. If
not, see <http://www.gnu.org/licenses/>.
"""

# reports best match
# * : best matching allele is not precise match
# -nLV : best matching ST is n-locus variant
# if an annotation column is provided (such as clonal complex) in the final column of the profiles file,
#	  this annotation will be reported in column 2 of the output table
# NOTE there is a bug with the culling_limit parameter in older versions of BLAST+...
# This code has been tested with BLAST+2.2.30. It does not work with BLAST2.2.25. Not sure about other versions.

import collections
import os
import sys
import subprocess
from optparse import OptionParser


def main():
	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	# options
	parser.add_option("-s", "--seqs", action="store", dest="seqs", default="",
					  help="rmpA and rmpA2 allele sequences file")
	parser.add_option("-d", "--rmpA_db", action="store", dest="database", default="",
					  help="rmpA profile database (col1 = allele, col2 = lineage)")
	parser.add_option("-m", "--minident", action="store", dest="minident", default="95",
					  help="Minimum percent identity (default 95)")
	return parser.parse_args()

if __name__ == "__main__":
	(options, args) = main()

	if options.database == "":
		sys.exit("No rmpA profiles databse provided (-d)")
	if options.seqs == "":
		sys.exit("No rmpA allele sequences provided (-s)")
	else:
		(path, fileName) = os.path.split(options.seqs)
		if not os.path.exists(options.seqs + ".nin"):
			with open(os.devnull, 'w') as devnull:
				subprocess.check_call("makeblastdb -dbtype nucl -in " + options.seqs,
									  stdout=devnull, shell=True)

	# read in rmpA database
	st_info = {}  # key = st, value = info relating to this ST, eg clonal group
	header = []
	info_title = "info"
	with open(options.database, "r") as f:
		for line in f:
			fields = line.rstrip().split("\t")
			if len(header)==0:
				header = fields
			else:
				st_info[fields[0]]=fields[1]

	# print header
	print("\t".join(["strain","rmpA_allele","rmpA_lineage","rmpA2_allele"]))

	# search input assemblies
	for contigs in args:
		(_, fileName) = os.path.split(contigs)
		(name, ext) = os.path.splitext(fileName)

		# blast against all rmpA and rmpA2 alleles
		f = os.popen("blastn -task blastn -db " + options.seqs + " -query " + contigs +
					 " -outfmt '6 sacc pident slen length score' -dust no -evalue 1E-20 -word_size 32"
					 " -max_target_seqs 10000 -culling_limit 1 -perc_identity " + options.minident)
					 
		rmpA_calls = []
		rmpA2_calls = []
		for line in f:
			fields = line.rstrip().split("\t")
			(gene_id, pcid, allele_length, length, score) = (fields[0], float(fields[1]), int(fields[2]),
															 int(fields[3]), float(fields[4]))
			
			if length > (allele_length/2):
				gene = gene_id.split("_")[0]
				if gene=="rmpA":
					info = " ("+st_info[gene_id.split("_")[1]]+")" # predict from best hit
					if pcid < 100.00 or allele_length < length:
						gene_id += "*" #indicate imprecise hit
					rmpA_calls.append(gene_id+info)
				else:
					if pcid < 100.00 or allele_length < length:
						rmpA2_calls.append(gene_id + "*") #indicate imprecise hit
					else:
						rmpA2_calls.append(gene_id)					
					
		f.close()
		
		if len(rmpA_calls) == 0:
			rmpA_calls.append("-")
		if len(rmpA2_calls) == 0:
			rmpA2_calls.append("-")
		
		print ("\t".join([name,",".join(rmpA_calls),",".join(rmpA2_calls)]))
