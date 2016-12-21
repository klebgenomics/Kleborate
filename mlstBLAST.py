# reports best match
# * : best matching allele is not precise match
# -nLV : best matching ST is n-locus variant
# if an annotation column is provided (such as clonal complex) in the final column of the profiles file,
#	 this annotation will be reported in column 2 of the output table
# NOTE there is a bug with the culling_limit parameter in older versions of BLAST+...
# This code has been tested with BLAST+2.2.30. It does not work with BLAST2.2.25. Not sure about other versions.

import string, re, collections
import os, sys, subprocess
from optparse import OptionParser
	
def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	# options
	parser.add_option("-s", "--seqs", action="store", dest="seqs", help="MLST allele sequences file", default="")
	parser.add_option("-d", "--database", action="store", dest="database", help="MLST profile database (col1=ST, other cols=loci, must have loci names in header)", default="")
	parser.add_option("-i", "--info", action="store", dest="info", help="Info (clonal group, lineage, etc) provided in las column of profiles (yes (default), no)", default="yes")
	parser.add_option("-m", "--minident", action="store", dest="minident", help="Minimum percent identity (default 95)", default="95")
	parser.add_option("-n", "--maxmissing", action="store", dest="maxmissing", help="Maximum missing/uncalled loci to still calculate closest ST (default 3)", default="3")
	return parser.parse_args()

if __name__ == "__main__":

	(options, args) = main()

	if options.database=="":
		DoError("No MLST profiles databse provided (-d)")
		
	if options.seqs=="":
		DoError("No MLST allele sequences provided (-s)")
	else:
		(path,fileName) = os.path.split(options.seqs)
		if not os.path.exists(options.seqs + ".nin"):
			os.system("makeblastdb -dbtype nucl -logfile blast.log -in " + options.seqs)
		(fileName,ext) = os.path.splitext(fileName)
		
	def getClosestLocusVariant(query, sts):
		closest = []
		min_dist = len(query)
		for index, item in enumerate(query):
			if item == "-":	
				query[index] = "0"
		for st in sts:
			d = sum(map(lambda x,y: bool(int(x)-int(y)),st.split(","),query))
			if d <= min_dist:
				closest.append(int(sts[st]))
				min_dist = d
		return (str(min(closest)), min_dist)
		
	sts = {} # key = concatenated string of alleles, value = st
	st_info = {} # key = st, value = info relating to this ST, eg clonal group
	max_st = 0 # changeable variable holding the highest current ST, incremented when novel combinations are encountered
	header = []
	info_title = "info"
	f = file(options.database,"r")
	for line in f:
		fields = line.rstrip().split("\t")
		if len(header)==0:
			header = fields
			header.pop(0) # remove st label
			if options.info=="yes":
				info_title = header.pop() # remove info label
		else:
			st = fields.pop(0)
			if options.info=="yes":
				info = fields.pop()
			sts[",".join(fields)] = st
			if int(st) > max_st:
				max_st = int(st)
			if options.info=="yes":
				st_info[st] = info
	f.close()
	
	best_match = collections.defaultdict(dict) # key1 = strain, key2 = locus, value = best match (clean for ST, annotated)
	perfect_match = collections.defaultdict(dict) # key1 = strain, key2 = locus, value = perfect match if available
	
	# print header
	if options.info=="yes":
		print "\t".join(["strain",info_title,"ST"]+header)
	else:
		print "\t".join(["strain","ST"]+header)
	
	for contigs in args:
		(dir,fileName) = os.path.split(contigs)
		(name,ext) = os.path.splitext(fileName)

		best_score = {} # key = locus, value = BLAST score for best matching allele encountered so far
		best_allele = {} # key = locus, value = best allele (* if imprecise match)

		# blast against all
		f = os.popen("blastn -task blastn -db " + options.seqs + " -query " + contigs + " -outfmt '6 sacc pident slen length score' -ungapped -dust no -evalue 1E-20 -word_size 32 -max_target_seqs 10000 -culling_limit 2 -perc_identity " + options.minident) 

		for line in f:
			fields = line.rstrip().split("\t")
			(gene_id,pcid,length,allele_length,score) = (fields[0],float(fields[1]),int(fields[2]),int(fields[3]),float(fields[4]))
			#print line
			if "__" in gene_id:
				# srst2 formated file
				gene_id_components = gene_id.split("__")
				locus = gene_id_components[1]
				allele = gene_id_components[2]
			else:
				allele = gene_id
				locus = gene_id.split("_")[0]
			if pcid < 100.00 or allele_length < length:
				allele = allele + "*" # imprecise match
			# store best match for each one locus
			if locus in best_score:
				if score > best_score[locus]:
					# update
					best_score[locus] = score
					best_allele[locus] = allele.split("_")[1] # store number only
			else:
				# initialise
				best_score[locus] = score
				best_allele[locus] = allele.split("_")[1] # store number only
		f.close()	

		best_st = []
		best_st_annotated = []
		
		for locus in header:
			if locus in best_allele:
				allele = best_allele[locus]
				allele_number = allele.replace("*","")
				best_st.append(allele_number)
				best_st_annotated.append(allele) # will still have * if imperfect match
			else:
				best_st.append("-")
				best_st_annotated.append("-")
			
		# assign ST
		bst = ",".join(best_st)
		
		mismatch_loci = 0
		
		if bst in sts:
			bst = sts[bst]
		elif bst.count("-") <= int(options.maxmissing):
			(bst, mismatch_loci) = getClosestLocusVariant(best_st, sts)
		else:
			bst = "0"

		# pull info column
		if options.info=="yes":
			info_final = ""
			if bst in st_info:
				info_final = st_info[bst]
		
		if best_st != best_st_annotated:
			bst += "*"
		if mismatch_loci > 0:
			bst += "-" + str(mismatch_loci) + "LV"
		
		if options.info=="yes":
			print "\t".join([name,info_final,bst] + best_st_annotated)
		else:
			print "\t".join([name,bst] + best_st_annotated)