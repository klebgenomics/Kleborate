# add reporting of clonal group
import string, re, collections
import os, sys, subprocess
from optparse import OptionParser
	
def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	# required qsub options
	parser.add_option("-s", "--summary", action="store", dest="summary", help="text file giving paths to allele sequences (one line/file per locus)", default="")
	parser.add_option("-d", "--database", action="store", dest="database", help="MLST profile database (col1=ST, other cols=loci, must have loci names in header)", default="")
	parser.add_option("-n", "--namesep", action="store", dest="namesep", help="separator for allele names (either '-' (default) or '_')", default="-")

	return parser.parse_args()

if __name__ == "__main__":

	(options, args) = main()

	if options.database=="":
		DoError("No MLST databse provided (-d)")
		
	locus_seqs = {} # key = id (file name before extension), value = path to sequences
	if options.summary=="":
		DoError("No query sequences provided (-s)")
	else:
		f = file(options.summary,"r")
		for line in f:
			line.rstrip()
			(path,fileName) = os.path.split(line)
			if not os.path.exists(fileName + ".nin"):
				os.system("makeblastdb -dbtype nucl -logfile blast.log -in " + fileName)
			(fileName,ext) = os.path.splitext(fileName)
			locus_seqs[fileName]=line
		f.close()	
		
	sts = {} # key = concatenated string of alleles, value = st
	max_st = 0 # changeable variable holding the highest current ST, incremented when novel combinations are encountered
	header = []
	f = file(options.database,"r")
	for line in f:
		fields = line.rstrip().split("\t")
		if len(header)==0:
			header = fields
			header.pop(0) # remove st label
		else:
			sts[",".join(fields[1:])] = fields[0]
			if int(fields[0]) > max_st:
				max_st = int(fields[0])
	f.close()
	
	# check the loci match up
	for locus in header:
		if locus not in locus_seqs:
			DoError("Locus "+locus+"in ST database file " + options.database + " but no matching sequence in " + options.summary)
	
	best_match = collections.defaultdict(dict) # key1 = strain, key2 = locus, value = best match (clean for ST, annotated)
	perfect_match = collections.defaultdict(dict) # key1 = strain, key2 = locus, value = perfect match if available
	
	# print header
	print "\t".join(["strain","perfectMatchST"]+header+["bestMatchST"]+header)
	
	for contigs in args:
		tmp = contigs + ".tmp"
		(dir,fileName) = os.path.split(contigs)
		(name,ext) = os.path.splitext(fileName)
		perfect_st = []
		best_st = []
		best_st_annotated = []
		for locus in header:
			# correct order to build up concatenated ST
			locus_seq = locus_seqs[locus]
			cmd = " ".join(["blastn","-query",contigs,"-db",locus_seq.rstrip(),"-max_target_seqs","1","-outfmt '6 qseqid sacc pident length slen qlen'",">",tmp,"\n"])
			os.system(cmd)
			# read results
			if os.stat(tmp)[6]!=0:
				# file is not empty, ie match found
				match = ""
				f = file(tmp,"r")
				fields = f.readline().rstrip().split("\t") # only has one line
				(contig,allele,pcid,length,allele_length,contig_length) = (fields[0],fields[1].split(options.namesep)[1],float(fields[2]),int(fields[3]),int(fields[4]),int(fields[5]))
				if pcid==100.00 and allele_length==length:
					match = ""
					perfect_st.append(allele)
					best_st.append(allele)
					best_st_annotated.append(allele)
				elif pcid>90 and float(length)/float(allele_length) > 0.9:
					match = "/" + str(pcid) + "," + str(float(length)/float(allele_length))
					perfect_st.append("-")
					best_st.append(allele)
					best_st_annotated.append(allele + match)					
				else:
					perfect_st.append("-") # no allele at all
					best_st.append("-")
					best_st_annotated.append("-")
				f.close()
			else:
				perfect_st.append("-") # no allele at all
				best_st.append("-")
				best_st_annotated.append("-")
			os.system("rm -rf "+tmp)
			
		# assign ST
		pst = ",".join(perfect_st)
		bst = ",".join(best_st)
		
		if pst in sts:
			pst = sts[pst]
		elif "-" not in pst:
			max_st += 1 # new combination
			sts[pst] = str(max_st)
			pst = str(max_st)
		else:
			pst = "0"
		
		if bst in sts:
			bst = sts[bst]
		elif "-" not in bst:
			max_st += 1 # new combination
			sts[bst] = str(max_st)
			bst = str(max_st)
		else:
			bst = "0"

		if best_st != best_st_annotated:
			bst = "*"+bst
			
		print "\t".join([name,pst] + perfect_st + [bst] + best_st_annotated)