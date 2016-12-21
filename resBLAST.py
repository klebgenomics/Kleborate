# blast for resistance genes, summarise by class (one class per column)
import string, re, collections
import os, sys, subprocess
from optparse import OptionParser
	
def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	# options
	parser.add_option("-s", "--seqs", action="store", dest="seqs", help="res gene sequences to screen for", default="ARGannot.r1.fasta")
	parser.add_option("-t", "--class", action="store", dest="res_class_file", help="res gene classes (CSV)", default="ARGannot_clustered80.csv")
	parser.add_option("-m", "--minident", action="store", dest="minident", help="Minimum percent identity (default 90)", default="90")
	parser.add_option("-c", "--mincov", action="store", dest="mincov", help="Minimum percent coverage (default 80)", default="80")
	
	return parser.parse_args()

if __name__ == "__main__":

	(options, args) = main()
		
	if options.seqs=="":
		DoError("No res gene sequences provided (-s)")
	else:
		(path,fileName) = os.path.split(options.seqs)
		if not os.path.exists(options.seqs + ".nin"):
			os.system("makeblastdb -dbtype nucl -logfile blast.log -in " + options.seqs)
		(fileName,ext) = os.path.splitext(fileName)
		
	# read table of genes and store classes
	
	gene_info = {} # key = sequence id (fasta header in seq file), value = (allele,class,Bla_Class)
	res_classes = []
	bla_classes = []

	if options.res_class_file=="":
		DoError("No res gene class file provided (-t)")
	else:
		f = file(options.res_class_file,"r")
		header=0
		for line in f:
			if header == 0:
				header = 1
				#seqID,clusterid,gene,allele,cluster_contains_multiple_genes,gene_found_in_multiple_clusters,idInFile,symbol,class,accession,positions,size,Lahey,Bla_Class
			else:
				fields = line.rstrip().split(",")
				(seqID, clusterID, gene, allele, allele_symbol, res_class, bla_class) = (fields[0], fields[1], fields[2], fields[3], fields[3], fields[8], fields[13])
				seq_header = "__".join([clusterID,gene,allele,seqID])
				if res_class == "Bla" and bla_class == "NA":
					bla_class = "Bla"
				gene_info[seq_header] = (allele_symbol, res_class, bla_class)
				if res_class not in res_classes:
					res_classes.append(res_class)
				if bla_class not in bla_classes:
					bla_classes.append(bla_class)
		f.close()
	
	res_classes.sort()
	res_classes.remove("Bla")
	bla_classes.sort()
	bla_classes.remove("NA")
		
	# print header
	print "\t".join(["strain"] + res_classes + bla_classes)

	for contigs in args:
		(dir,fileName) = os.path.split(contigs)
		(name,ext) = os.path.splitext(fileName)

		# blast against all
		f = os.popen("blastn -task blastn -db " + options.seqs + " -query " + contigs + " -outfmt '6 sacc pident slen length score' -ungapped -dust no -evalue 1E-20 -word_size 32 -max_target_seqs 10000 -culling_limit 1 -perc_identity " + options.minident)

		# list of genes in each class with hits
		hits_dict = {} # key = class, value = list
		
		for line in f:
			fields = line.rstrip().split("\t")
			(gene_id,pcid,length,allele_length,score) = (fields[0],float(fields[1]),float(fields[2]),float(fields[3]),float(fields[4]))
			if (allele_length/length*100) > float(options.mincov):
				(hit_allele, hit_class, hit_bla_class) = gene_info[gene_id]
				if hit_class == "Bla":
					hit_class = hit_bla_class
				if pcid < 100.00:
					hit_allele += "*" # imprecise match
				if allele_length < length:
					hit_allele += "?" # partial match
				if hit_class in hits_dict:
					hits_dict[hit_class].append(hit_allele)
				else:
					hits_dict[hit_class] = [hit_allele]
		f.close()
		
		hit_string = [name]
		for res_class in (res_classes + bla_classes):
			if res_class in hits_dict:
				hit_string.append(";".join(hits_dict[res_class]))
			else:
				hit_string.append("-")
				
		print "\t".join(hit_string)