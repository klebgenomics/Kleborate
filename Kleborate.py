# run chromosome, yersiniabactin and colibactin MLST on a Klebs genome
import string, re, collections
import os, sys, subprocess
from optparse import OptionParser
	
def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	# options
	parser.add_option("-p", "--path", action="store", dest="repo_path", help="Path to Kleborate directory (default Kleborate)", default="Kleborate")
	parser.add_option("-o", "--outfile", action="store", dest="outfile", help="File for detailed output (default Kleborate_results.txt)", default="Kleborate_results.txt")
	
	return parser.parse_args()

if __name__ == "__main__":

	(options, args) = main()

	print "\t".join(["strain","ST","Yersiniabactin","YbST","Colibactin","CbST","aerobactin","salmochelin","hypermucoidy"])
	
	o = file(options.outfile, "w")
	o.write("\t".join(["strain","ST","Yersiniabactin","YbST","Colibactin","CbST","aerobactin","salmochelin","hypermucoidy","Chr_ST","gapA","infB","mdh","pgi","phoE","rpoB","tonB","YbST","ybtS","ybtX","ybtQ","ybtP","ybtA","irp2","irp1","ybtU","ybtT","ybtE","fyuA","CbST","clbA","clbB","clbC","clbD","clbE","clbF","clbG","clbH","clbI","clbL","clbM","clbN","clbO","clbP","clbQ"]))
	o.write("\n")

	for contigs in args:
		(dir,fileName) = os.path.split(contigs)
		(name,ext) = os.path.splitext(fileName)
		
		f = os.popen("python "+ options.repo_path + "/mlstBLAST.py -s "+ options.repo_path + "/data/Klebsiella_pneumoniae.fasta -d "+ options.repo_path + "/data/kpneumoniae.txt -i no --maxmissing 3 " + contigs) 

		# run chromosome MLST
		chr_ST = ""
		chr_ST_detail = []
		
		for line in f:
			fields = line.rstrip().split("\t")
			if fields[1] != "ST":
				# skip header
				(strain, chr_ST) = (fields[0], fields[1])
				chr_ST_detail = fields[2:]
		f.close()
		
		# run ybt MLST
		
		f = os.popen("python "+ options.repo_path + "/mlstBLAST.py -s "+ options.repo_path + "/data/ybt_alleles.fasta -d "+ options.repo_path + "/data/YbST_profiles.txt -i yes --maxmissing 7 " + contigs) 

		Yb_ST = ""
		Yb_group = ""
		Yb_ST_detail = []
		
		for line in f:
			fields = line.rstrip().split("\t")
			if fields[2] != "ST":
				# skip header
				(strain,Yb_ST, Yb_group) = (fields[0],fields[2], fields[1])
				Yb_ST_detail = fields[3:]
				if Yb_group == "":
					Yb_group = "-"
		f.close()
		
		# run colibactin MLST
		
		f = os.popen("python "+ options.repo_path + "/mlstBLAST.py -s "+ options.repo_path + "/data/colibactin_alleles.fasta -d "+ options.repo_path + "/data/CbST_profiles.txt -i yes --maxmissing 10 " + contigs) 

		Cb_ST = ""
		Cb_group = ""
		Cb_ST_detail = []
		
		for line in f:
			fields = line.rstrip().split("\t")
			if fields[2] != "ST":
				# skip header
				(strain,Cb_ST, Cb_group) = (fields[0],fields[2], fields[1])
				Cb_ST_detail = fields[3:]
				if Cb_group == "":
					Cb_group = "-"
		f.close()
		
		# screen for other virulence genes (binary calls)

		f = os.popen("python "+ options.repo_path + "/clusterBLAST.py -s "+ options.repo_path + "/data/other_vir_clusters.fasta " + contigs) 
		
		for line in f:
			fields = line.rstrip().split("\t")
			if fields[1] != "aerobactin":
				# skip header
				(strain,vir_hits) = (fields[0],"\t".join(fields[1:]))
		f.close()	
		
		print "\t".join([name,chr_ST,Yb_group,Yb_ST,Cb_group,Cb_ST,vir_hits])
		
		o.write("\t".join([name,chr_ST,Yb_group,Yb_ST,Cb_group,Cb_ST,vir_hits,chr_ST]+chr_ST_detail+[Yb_ST]+Yb_ST_detail + [Cb_ST] + Cb_ST_detail))
		o.write("\n")

		# run Kaptive
		
	o.close()