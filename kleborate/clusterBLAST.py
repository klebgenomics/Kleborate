# blast for sets of genes that make up operons for screening
# summarise hits for each operon
import os
import sys
from optparse import OptionParser


def main():

    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)

    # options
    parser.add_option("-s", "--seqs", action="store", dest="seqs", default="",
                      help="operon sequences to screen for")
    parser.add_option("-m", "--minident", action="store", dest="minident", default="90",
                      help="Minimum percent identity (default 90)")
    parser.add_option("-c", "--mincov", action="store", dest="mincov", default="80",
                      help="Minimum percent coverage (default 80)")
    return parser.parse_args()

if __name__ == "__main__":

    (options, args) = main()

    def check_dup(x):
        once = []
        twice = []
        for i in x:
            if i not in once:
                once.append(i)
            else:
                twice.append(i)
        return once, twice

    if options.seqs == "":
        sys.exit("No operon sequences provided (-s)")
    else:
        (path, fileName) = os.path.split(options.seqs)
        if not os.path.exists(options.seqs + ".nin"):
            os.system("makeblastdb -dbtype nucl -logfile blast.log -in " + options.seqs)
        (fileName, ext) = os.path.splitext(fileName)

    # print header
    print "\t".join(["strain", "aerobactin", "salmochelin", "hypermucoidy"])

    for contigs in args:
        (_, fileName) = os.path.split(contigs)
        (name, ext) = os.path.splitext(fileName)

        # blast against all
        f = os.popen("blastn -task blastn -db " + options.seqs + " -query " + contigs +
                     " -outfmt '6 sacc pident slen length score' -ungapped -dust no -evalue 1E-20 -word_size 32"
                     " -max_target_seqs 10000 -culling_limit 1 -perc_identity " + options.minident)

        # list of genes in each locus with hits
        iro = []
        rmpA = []
        iuc = []
        for line in f:
            fields = line.rstrip().split("\t")
            (gene_id, pcid, length, allele_length, score) = (fields[0], float(fields[1]), float(fields[2]),
                                                             float(fields[3]), float(fields[4]))
            if gene_id.startswith("iro"):
                if (allele_length / length * 100) > float(options.mincov):
                    iro.append(gene_id[3])  # 4th character is the gene letter
            if gene_id.startswith("iuc"):
                if (allele_length / length * 100) > float(options.mincov):
                    iuc.append(gene_id[3])  # 4th character is the gene letter
            if gene_id.startswith("rmpA"):
                if (allele_length / length * 100) > float(options.mincov):
                    rmpA.append(gene_id)
        f.close()

        iro.sort()
        iuc.sort()
        rmpA.sort()

        iro_string = "-"
        if len(iro) > 0:
            (iro, iro_dup) = check_dup(iro)
            iro_string = "iro" + "".join(iro)
            if len(iro_dup) > 0:
                iro_dup = list(set(iro_dup))  # converting to set removes duplicates
                iro_dup.sort()
                iro_string += ";iro" + "".join(iro_dup)

        iuc_string = "-"
        if len(iuc) > 0:
            (iuc, iuc_dup) = check_dup(iuc)
            iuc_string = "iuc" + "".join(iuc)
            if len(iuc_dup) > 0:
                iuc_dup = list(set(iuc_dup))  # converting to set removes duplicates
                iuc_dup.sort()
                iuc_string += ";iuc" + "".join(iuc_dup)

        rmpA_string = "-"
        if len(rmpA) > 0:
            rmpA_string = ";".join(rmpA)

        print "\t".join([name, iuc_string, iro_string, rmpA_string])
