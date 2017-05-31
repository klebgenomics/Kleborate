# reports best match
# * : best matching allele is not precise match
# -nLV : best matching ST is n-locus variant
# if an annotation column is provided (such as clonal complex) in the final column of the profiles file,
#     this annotation will be reported in column 2 of the output table
# NOTE there is a bug with the culling_limit parameter in older versions of BLAST+...
# This code has been tested with BLAST+2.2.30. It does not work with BLAST2.2.25. Not sure about other versions.

import collections
import os
import sys
from optparse import OptionParser


def main():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)

    # options
    parser.add_option("-s", "--seqs", action="store", dest="seqs", default="",
                      help="MLST allele sequences file")
    parser.add_option("-d", "--database", action="store", dest="database", default="",
                      help="MLST profile database (col1=ST, other cols=loci, must have loci names in header)")
    parser.add_option("-i", "--info", action="store", dest="info", default="yes",
                      help="Info (clonal group, lineage, etc) provided in last column of profiles (yes (default), no)")
    parser.add_option("-m", "--minident", action="store", dest="minident", default="95",
                      help="Minimum percent identity (default 95)")
    parser.add_option("-n", "--maxmissing", action="store", dest="maxmissing", default="3",
                      help="Maximum missing/uncalled loci to still calculate closest ST (default 3)")
    return parser.parse_args()

if __name__ == "__main__":
    (options, args) = main()

    if options.database == "":
        sys.exit("No MLST profiles databse provided (-d)")

    if options.seqs == "":
        sys.exit("No MLST allele sequences provided (-s)")
    else:
        (path, fileName) = os.path.split(options.seqs)
        if not os.path.exists(options.seqs + ".nin"):
            os.system("makeblastdb -dbtype nucl -logfile blast.log -in " + options.seqs)
        (fileName, ext) = os.path.splitext(fileName)

    def get_closest_locus_variant(query, annotated_query, sts):
        annotated_query = list(annotated_query)  # copy the list so we don't change the original
        closest = []
        closest_alleles = {} # key = st, value = list
        min_dist = len(query) # number mismatching loci, ignoring SNPs

        for index, item in enumerate(query):
            if item == "-":
                query[index] = "0"

        # get distance from closest ST, ignoring SNPs (*)
        for st in sts:
            d = sum(map(lambda x, y: bool(int(x)-int(y)), st.split(","), query))
            if d == min_dist:
                closest.append(int(sts[st]))
                closest_alleles[sts[st]] = st
            elif d < min_dist:
                # reset
                closest = [int(sts[st])]
                closest_alleles[sts[st]] = st
                min_dist = d # distance from closest ST, ignoring SNPs (*)

        closest_st = str(min(closest))

        for index, item in enumerate(annotated_query):
            if item == "-" or item.endswith("*"):
                annotated_query[index] = "0"

        # get distance from closest ST, including SNPs (*)
        min_dist_incl_snps = sum(map(lambda x, y: bool(int(x)-int(y)),
                                     closest_alleles[closest_st].split(","), annotated_query))

        return closest_st, min_dist, min_dist_incl_snps

    sts = {}  # key = concatenated string of alleles, value = st
    st_info = {}  # key = st, value = info relating to this ST, eg clonal group
    max_st = 0  # changeable variable holding the highest current ST, incremented when novel combinations are encountered
    header = []
    info_title = "info"
    f = file(options.database, "r")
    for line in f:
        fields = line.rstrip().split("\t")
        if len(header)==0:
            header = fields
            header.pop(0)  # remove st label
            if options.info == "yes":
                info_title = header.pop() # remove info label
        else:
            st = fields.pop(0)
            if options.info == "yes":
                info = fields.pop()
            sts[",".join(fields)] = st
            if int(st) > max_st:
                max_st = int(st)
            if options.info == "yes":
                st_info[st] = info
    f.close()

    best_match = collections.defaultdict(dict)  # key1 = strain, key2 = locus, value = best match (clean for ST, annotated)
    perfect_match = collections.defaultdict(dict)  # key1 = strain, key2 = locus, value = perfect match if available

    # print header
    if options.info == "yes":
        print "\t".join(["strain", info_title, "ST"] + header)
    else:
        print "\t".join(["strain", "ST"] + header)

    for contigs in args:
        (_, fileName) = os.path.split(contigs)
        (name, ext) = os.path.splitext(fileName)

        best_score = {}  # key = locus, value = BLAST score for best matching allele encountered so far
        best_allele = {}  # key = locus, value = best allele (* if imprecise match)

        # blast against all
        f = os.popen("blastn -task blastn -db " + options.seqs + " -query " + contigs +
                     " -outfmt '6 sacc pident slen length score' -ungapped -dust no -evalue 1E-20 -word_size 32"
                     " -max_target_seqs 10000 -culling_limit 2 -perc_identity " + options.minident)
        for line in f:
            fields = line.rstrip().split("\t")
            (gene_id, pcid, length, allele_length, score) = (fields[0], float(fields[1]), int(fields[2]),
                                                             int(fields[3]), float(fields[4]))
            if "__" in gene_id:
                # srst2 formatted file
                gene_id_components = gene_id.split("__")
                locus = gene_id_components[1]
                allele = gene_id_components[2]
            else:
                allele = gene_id
                locus = gene_id.split("_")[0]
            if pcid < 100.00 or allele_length < length:
                allele += "*"  # imprecise match
            # store best match for each one locus
            if locus in best_score:
                if score > best_score[locus]:
                    # update
                    best_score[locus] = score
                    best_allele[locus] = allele.split("_")[1]  # store number only
            else:
                # initialise
                best_score[locus] = score
                best_allele[locus] = allele.split("_")[1]  # store number only
        f.close()

        best_st = []
        best_st_annotated = []

        mismatch_loci = 0
        mismatch_loci_including_SNPs = 0

        for locus in header:
            if locus in best_allele:
                allele = best_allele[locus]
                allele_number = allele.replace("*", "")
                if allele.endswith("*"):
                    mismatch_loci_including_SNPs += 1
                best_st.append(allele_number)
                best_st_annotated.append(allele) # will still have * if imperfect match
            else:
                best_st.append("-")
                best_st_annotated.append("-")
                mismatch_loci += 1
                mismatch_loci_including_SNPs += 1

        # assign ST
        bst = ",".join(best_st)

        if mismatch_loci_including_SNPs <= int(options.maxmissing):
            # only report ST if enough loci are precise matches
            if bst in sts:
                bst = sts[bst]  # note may have mismatching alleles due to SNPs, this will be recorded in mismatch_loci_including_SNPs
            else:
                # determine closest ST
                (bst, mismatch_loci, mismatch_loci_including_SNPs) = get_closest_locus_variant(best_st,
                                                                                               best_st_annotated, sts)
        else:
            bst = "0"

        # pull info column
        if options.info == "yes":
            info_final = ""
            if bst in st_info:
                info_final = st_info[bst]

        if mismatch_loci_including_SNPs > 0 and bst != "0":
            bst += "-" + str(mismatch_loci_including_SNPs) + "LV"

        if options.info == "yes":
            print "\t".join([name, info_final, bst] + best_st_annotated)
        else:
            print "\t".join([name, bst] + best_st_annotated)
