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

# blast for sets of genes that make up operons for screening
# summarise hits for each operon

import os
import sys
import subprocess
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
            with open(os.devnull, 'w') as devnull:
                subprocess.check_call("makeblastdb -dbtype nucl -in " + options.seqs,
                                      stdout=devnull, shell=True)
        (fileName, ext) = os.path.splitext(fileName)

    # print header
    print("\t".join(["strain", "hypermucoidy"]))

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
            if gene_id.startswith("rmpA"):
                if (allele_length / length * 100) > float(options.mincov):
                    rmpA.append(gene_id)
        f.close()

        rmpA.sort()

        rmpA_string = "-"
        if len(rmpA) > 0:
            rmpA_string = ";".join(rmpA)

        print("\t".join([name, rmpA_string]))
