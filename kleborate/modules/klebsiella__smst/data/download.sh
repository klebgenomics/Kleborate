wget -O profiles.tsv "https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl?db=pubmlst_klebsiella_seqdef&page=downloadProfiles&scheme_id=21"
wget -O iroB.fasta "https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl?db=pubmlst_klebsiella_seqdef&page=downloadAlleles&locus=iroB"
wget -O iroC.fasta "https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl?db=pubmlst_klebsiella_seqdef&page=downloadAlleles&locus=iroC"
wget -O iroD.fasta "https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl?db=pubmlst_klebsiella_seqdef&page=downloadAlleles&locus=iroD"
wget -O iroN.fasta "https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl?db=pubmlst_klebsiella_seqdef&page=downloadAlleles&locus=iroN"

for f in *.fasta; do
    seqtk seq "$f" > temp.fasta && mv temp.fasta "$f"  # make FASTA files one-line-per-sequence
done
