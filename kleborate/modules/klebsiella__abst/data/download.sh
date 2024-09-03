wget -O profiles.tsv "https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl?db=pubmlst_klebsiella_seqdef&page=downloadProfiles&scheme_id=20"
wget -O iucA.fasta "https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl?db=pubmlst_klebsiella_seqdef&page=downloadAlleles&locus=iucA"
wget -O iucB.fasta "https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl?db=pubmlst_klebsiella_seqdef&page=downloadAlleles&locus=iucB"
wget -O iucC.fasta "https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl?db=pubmlst_klebsiella_seqdef&page=downloadAlleles&locus=iucC"
wget -O iucD.fasta "https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl?db=pubmlst_klebsiella_seqdef&page=downloadAlleles&locus=iucD"
wget -O iutA.fasta "https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl?db=pubmlst_klebsiella_seqdef&page=downloadAlleles&locus=iutA"

for f in *.fasta; do
    seqtk seq "$f" > temp.fasta && mv temp.fasta "$f"  # make FASTA files one-line-per-sequence
done
