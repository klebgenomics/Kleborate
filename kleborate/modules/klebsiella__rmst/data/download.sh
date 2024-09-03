wget -O profiles.tsv "https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl?db=pubmlst_klebsiella_seqdef&page=downloadProfiles&scheme_id=23"
wget -O rmpA.fasta "https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl?db=pubmlst_klebsiella_seqdef&page=downloadAlleles&locus=rmpA"
wget -O rmpC.fasta "https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl?db=pubmlst_klebsiella_seqdef&page=downloadAlleles&locus=rmpC"
wget -O rmpD.fasta "https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl?db=pubmlst_klebsiella_seqdef&page=downloadAlleles&locus=rmpD"

for f in *.fasta; do
    seqtk seq "$f" > temp.fasta && mv temp.fasta "$f"  # make FASTA files one-line-per-sequence
done
