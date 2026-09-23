wget -O profiles.tsv "https://bigsdb.pasteur.fr/api/db/pubmlst_klebsiella_seqdef/schemes/1/profiles_csv"
wget -O gapA.fasta "https://bigsdb.pasteur.fr/api/db/pubmlst_klebsiella_seqdef/loci/gapA/alleles_fasta"
wget -O infB.fasta "https://bigsdb.pasteur.fr/api/db/pubmlst_klebsiella_seqdef/loci/infB/alleles_fasta"
wget -O mdh.fasta "https://bigsdb.pasteur.fr/api/db/pubmlst_klebsiella_seqdef/loci/mdh/alleles_fasta"
wget -O pgi.fasta "https://bigsdb.pasteur.fr/api/db/pubmlst_klebsiella_seqdef/loci/pgi/alleles_fasta"
wget -O phoE.fasta "https://bigsdb.pasteur.fr/api/db/pubmlst_klebsiella_seqdef/loci/phoE/alleles_fasta"
wget -O rpoB.fasta "https://bigsdb.pasteur.fr/api/db/pubmlst_klebsiella_seqdef/loci/rpoB/alleles_fasta"
wget -O tonB.fasta "https://bigsdb.pasteur.fr/api/db/pubmlst_klebsiella_seqdef/loci/tonB/alleles_fasta"

for f in *.fasta; do
    seqtk seq "$f" > temp.fasta && mv temp.fasta "$f"  # make FASTA files one-line-per-sequence
done
