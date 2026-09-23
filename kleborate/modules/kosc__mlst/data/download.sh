wget -O profiles.tsv "https://pubmlst.org/bigsdb?db=pubmlst_koxytoca_seqdef&page=downloadProfiles&scheme_id=1"
wget -O gapA.fasta "https://pubmlst.org/bigsdb?db=pubmlst_koxytoca_seqdef&page=downloadAlleles&locus=gapA"
wget -O infB.fasta "https://pubmlst.org/bigsdb?db=pubmlst_koxytoca_seqdef&page=downloadAlleles&locus=infB"
wget -O mdh.fasta "https://pubmlst.org/bigsdb?db=pubmlst_koxytoca_seqdef&page=downloadAlleles&locus=mdh"
wget -O pgi.fasta "https://pubmlst.org/bigsdb?db=pubmlst_koxytoca_seqdef&page=downloadAlleles&locus=pgi"
wget -O phoE.fasta "https://pubmlst.org/bigsdb?db=pubmlst_koxytoca_seqdef&page=downloadAlleles&locus=phoE"
wget -O rpoB.fasta "https://pubmlst.org/bigsdb?db=pubmlst_koxytoca_seqdef&page=downloadAlleles&locus=rpoB"
wget -O tonB.fasta "https://pubmlst.org/bigsdb?db=pubmlst_koxytoca_seqdef&page=downloadAlleles&locus=tonB"

for f in *.fasta; do
    seqtk seq "$f" > temp.fasta && mv temp.fasta "$f"  # make FASTA files one-line-per-sequence
done
