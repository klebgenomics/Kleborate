wget -O profiles.tsv "https://pubmlst.org/bigsdb?db=pubmlst_escherichia_seqdef&page=downloadProfiles&scheme_id=2"
wget -O dinB.fasta "https://pubmlst.org/bigsdb?db=pubmlst_escherichia_seqdef&page=downloadAlleles&locus=dinB"
wget -O icdA.fasta "https://pubmlst.org/bigsdb?db=pubmlst_escherichia_seqdef&page=downloadAlleles&locus=icdA"
wget -O pabB.fasta "https://pubmlst.org/bigsdb?db=pubmlst_escherichia_seqdef&page=downloadAlleles&locus=pabB"
wget -O polB.fasta "https://pubmlst.org/bigsdb?db=pubmlst_escherichia_seqdef&page=downloadAlleles&locus=polB"
wget -O putP.fasta "https://pubmlst.org/bigsdb?db=pubmlst_escherichia_seqdef&page=downloadAlleles&locus=putP"
wget -O trpA.fasta "https://pubmlst.org/bigsdb?db=pubmlst_escherichia_seqdef&page=downloadAlleles&locus=trpA"
wget -O trpB.fasta "https://pubmlst.org/bigsdb?db=pubmlst_escherichia_seqdef&page=downloadAlleles&locus=trpB"
wget -O uidA.fasta "https://pubmlst.org/bigsdb?db=pubmlst_escherichia_seqdef&page=downloadAlleles&locus=uidA"

for f in *.fasta; do
    seqtk seq "$f" > temp.fasta && mv temp.fasta "$f"  # make FASTA files one-line-per-sequence
done
