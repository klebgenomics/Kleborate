wget -O profiles.tsv "https://pubmlst.org/bigsdb?db=pubmlst_escherichia_seqdef&page=downloadProfiles&scheme_id=1"
wget -O adk.fasta "https://pubmlst.org/bigsdb?db=pubmlst_escherichia_seqdef&page=downloadAlleles&locus=adk"
wget -O fumC.fasta "https://pubmlst.org/bigsdb?db=pubmlst_escherichia_seqdef&page=downloadAlleles&locus=fumC"
wget -O gyrB.fasta "https://pubmlst.org/bigsdb?db=pubmlst_escherichia_seqdef&page=downloadAlleles&locus=gyrB"
wget -O icd.fasta "https://pubmlst.org/bigsdb?db=pubmlst_escherichia_seqdef&page=downloadAlleles&locus=icd"
wget -O mdh.fasta "https://pubmlst.org/bigsdb?db=pubmlst_escherichia_seqdef&page=downloadAlleles&locus=mdh"
wget -O purA.fasta "https://pubmlst.org/bigsdb?db=pubmlst_escherichia_seqdef&page=downloadAlleles&locus=purA"
wget -O recA.fasta "https://pubmlst.org/bigsdb?db=pubmlst_escherichia_seqdef&page=downloadAlleles&locus=recA"

for f in *.fasta; do
    seqtk seq "$f" > temp.fasta && mv temp.fasta "$f"  # make FASTA files one-line-per-sequence
done
