# NTUH-K2044 (ST23, ybt 2; ICEKp1)
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/009/885/GCA_000009885.1_ASM988v1/GCA_000009885.1_ASM988v1_genomic.fna.gz
gunzip GCA_000009885.1_ASM988v1_genomic.fna.gz
mv GCA_000009885.1_ASM988v1_genomic.fna NTUH-K2044.fna

# Kp1084 (ST23, ybt 1; ICEKp10, clb 2)
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/294/365/GCA_000294365.1_ASM29436v1/GCA_000294365.1_ASM29436v1_genomic.fna.gz
gunzip GCA_000294365.1_ASM29436v1_genomic.fna.gz
mv GCA_000294365.1_ASM29436v1_genomic.fna Klebs_Kp1084.fna

# HS11286 (ST11, ybt 9; ICEKp3)
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/240/185/GCA_000240185.2_ASM24018v2/GCA_000240185.2_ASM24018v2_genomic.fna.gz
gunzip GCA_000240185.2_ASM24018v2_genomic.fna.gz
mv GCA_000240185.2_ASM24018v2_genomic.fna Klebs_HS11286.fna

# MGH 78578 (ST38, no ICE)
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/016/305/GCA_000016305.1_ASM1630v1/GCA_000016305.1_ASM1630v1_genomic.fna.gz
gunzip GCA_000016305.1_ASM1630v1_genomic.fna.gz
mv GCA_000016305.1_ASM1630v1_genomic.fna MGH78578.fna

# run typing
# NOTE: -p must point to the Kleborate directory
kleborate -p . -o details.txt *.fna
