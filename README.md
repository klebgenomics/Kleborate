# Kleborate

_Klebsiella pneumoniae_ (_Kp_) is a commensal bacterium that causes opportunistic infections, with a handful of hypervirulent lineages recognised as true human pathogens. Evidence is now mounting that other _Kp_ strains carrying acquired siderophores (yersiniabactin, salmochelin and aerobactin) and/or the genotoxin colibactin are also highly pathogenic and can cause invasive disease. 

We recently explored the diversity of the _Kp_ integrative conjugative element (ICEKp), which mobilises the yersiniabactin locus _ybt_, using genomic analysis of a diverse set of 2499 _Kp_. Overall, we found _ybt_ in about a third of all _Kp_ genomes and _clb_ in about 5%.

We identified 17 distinct lineages of _ybt_ embedded within 14 structural variants of ICEKp (some of which include the colibactin _clb_ or salmochelin _iro_ synthesis loci, annotated reference sequences for each ICEKp variant are included in the /data directory of this repository) that can integrate at any of four tRNA-Asn sites in the chromosome. Our analyses reveal hundreds of ICEKp transmission events affecting hundreds of chromosomal _Kp_ lineages, including nearly two dozen transfers into the globally disseminated carbapenem-resistant clonal group 258. Additionally, we identify a lineage of _ybt_ that is plasmid-encoded, representing a new mechanism for _ybt_ dispersal in _Kp_ populations. 

Based on this analysis, we developed a MLST-style approach for identifying _ybt_ and _clb_ variants from genome data. 

Our goal is to help identify emerging pathogenic _Kp_ lineages, and to make it easy for people who are using genomic surveillance to monitor for antibiotic resistance to also look out for the convergence of antibiotic resistance and virulence.

To help facilitate that, in this repo we share the new _ybt_ and _clb_ schemes (/data), annotated ICEKp structures (/ICEKp_references) and code for genotyping virulence and resistance genes in _K. pneumoniae_. A table of pre-computed results for 2500 public Klebs genomes is also provided in the /data directory.

If you are interested in inferring capsule types from genome data, see the [Kaptive](https://github.com/katholt/Kaptive) repo.

## Let's get genotyping!
Just want to get cracking with screening a bunch of _K. pneumoniae_ genome assemblies? Use the Kleborate.py script. This will detect the MLST sequence type of the strain, genotype the _ybt_ and _clb_ loci, determine the _wzi_ (capsule synthesis gene) allele and also check for presence/absence of the acquired siderophores salmochelin (_iro_) and aerobactin (_iuc_) loci and the hypermucoidy genes _rmpA_ and _rmpA2_ (allelic typing of these should be available soon). For convenience, we provide code for screening for acquired resistance genes (resBLAST.py) and quinolone-resistance determining mutations in _gyrA_ and _parC_, which can optionally be called when you run Kleborate.py or as a standalone script.

#### Basic usage:

```
# get the code
git clone https://github.com/katholt/Kleborate 

# screen some genomes for MLST + virulence loci
python Kleborate.py -p Kleborate -o detailed_results.txt *.fasta

# screen some genomes for MLST + virulence loci + acquired resistance genes
python Kleborate.py -p Kleborate -o detailed_results.txt -r on *.fasta
```

See below for more details, examples and outputs.

#### Dependencies

* Python v2

* BLAST+ v2.2.30 (note earlier versions have a bug with the culling_limit parameter)

* EMBOSS (optional, to help identify quinolone resistant SNPs)


# About the MLST schemes

We have created two separate schemes: one for yersiniabactin sequence types (YbST) and one for colibactin sequence types (CbST).

MLST-style schemes are included in the [_Klebsiella pneumoniae_ BIGSdb hosted at the Pasteur Institute](http://bigsdb.pasteur.fr/klebsiella/klebsiella.html), and are also included in the /data directory of this repository. The schemes include all known alleles for the genes that make up the yersiniabactin and colibactin synthesis loci that are mobilised by the _Klebsiella_ ICEKp.

### YbST - Yersiniabactin Sequence Types

YbST sequences cluster into 17 distinct lineages of _ybt_, which are each associated with a particular ICEKp structure (see the paper for details). Three lineages (_ybt_ 1, _ybt_ 12, _ybt_ 17) are associated with the ICEKp10 structure in which the _clb_ colibactin locus is found; _ybt_ 4 is plasmid-encoded; and all other lineages correspond to a single ICEKp.

### CbST - Colibactin Sequence Types

CbST sequences cluster into 3 lineages, which are each associated with a single _ybt_ lineage (_clb_ 1 / _ybt_ 12; _clb_ 2, _ybt_ 1; ; _clb_ 3, _ybt_ 17) and the ICEKp10 structure.

## Genotypes of publicly available strains

A table of pre-computed yersiniabactin, colibactin,  capsule locus and chromosomal MLST assignments for 2500 public Klebs genomes is provided in the /data directory.

# Typing genome assemblies using Kleborate

The Kleborate.py script in this repo can be used to determine chromosomal, yersiniabactin and colibactin genotypes from assembled draft or complete genomes. It also reports presence/absence of acquired siderophores salmochelin (_iro_) and aerobactin (_iuc_) loci and the hypermucoidy genes _rmpA_ and _rmpA2_ (allelic typing of these should be available soon, currently we just screen for the alleles in the virulence plasmid pLVPK). We also extract the _wzi_ gene allele to give an idea of the capsule type, but these are not totally predictive so we suggest you use our dedicated capsule typing tool [Kaptive](https://github.com/katholt/Kaptive) for this. For convenience, Kleborate.py can optionally screen for acquired resistance genes as well (using the SRST2-formatted version of the ARG-Annot database, with the core Klebs genes _oqxA_ and _oxqB_ removed).

A summary of sequence types and ICE/lineage information is printed to standard out; full allele calls are saved to a file specified by `-o`. See below for details of output formats.

#### Dependencies

* Python v2

* BLAST+ v2.2.30 (note earlier versions have a bug with the culling_limit parameter)

* EMBOSS (optional, to help identify quinolone resistant SNPs)

#### To download:

````
git clone https://github.com/katholt/Kleborate 
```

#### Usage:

```
python Kleborate.py -h
Usage: Kleborate.py [options]

Options:
  -h, --help            show this help message and exit
  -p REPO_PATH, --path=REPO_PATH
                        Path to Kleborate directory (default Kleborate)
  -o OUTFILE, --outfile=OUTFILE
                        File for detailed output (default Kleborate_results.txt)
```

#### Example command:

```
python Kleborate.py -p Kleborate -o detailed_results.txt genome.fasta
```

#### Test on well known genomes:

```
## GET CODE

git clone https://github.com/katholt/Kleborate 

## GET DATA FROM NCBI

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

# RUN KLEBORATE
# NOTE: -p must point to the Kleborate directory

python Kleborate.py -p Kleborate -o details.txt *.fna

## EXPECTED OUTPUT
strain	ST	Yersiniabactin	YbST	Colibactin	CbST	aerobactin	salmochelin	hypermucoidy	wzi 
Klebs_HS11286	11	ybt 9; ICEKp3	15	-	0	-	-	-	wzi74 
Klebs_Kp1084	23	ybt 1; ICEKp10	47	clb 2	37	-	-	-	wzi172 
MGH78578	38	-	0	-	0	-	-	-	wzi50 
NTUH-K2044	23	ybt 2; ICEKp1	326	-	0	iucABCD	iroBCDN;iroBCDN	rmpA;rmpA	wzi1

# NOTE: NTUH-K2044 has two copies each of the iro and rmpA loci (one on the virulence plasmid with iuc, and one in ICEKp1).

# RUN KLEBORATE with resistance gene screening turned on:

python Kleborate.py -p Kleborate -o details.txt *.fna -r on

```
## Output

A tabulated summary is printed to standard out; details of the MLST analysis are printed to the file specified by -o.

#### Klebs, yersiniabactin and colibactin MLST
* Kleborate makes an effort to report the closest matching ST / clonal group if a precise match is not found.
* Imprecise allele matches are indicated with a \*
* Imprecise ST calls are indicated with -nLV, where n indicates the number of loci that disagree with the ST reported. So 258-1LV indicates a single-locus variant of (SLV) of ST258, i.e. 6/7 loci match the ST258 alleles.
* For _ybt_ and _clb_ schemes, the associated lineage and/or ICEKp structure are reported along with the ST. If these loci are not detected, the ST reported is 0.

#### Virulence locus detection
* Results in these columns should be interpreted as simply binary yes/no calls.
* By default, a gene is called as present if it is detected in a single sequence with >90% identity and >80% coverage of the allele sequence from the virulence plasmid pLVPK. Note that rmpA and rmpA2 are ~85% homologous and are reported separately.
* If multiple hits to the same query sequence are found in a given assembly, we attempt to report these separately. The NTUH-K2044 genome carries iroBCDN and rmpA in two locations (one on the virulence plasmid, one in the ICEKp1), and this should be reported as iroBDCN;iroBCDN and rmpA;rmpA as seen in the example above.

#### Wzi gene allele (marker of capsule type)
* The closest match amongst the _wzi_ alleles in the BIGSdb will be reported.
* This is a marker of capsule (K) type, although there is not a 1-1 relationship between wzi allele and K type
* This can a handy way of spotting the virulence-associated types (wzi=K1, wzi2=K2, wzi5=K5); or spotting capsule switching within clones, e.g. you can tell which ST258 lineage you have from the wzi type (wzi154: the main lineage II; wzi29: recombinant lineage I; others: probably other recombinant lineages)
* Note for proper capsule type prediction you should use our dedicated capsule typing tool [Kaptive](https://github.com/katholt/Kaptive)

#### Resistance gene detection
* Here we are screening against the ARG-Annot database of acquired resistance genes ([SRST2](https://github.com/katholt/srst2) version), which includes allelic variants.
* Kleborate attempts to report the best matching variant for each locus in the genome
* Imprecise allele matches are indicated with \*
* If the length of match is less than the length of the reported allele (ie a partial match), this is indicated with ?
* Note that the beta-lactamases ampH and SHV (narrow spectrum) are core genes in _K. pneumoniae_ so should be detected in most genomes
* If you see LEN or OKP beta-lactamases rather than SHV, you probably have _K. variicola_ (LEN) or _K. quasipneumoniae_ (OKP) rather than _K. pneumoniae_ (see [this paper](http://www.pnas.org/content/112/27/E3574.long) for clarification)
* Note that _oqxAB_ are also core genes  in _K. pneumoniae_, but have been removed from this version of the ARG-Annot DB as they don't actually confer resistance to fluoroquinolones
* In addition to acquired genes, we also check for the known quinolone resistance SNPs (GyrA 83 & 87; ParC 80 & 84)
* Results are grouped by drug class (according to the [ARG-Annot](https://www.ncbi.nlm.nih.gov/pubmed/24145532) DB), with beta-lactamases broken down into Lahey classes, as follows: 
  * AGly (aminoglycosides)
  * Bla (beta-lactamases)
  * Bla_broad (broad spectrum beta-lactamases)
  * Bla_broad_inhR (broad spectrum beta-lactamases with resistance to beta-lactamase inhibitors)
  * Bla_Carb (carbapenemase)
  * Bla_ESBL (extended spectrum beta-lactamases)
  * Bla_ESBL_inhR (extended spectrum beta-lactamases with resistance to beta-lactamase inhibitors)
  * Fcyn (fosfomycin)
  * Flq (fluoroquinolones)
  * Gly (glycopeptides)
  * MLS (macrolides)
  * Phe (phenicols)
  * Rif (rifampin)
  * Sul (sulfonamides)
  * Tet (tetracyclines)
  * Tmt (trimethoprim)


# Typing direct from Illumina reads

MLST assignment can also be achieved direct from reads using [SRST2](https://github.com/katholt/srst2). Steps are 

* download the YbST and CbST allele sequences and profile tables from the /data directory in _this_ repository
* Install [SRST2](https://github.com/katholt/srst2) if you don't already have it (`git clone https://github.com/katholt/srst2`); 
* Run SRST2, setting the `--mlst_scheme` and `--mlst_definitions` to point to the YbST or CbST allele sequences and profile tables like this:

```
srst2 --input_pe reads_1.fastq.gz reads_2.fastq.gz --output YbST --log --mlst_db ybt_alleles.fasta --mlst_definitions YbST_profiles.txt

srst2 --input_pe reads_1.fastq.gz reads_2.fastq.gz --output CbST --log --mlst_db colibactin_alleles.fasta --mlst_definitions CbST_profiles.txt
```

Note that currently you can only run SRST2 with one MLST scheme at a time, so in order to type MLST, YbST and CbST you will need to run three separate commands:

```
srst2 --input_pe reads_1.fastq.gz reads_2.fastq.gz --output YbST --log --mlst_db ybt_alleles.fasta --mlst_definitions YbST_profiles.txt

srst2 --input_pe reads_1.fastq.gz reads_2.fastq.gz --output CbST --log --mlst_db colibactin_alleles.fasta --mlst_definitions CbST_profiles.txt

srst2 --input_pe reads_1.fastq.gz reads_2.fastq.gz --output Klebs --log --mlst_db Klebsiella_pneumoniae.fasta --mlst_definitions kpnuemoniae.txt
```


## Kleboration

Kleborate is under active development with many other Klebs genomic analysis tools and projects in progress. 

Please get in touch via the issues tracker if you have any issues, questions or ideas.

-------------

Stop! Kleborate and listen

ICEKp is back with with my brand new invention

If there was a problem, Klebs'll solve it

Check out the hook while Klebs evolves it
