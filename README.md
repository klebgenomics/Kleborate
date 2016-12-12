# Kleborate

_Klebsiella pneumoniae_ (_Kp_) is a commensal bacterium that causes opportunistic infections, with a handful of hypervirulent lineages recognised as true human pathogens. Evidence is now mounting that other _Kp_ strains carrying acquired siderophores (yersiniabactin, salmochelin and aerobactin) and/or the genotoxin colibactin are also highly pathogenic and can cause invasive disease. 

We recently explored the diversity of the _Kp_ integrative conjugative element (ICEKp), which mobilises the yersiniabactin locus _ybt_, using genomic analysis of a diverse set of 2499 _Kp_. Overall, we found _ybt_ in about a third of all _Kp_ genomes and _clb_ in about 5%.

We identified 17 distinct lineages of _ybt_ embedded within 14 structural variants of ICEKp (some of which include the colibactin _clb_ or salmochelin _iro_ synthesis loci) that can integrate at any of four tRNA-Asn sites in the chromosome. Our analyses reveal hundreds of ICEKp transmission events affecting hundreds of chromosomal _Kp_ lineages, including nearly two dozen transfers into the globally disseminated carbapenem-resistant clonal group 258. Additionally, we identify a lineage of _ybt_ that is plasmid-encoded, representing a new mechanism for _ybt_ dispersal in _Kp_ populations. 

Based on this analysis, we developed a MLST-style approach for identifying _ybt_ and _clb_ variants based on genome data. Our goal is to help identify emerging pathogenic _Kp_ lineages, and to make it easy for people who are using genomic surveillance to monitor for antibiotic resistance to also look out for the convergence of antibiotic resistance and virulence.


## MLST schemes

We have created two separate schemes: one for yersiniabactin sequence types (YbST) and one for colibactin sequence types (CbST).

MLST-style schemes are included in the [_Klebsiella pneumoniae_ BIGSdb hosted at the Pasteur Institute](http://bigsdb.pasteur.fr/klebsiella/klebsiella.html), and are also included in the /data directory of this repository. The schemes include all known alleles for the genes that make up the yersiniabactin and colibactin synthesis loci that are mobilised by the _Klebsiella_ ICE.

### YbST - Yersiniabactin Sequence Types

YbST sequences cluster into 17 distinct lineages of _ybt_, which are each associated with a particular ICEKp structure (see the paper for details). Three lineages (_ybt_ 1, _ybt_ 12, _ybt_ 17) are associated with the ICEKp10 structure in which the _clb_ colibactin locus is found; _ybt_ 4 is plasmid-encoded; and all other lineages correspond to a single ICEKp.

### CbST - Colibactin Sequence Types

CbST sequences cluster into 3 lineages, which are each associated with a single _ybt_ lineage (_clb_ 1 / _ybt_ 12; _clb_ 2, _ybt_ 1; ; _clb_ 3, _ybt_ 17) and the ICEKp10 structure.

## Genotypes of publicly available strains

A table of pre-computed yersiniabactin, colibactin,  capsule locus and chromosomal MLST assignments for 2500 public Klebs genomes is provided in the /data directory.

## Typing genome assemblies

Assembled draft or complete genomes can be analysed using the script mlstBLAST.py in this repository

Dependencies: Python v2; BLAST+ v2.2.30 (note earlier versions have a bug with the culling_limit parameter)

Example command:

````
python mlstBLAST.py -s Kleborate/data/ybt_alleles.fasta -d Kleborate/data/YbST_profiles.txt *.fasta
````

## Typing direct from Illumina reads

MLST assignment can be achieved direct from reads using [SRST2](https://github.com/katholt/srst2). Simply 

* download the YbST and CbST allele sequences and profile tables from the /data directory in _this_ repository
* download the [SRST2](https://github.com/katholt/srst2) distribution; 
* Run SRST2, setting the --mlst_scheme and --mlst_definitions to point to the YbST or CbST allele sequences and profile tables like this:

```
srst2 --input_pe reads_1.fastq.gz reads_2.fastq.gz --output YbST --log --mlst_db YbST.fasta --mlst_definitions YbST.txt

srst2 --input_pe reads_1.fastq.gz reads_2.fastq.gz --output CbST --log --mlst_db CbST.fasta --mlst_definitions CbST.txt
```

Note that currently you can only run SRST2 with one MLST scheme at a time, so in order to type MLST, YbST and CbST you will need to run three separate commands:

```
srst2 --input_pe reads_1.fastq.gz reads_2.fastq.gz --output YbST --log --mlst_db YbST.fasta --mlst_definitions YbST.txt

srst2 --input_pe reads_1.fastq.gz reads_2.fastq.gz --output CbST --log --mlst_db CbST.fasta --mlst_definitions CbST.txt

srst2 --input_pe reads_1.fastq.gz reads_2.fastq.gz --output Klebs --log --mlst_db Klebsiella_pneumoniae.fasta --mlst_definitions kpnuemoniae.txt
```


## Other Klebs typing tools

Capsule typing: [Kaptive](https://github.com/katholt/Kaptive) (input = assemblies, command line)

MLST: [BIGSdb](http://bigsdb.pasteur.fr/klebsiella/klebsiella.html) (input = assemblies, web-based)

Resistance and plasmid loci: [SRST2](https://github.com/katholt/srst2) (input = reads, command line)


## Kleboration

Kleborate is under active development with many other Klebs genomic analysis tools and projects in progress. 

Please get in touch via the issues tracker if you have any issues, questions or ideas.

-------------

Stop! Kleborate and listen

ICEKp is back with with my brand new invention

If there was a problem, Klebs'll solve it

Check out the hook while Klebs evolves it
