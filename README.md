<p align="center"><img src="logo.png" alt="Kleborate" width="400"></p>

Kleborate is a tool to screen genome assemblies of _Klebsiella pneumoniae_ and the _Klebsiella pneumoniae_ species complex (KpSC) for:
 * MLST sequence type
 * species (e.g. _K. pneumoniae_, _K. quasipneumoniae_, _K. variicola_, etc.)
 * ICE<i>Kp</i> associated virulence loci: yersiniabactin (_ybt_), colibactin (_clb_)
 * virulence plasmid associated loci: salmochelin (_iro_), aerobactin (_iuc_), hypermucoidy (_rmpA_, _rmpA2_)
 * antimicrobial resistance genes, including quinolone resistance SNPs and colistin resistance truncations
 * K (capsule) and O antigen (LPS) serotype prediction, via _wzi_ alleles and [Kaptive](https://github.com/katholt/Kaptive)
 
For _Klebsiella_ outside of the KpSC, Kleborate will accurately determine the species and will report the presence of any accessory genes detected (AMR, virulence, K & O types); however species-focused markers (mutational resistance, MLST) will not be reported.

A manuscript describing the Kleborate software in full is currently in preparation. (Note that the BLAST logic has been checked in the light of [this article](https://doi.org/10.1093/bioinformatics/bty833) describing a common misconception regarding the BLAST parameter -max_target_seqs.)

In the meantime, if you use Kleborate, please cite the component schemes that you report:<br>
> Yersiniabactin and colibactin (ICE<i>Kp</i>) [Lam, MMC. et al. Genetic diversity, mobilisation and spread of the yersiniabactin-encoding mobile element ICE<i>Kp</i> in _Klebsiella pneumoniae_ populations. _Microbial Genomics_ (2018).](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000196)

> Aerobactin and salmochelin:
[Lam, MMC. et al. Tracking key virulence loci encoding aerobactin and salmochelin siderophore synthesis in _Klebsiella pneumoniae_. _Genome Medicine_ (2018).](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-018-0587-5)

> Kaptive for capsule (K) serotyping:
[Wyres, KL. et al. Identification of _Klebsiella_ capsule synthesis loci from whole genome data. _Microbial Genomics_ (2016).](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000102)

> Kaptive for O antigen (LPS) serotyping:
[Wick, RR et. al. Kaptive Web: user-friendly capsule and lipopolysaccharide serotype prediction for _Klebsiella_ genomes. _Journal of Clinical Microbiology_ (2018).](http://jcm.asm.org/content/56/6/e00197-18)


## Table of Contents

* [Background](#background)
* [Requirements](#requirements)
* [Installation](#installation)
   * [Update the MLST database](#updating-the-mlst-database)
* [Basic usage](#basic-usage)
* [Full usage](#full-usage)
* [Screening functions & outputs](#screening-details)
   * [Genome assembly quality metrics](#assembly-quality-metrics)
   * [<em>Klebsiella</em> species identification](#klebsiella-species)
   * [Multi-locus sequence typing (MLST)](#mlst)
   * [Acquired virulence loci](#acquired-virulence-loci)
   * [Resistance gene detection](#antimicrobial-resistance-determinants)
   * [Resistance & virulence scores and counts](#scores-and-counts)
   * [Serotype prediction](#serotype-prediction)
* [Example output](#example-output)
   * [Test commands](#test-commands)
   * [Concise results (stdout)](#concise-results-stdout)
   * [Full results (file)](#full-results-file)
* [Using Kleborate genotyping schemes with raw Illumina reads](#typing-from-illumina-reads)
* [Contact us](#contact-us)
* [License](#license)



## Background

_Klebsiella pneumoniae_ (_Kp_) is a commensal bacterium that causes opportunistic infections in hospitals. _K. pneumoniae_ are intrinsically resistant to ampicillin and frequently acquire additional antimicrobial resistances through horizontal gene transfer and chromosomal mutations. A handful of hypervirulent lineages are also recognised, which encode a constellation of acquired virulence factors and can cause invasive disease outside the hospital setting. Evidence is now mounting that other _K. pneumoniae_ strains carrying one or more of these acquired factors – including siderophores (yersiniabactin, salmochelin and aerobactin), regulators of hypermucoidy (_rmpA/rmpA2_ genes) and/or the genotoxin colibactin – can also be highly pathogenic and cause more severe disease both inside and outside hospitals. Capsule (K) and LPS (O) antigen variation is also of great interest to the research community due to its importance in host-pathogen and phage interactions, and thus potential relevance to alternative control measures such as vaccines, immunotherapy and phage therapy. _K. pneumoniae_ has six close relatives (species and subspecies) known as the [_K. pneumoniae_ species complex (KpSC)](#klebsiella-species); these are difficult to distinguish from one another in clinical labs using biotyping or MALDI-TOF and are often confused for _K. pneumoniae sensu stricto_.

To make it easier to extract clinically relevant genotyping information on _K. pneumoniae_ and the species complex suing genome data, we have developed *Kleborate*, a genomic surveillance tool designed to (a) accurately identify species and sequence types, and (b) identify the key *acquired* genetic features for which there is strong evidence of association with either antibiotic resistance or hypervirulence in _K. pneumoniae sensu stricto_. While many generic tools can be used to identify sequence types or resistance determinants from bacterial genomes, we hope that this organism-specific tool will help avoid many of the common confusions faced by people working with _K. pneumoniae_ genomes and also facilitate monitoring for the convergence of antibiotic resistance with the hypervirulence-associated factors noted above. 


## Requirements

Software requirements:
* [Python 3](https://www.python.org/) (v3.5 or later should work)
* [setuptools](https://pypi.python.org/pypi/setuptools) (required to install Kleborate)
  * To install: `pip install setuptools`
* [Biopython](https://biopython.org/)
  * To install: `pip install biopython`
* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
  * Version 2.7.1 or later is needed, as earlier versions have a bug with the `culling_limit` parameter and/or tblastx results.
  * We test Kleborate on BLAST+ v2.7.1. Later versions will probably also work but stick to v2.7.1 if you want to play it safe.
* [Mash](https://github.com/marbl/Mash)
  * You can download a pre-compiled version from the [Mash releases page](https://github.com/marbl/Mash/releases) (both Mac and Linux binaries are available) and copy the executable somewhere into your PATH (e.g. `/usr/local/bin`).
  * Alternatively, you can install Mash on a Mac with Homebrew: `brew install mash`


Input files:
Kleborate takes _Klebsiella_ genome assemblies (either completed or draft) in fasta format (can be gzipped). If you have unassembled reads, try assembling them with our [Unicycler](https://github.com/rrwick/Unicycler) assembler which works great on Illumina or hybrid Illumina + Nanopore/PacBio reads).



## Installation

Kleborate can be installed to your system for easy usage:

```bash
git clone --recursive https://github.com/katholt/Kleborate.git
cd Kleborate
python setup.py install
kleborate -h
```


Alternatively, you can clone and run Kleborate without installation directly from its source directory:

```bash
git clone --recursive https://github.com/katholt/Kleborate.git
Kleborate/kleborate-runner.py -h
```

See [examples below](#example-output) to test out your installation on some public genome data. And if you'd like to thoroughly check that everything works as intended, you can also run this repo's [automated tests](test) after installation.

Note that Kleborate depends on a git submodule ([Kaptive](https://github.com/katholt/Kaptive)) which is why `--recursive` is required when cloning. If you update your local copy of Kleborate using `git pull`, you should also run `git submodule update` to ensure that its Kaptive is also up-to-date.


### Updating the MLST database

Each Kleborate release includes a copy of the [_K. pneumoniae_ species complex MLST database](http://bigsdb.pasteur.fr/klebsiella/klebsiella.html) to screen against. The version included is current at the time of the release, however the _K. pneumoniae_ species complex BIGSdb is being updated all the time with new STs, so Kleborate users may wish to update their copy of Kleborate regularly with the latest MLST database.

The MLST database is made up of 2 files, which are located in the `Kleborate/kleborate/data` directory:
* `Klebsiella_pneumoniae.fasta` (allele seuqences)
* `kpneumoniae.txt` (sequence type definitions)

A python3 script to download the latest versions of these 2 files is provided in the `Kleborate/scripts` directory, the downloaded files can then just be copied into the `Kleborate/kleborate/data` directory

```
cd Kleborate/scripts
python getmlst.py --species "Klebsiella pneumoniae"
mv Klebsiella_pneumoniae.fasta ../kleborate/data
mv kpneumoniae.txt ../kleborate/data
```


## Basic usage

__Screen some genomes for MLST and virulence loci:__<br>
`kleborate -o results.txt -a *.fasta`

__Also screen for resistance genes:__<br>
`kleborate --resistance -o results.txt -a *.fasta`

__Turn on all of Kleborate's optional screens (resistance genes, species check and both K and O loci):__<br>
`kleborate --all -o results.txt -a *.fasta`

__Screen everything in a set of gzipped assemblies:__<br>
`kleborate --all -o results.txt -a *.fasta.gz`



## Full usage

```
usage: kleborate -a ASSEMBLIES [ASSEMBLIES ...] [-r] [-s] [--kaptive_k]
                 [--kaptive_o] [-k] [--all] [-o OUTFILE]
                 [--kaptive_k_outfile KAPTIVE_K_OUTFILE]
                 [--kaptive_o_outfile KAPTIVE_O_OUTFILE] [-h] [--version]

Kleborate: a tool for characterising virulence and resistance in Klebsiella

Required arguments:
  -a ASSEMBLIES [ASSEMBLIES ...], --assemblies ASSEMBLIES [ASSEMBLIES ...]
                        FASTA file(s) for assemblies, can be gzipped (.gz)

Screening options:
  -r, --resistance      Turn on resistance genes screening (default: no
                        resistance gene screening)
  --kaptive_k           Turn on Kaptive screening of K loci (default: do not
                        run Kaptive for K loci)
  --kaptive_o           Turn on Kaptive screening of O loci (default: do not
                        run Kaptive for O loci)
  -k, --kaptive         Equivalent to --kaptive_k --kaptive_o
  --all                 Equivalent to --resistance --kaptive

Output options:
  -o OUTFILE, --outfile OUTFILE
                        File for detailed output (default:
                        Kleborate_results.txt)
  --kaptive_k_outfile KAPTIVE_K_OUTFILE
                        File for full Kaptive K locus output (default: do not
                        save Kaptive K locus results to separate file)
  --kaptive_o_outfile KAPTIVE_O_OUTFILE
                        File for full Kaptive O locus output (default: do not
                        save Kaptive O locus results to separate file)

Help:
  -h, --help            Show this help message and exit
  --version             Show program's version number and exit
```



## Screening details

### Assembly quality metrics

The quality and completeness of Kleborate results depends on the quality of the input genome assemblies. We provide some basic assembly statistics (contig count,	N50,	largest contig size, detection of ambiguous bases) to help users understand their Kleborate results in the context of assembly quality, but we recommend users conduct more comprehensive QC themselves before running Kleborate (e.g. screen for contamination, etc).


### _Klebsiella_ species

Kleborate will attempt to identify the species of each input assembly. It does this by comparing the assembly using Mash to a curated set of _Klebsiella_ assemblies [from NCBI](https://www.ncbi.nlm.nih.gov/assembly) and reporting the species of the closest match. Kleborate considers a Mash distance of ≤ 0.01 to be a strong species match. A distance of > 0.01 and ≤ 0.03 is a weak match and might indicate that your sample is a novel lineage or a hybrid between multiple _Klebsiella_ species.

Here is an annotated tree of the reference assemblies, made by [mashtree](https://github.com/lskatz/mashtree):
<p align="center"><img src="images/species_tree.png" alt="Klebsiella species tree" width="90%"></p>

Kleborate is designed for detailed genotyping of the well-studied _K. pneumoniae_ species complex (KpSC) labelled on the tree, which includes the seven species listed in the table below. These were previously considered as phylogroups within _K. pneumoniae_. We've included the phylogroup numbers in the table below to allow backwards compatibility, but these are not reported in the Kleborate output. 

| Species                                       | Kp phylogroup<sup>a</sup> | Kp phylogroup (alternative)<sup>b</sup> | Reference |
| --------------------------------------------- | ---------------------- | -------------------------------- | --------- |
| _K. pneumoniae_                               | Kp1                    | KpI                              | [Brenner, D.J. 1979 Int J Syst Evol Microbiol 29: 38-41](https://ijs.microbiologyresearch.org/content/journal/ijsem/10.1099/00207713-29-1-38) | 
| _K. quasipneumoniae_ subsp _quasipneumoniae_    | Kp2                    | KpIIa                            | [Brisse et al. 2014 Int J Syst Evol Microbiol 64:3146-52](https://ijs.microbiologyresearch.org/content/journal/ijsem/10.1099/ijs.0.062737-0#tab2) | 
| _K. quasipneumoniae_ subsp _similipneumoniae_   | Kp4                    | KpIIb                            | [Brisse et al. 2014 Int J Syst Evol Microbiol 64:3146-52](https://ijs.microbiologyresearch.org/content/journal/ijsem/10.1099/ijs.0.062737-0#tab2) | 
| _K. variicola_ subsp _variicola_                | Kp3                    | KpIII                            | [Rosenblueth et al. 2004 Syst Appl Microbiol 27:27-35](https://www.sciencedirect.com/science/article/abs/pii/S0723202004702349?via%3Dihub) (described as subsp _variicola_ in this paper) |
| _K. variicola_ subsp _tropica_            | Kp5                    | -                                | [Rodrigues et al. 2019 Res Microbiol ﻿S0923-2508:﻿30019-1](https://www.sciencedirect.com/science/article/pii/S0923250819300191?via%3Dihub) (described as subsp _tropicalensis_ in this paper) | 
| _K. quasivariicola_                           | Kp6                    | -                                | [Long et al. 2017 Genome Announc 5: ﻿e01057-17](https://mra.asm.org/content/5/42/e01057-17) | 
| _K. africana_                             | Kp7                    | -                                | [Rodrigues et al. 2019 Res Microbiol ﻿S0923-2508:﻿30019-1](https://www.sciencedirect.com/science/article/pii/S0923250819300191?via%3Dihub) (described as _africanensis_ in this paper) | 

<sup>a</sup> Kp phylogroup numbers as described in [Rodrigues et al. 2019](https://www.sciencedirect.com/science/article/pii/S0923250819300191?via%3Dihub)

<sup>b</sup> alternative (older) Kp phylogroup numbers as described in [Brisse et al. 2001](https://ijs.microbiologyresearch.org/content/journal/ijsem/10.1099/00207713-51-3-915#tab2) and [Fevre et al. 2005](https://aac.asm.org/content/49/12/5149) prior to the identification of _K. variicola_ subsp _tropica_, _K. quasivariicola_ and _K. africana_.

More distant _Klebsiella_ species (_oxytoca_, _michiganensis_, _grimontii_ and _aerogenes_) will be accurately identified by Kleborate, although please note that the diversty and relevance of _K. pneumoniae_ virulence factors in these species is not yet well understood.

Kleborate will also yield reliable species identifications across the family Enterobacteriaceae, as different species sometimes end up in _Klebsiella_ collections. These names are again assigned based on the clades in a mashtree, but were not as carefully curated as the _Klebsiella_ species (so take them with a grain of salt).



### MLST

Genomes identified by Kleborate as belonging to the _K. pneumoniae_ species complex are then subjected to multi-locus sequence typing (MLST) using the 7-locus scheme described at the [_K. pneumoniae_ BIGSdb hosted at the Pasteur Institute](http://bigsdb.pasteur.fr/klebsiella/klebsiella.html). Note that this scheme is not specific to _K. pneumoniae sensu stricto_ but covers the whole [_K. pneumoniae_ species complex](#klebsiella-species). A copy of the MLST alleles and ST definitions is stored in the [data directory](https://github.com/katholt/Kleborate/tree/master/kleborate/data) of this repository. See [above](#updating-the-mlst-database) for instructions on how to update the MLST database in your copy of Kleborate.

Notes on Kleborate's MLST calls:
* Kleborate makes an effort to report the closest matching ST if a precise match is not found.
* Imprecise allele matches are indicated with a `*`.
* Imprecise ST calls are indicated with `-nLV`, where n indicates the number of loci that disagree with the ST reported. So `258-1LV` indicates a single-locus variant of (SLV) of ST258, i.e. 6/7 loci match ST258.

*Note that allele definitions for ST1047 and ST1078 were changed in the MLST database in Feburary 2018, and these new allele combinations are incorporated in Kleborate since v0.4.0. This is highly unusual and other allele and ST assignment should be stable across versions.*

|allele         |ST1047 old     |ST1047 current |ST1078 old     |ST1078 current |
| ------------- | ------------- |---------------|---------------|---------------|
|_gapA_         |10             |2              |16             |4              |
|_infB_         |20             |1              |18             |5              |
|_mdh_          |1              |2              |1              |1              |
|_pgi_          |1              |20             |76             |3              |
|_phoE_         |9              |7              |47             |12             |
|_rpoB_         |11             |1              |1              |4              |
|_tonB_         |14             |4              |124            |46             |


### Acquired virulence loci

Kleborate examines four key acquired virulence gene clusters that contribute to hypervirulence in _K. pneumoniae_: the siderophores yersiniabactin (_ybt_), aerobactin (_iuc_) and salmochelin (_iro_), and the genotoxin colibactin (_clb_). (We also screen for the hypermucoidy genes _rmpA_ and _rmpA2_, details below).
* For each of these loci, Kleborate will call a sequence type using the same logic as the MLST described above, using the locus-specific schemes defined in the [BIGSdb](http://bigsdb.pasteur.fr/klebsiella/klebsiella.html).
* Kleborate will also report the lineage associated with the virulence sequence types, as outlined below and detailed in the corresponding papers (for yersiniabactin, we also report the predicted ICE<i>Kp</i> structure based on the _ybt_ lineage assignment).
* If the locus is not detected, Kleborate reports the ST as `0` and the lineage as `-`.

#### Yersiniabactin and colibactin (primarily mobilised by ICE<i>Kp</i>)

We recently explored the diversity of the _K. pneumoniae_ integrative conjugative element (ICE<i>Kp</i>), which mobilises the yersiniabactin locus _ybt_, using genomic analysis of a diverse set of 2498 _Klebsiella_ (see [this paper](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000196)). Overall, we found _ybt_ in about a third of all _K. pneumoniae_ genomes and _clb_ in about 14%. We identified 17 distinct lineages of _ybt_ (see figure) embedded within 14 structural variants of ICE<i>Kp</i> that can integrate at any of four tRNA-Asn sites in the chromosome. Three of the 17 _ybt_ lineages were associated with three lineages of colibactin, with which they are co-located in the same ICE structure designated ICE<i>Kp10</i>. One ICE structure (ICE<i>Kp1</i>) carries the salmochelin synthesis locus _iro_ and _rmpA_ hypermucoidy gene in addition to _ybt_ (lineage 2). Additionally, we identified a lineage of _ybt_ that is plasmid-encoded, representing a new mechanism for _ybt_ dispersal in _K. pneumoniae_ populations. Based on this analysis, we developed a MLST-style approach for assigning yersiniabactin sequence types (YbST) and colibactin sequence types (CbST), which is implemented in Kleborate. Annotated reference sequences for each ICE<i>Kp</i> variant are included in the [data directory](https://github.com/katholt/Kleborate/tree/master/kleborate/data) of this repository). 

ICE<i>Kp</i> is occasionally found in other species within the KpSC, and even in other genera of Enterobacteriaceae (see [paper](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000196)), however most of the known variation included in the database is derived from _K. pneumoniae_.

<p align="left"><img src="images/ybt_trees.png" alt="ybt tree" width="70%"></p>

#### Aerobactin and salmochelin (primarily mobilised by virulence plasmids)
We further explored the genetic diversity of the aerobactin (_iuc_) and salmochelin (_iro_) loci among a dataset of 2733 _Klebsiella_ genomes (see [this paper](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-018-0587-5)). We identified five _iro_ and six _iuc_ lineages (see figure), each of which was associated with a specific location within _K. pneumoniae_ genomes. The most common lineages were _iuc1_ and _iro1_, which are found together on the virulence plasmid KpVP-1 (typified by pK2044 or pLVPK common to the hypervirulent clones ST23, ST86, etc). _iuc2_ and _iro2_ lineages were associated with the alternative virulence plasmid KpVP-2 (typified by Kp52.145 plasmid II from the K2 ST66 lab strain known as Kp52.145 or B5055). _iuc5_ and _iro5_ originate from _E. coli_ and are carried (often together) on _E. coli_ plasmids that can transfer to _K. pneumoniae_. The lineages _iuc2A_, _iuc3_ and _iro4_ were associated with other novel plasmids that had not been previously described in _K. pneumoniae_ but sequences for which are included in [the paper](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-018-0587-5). In addition, we found the salmochelin locus present in ICE<i>Kp1</i> constitutes its own lineage _iro3_, and the aerobactin locus present in the chromosome of ST67 _K. pneumoniae_ subsp _rhinoscleromatis_ strains constitutes its own lineage _iuc4_. Based on this analysis, we developed a MLST-style approach for assigning aerobactin sequence types (AbST) and salmochelin sequence types (SmST) which is implemented in Kleborate.

<p align="center"><img src="images/iuc_iro_trees.png" alt="iuc and iro trees" width="70%"></p>

Please note that the aerobactin _iuc_ and salmochelin _iro_ lineage names have been updated between Kleborate version 0.2.0 and 0.3.0 to match the nomenclature used in [the paper](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-018-0587-5). The AbST and SmST allele numbers are unchanged. Lineage name re-assignments are:

| v0.2.0        | v0.3.0        | location (see [paper](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-018-0587-5) for details)
| ------------- | ------------- |------------------
| *iuc 2*         | *iuc 1*         | KpVP-1 (e.g. pLVPK) |
| *iuc 3B*        | *iuc 2*         | KpVP-2 (e.g. Kp52.145 plasmid II) |
| *iuc 3A*        | *iuc 2A*        | other plasmids |
| *iuc 4*         | *iuc 3*         | other plasmids |
| *iuc 5*         | *iuc 4*         | rhinoscleromatis chromosome |
| *iuc 1*         | *iuc 5*         | _E. coli_ variant |
| *iro 3*         | *iro 1*         | KpVP-1 (e.g. pLVPK) |
| *iro 4*         | *iro 2*         | KpVP-2 |
| *iro 5*         | *iro 3*         | ICEKp1 |
| *iro 2*         | *iro 4*         | _Enterobacter_ variant |
| *iro 1*         | *iro 5*         | _E. coli_ variant |



#### Hypermucoidy genes

Kleborate screens for alleles of the _rmpA_ and _rmpA2_ genes which can result in a hypermucoid phenotype by upregulating capsule production.

* The two genes share ~83% nucleotide identity so are easily distinguished, and are reported in separate columns.
* Alleles for each gene are sourced from the [BIGSdb](http://bigsdb.pasteur.fr/klebsiella/klebsiella.html). For _rmpA_, we have also mapped these alleles to the various known locations for _rmpA_ in _Klebsiella_ (i.e. major virulence plasmids KpVP-1 and KpVP-2; other virulences plasmids simply designated as VP; ICE<i>Kp1</i> and the chromosome in rhinoscleromatis).
* Unique (non-overlapping) nucleotide BLAST hits with >95% identity and >50% coverage are reported. Note multiple hits to the same gene are reported if found (e.g. the NTUH-K2044 genome carries _rmpA_ in the virulence plasmid and also in ICE<i>Kp1</i>, which is reported in the _rmpA_ column as rmpA_11(ICEKp1),rmpA_2(KpVP-1)).
* Truncations in the _rmpA_ and _rmpA2_ genes are expressed as a percentage of the amino acid length from the start codon, e.g. rmpA_5-54% indicates the RmpA protein is truncated after 54% length of the intact amino acid sequence. These truncations appear to be common, due to insertions and deletions within a poly-G tract, and almost certainly result in loss of protein function.


### Antimicrobial resistance determinants

By using the `--resistance` option, Kleborate will screen for acquired resistance genes against the ARG-Annot database of acquired resistance genes (updated version from [SRST2](https://github.com/katholt/srst2)), which includes allelic variants. It attempts to report the best matching variant for each locus in the genome:
* Exact nucleotide matches are reported with no further annotation (e.g. "TEM-15"). 
* If no exact nucleotide match is found, Kleborate searches for an exact amino acid match, and will report this with a "^" symbol (e.g. "TEM-15^" indicates an exact match to the TEM-15 protein sequence but with 1 or more nucleotide differences). If no exact amino acid match is found, the closest nucleotide match is reported with "\*" symbol (e.g. "TEM-30\*" indicates no precise nucleotide or amino acid match is found, but the closest nucleotide match is to TEM-30).
* If the length of match is less than the length of the reported allele (i.e. a partial match), this is indicated with `?`.
* Note that KpSC carry a core beta-lactamase gene (SHV in _K. pneumoniae_, LEN in _K. variicola_, OKP in _K. quasipneumoniae_) that confers clinically significant resistance to ampicillin. As these are present in all genomes, non-ESBL alleles of these genes are not included in the count of acquired resistance genes or drug classes.
* ESBL alleles of SHV are almost always carried on plasmids (in addition to the intrinsic narrow-spectrum SHV/LEN/OKP allele in the chromosome). However it is possible to have a mutation in a chromosomal SHV gene that gives a match to an ESBL allele, which would be reported in the ESBL column and counted as an acquired gene (and it is very hard to tell the difference without manual exploration of the genetic context).
  * See [this paper](http://www.pnas.org/content/112/27/E3574.long) for more information.
* Note that _oqxAB_ and _fosA_ are also core genes in _K. pneumoniae_ and don't confer clinical resistance to fluoroquinolones or fosfomycin, hence Kleborate does not report them.

Using the `--resistance` option also turns on screening for chromosomal mutations for which there is strong evidence of an association with clinical resistance in KpSC (note these are ONLY reported if the genome was recognised as part of the KpSC):
* Fluoroquinolone resistance SNPs: GyrA 83 & 87 and ParC 80 & 84.
* Colistin resistance due to truncation or loss of MgrB or PmrB (truncations are expressed as % amino acid length from the start codon).
* OmpK35 and OmpK36 truncations and mutations resulting in reduced susceptibility to beta-lactamases (truncations are expressed as % amino acid length from the start codon). See [this paper](https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1007218) for more information.

Note these do not count towards acquired resistance gene counts, but do count towards drug classes (with the exception of Omp mutations, whose spectrum of effects depends on the presence of acquired beta-lactamases and thus their impact on specific beta-lactam drug classes is hard to predict).

All resistance results (both for the gene screen and mutation screen) are grouped by drug class (according to the [ARG-Annot](https://www.ncbi.nlm.nih.gov/pubmed/24145532) DB), with beta-lactamases broken down into [Lahey](https://www.lahey.org/Studies/) classes, as follows: 
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
* Ntmdz (nitroimidazole, e.g. metronidazole)
* Phe (phenicols)
* Rif (rifampin)
* Sul (sulfonamides)
* Tet (tetracyclines)
* Tmt (trimethoprim)
* Tgc (tigecycline)

Note there is a separate column 'Omp' reporting known resistance-related mutations in the OmpK35 and OmpK36 osmoporins. 

Note that Kleborate reports resistance results for all antimicrobial classes with confidently attributable resistance mechanisms in KpSC. Not all of these are actually used clinically for treatment of KpSC infections (e.g. Ntmdz, MLS, Rif) but they are still reported here as the presence of acquired resistance determinants to these classes is of interest to researchers for other reasons (e.g. these genes can be useful markers of MGEs and MGE spread; there is potential for use of these drugs against other organisms to select for KpSC in co-infected patients or in the environment). For an overview of antimicrobial resistance and consensus definitions of multidrug resistance (MDR), extreme drug resistance (XDR) and pan drug resistance in Enterobacteriaceae, see [Magiorakos 2012](https://www.clinicalmicrobiologyandinfection.com/article/S1198-743X(14)61632-3/fulltext).



### Scores and counts

Kleborate outputs a simple categorical virulence score, and if resistance screening is enabled, an antimicrobial resistance score as well. These scores provide a rough categorisation of the strains to facilitate monitoring resistance-virulence convergence:

* The virulence score ranges from 0 to 5:
  * 0 = none of the acquired virulence loci (i.e. negative for all of yersiniabactin, colibactin, aerobactin, salmochelin)
  * 1 = yersiniabactin only
  * 2 = yersiniabactin and colibactin, or colibactin only 
  * 3 = aerobactin and/or salmochelin only (without yersiniabactin or colibactin)
  * 4 = aerobactin and/or salmochelin with yersiniabactin (without colibactin)
  * 5 = yersiniabactin, colibactin and aerobactin and/or salmochelin
  
* The resistance score ranges from 0 to 3:
  * 0 = no ESBL, no carbapenemase (regardless of colistin resistance)
  * 1 = ESBL, no carbapenemase (regardless of colistin resistance)
  * 2 = Carbapenemase without colistin resistance (regardless of ESBL, OmpK mutations not considered)
  * 3 = Carbapenemase with colistin resistance (regardless of ESBL, OmpK mutations not considered)

When resistance screening is enabled, Kleborate also quantifies how many acquired resistance genes are present and how many drug classes (in _addition_ to Bla/ampicillin) have at least one resistance determinant detected. A few things to note:
  * The presence of resistance _mutations_, and non-ESBL forms of core genes SHV/LEN/OKP, do not contribute to the resistance _gene_ count.
  * Mutations do contribute to the drug class count, e.g. fluoroquinolone resistance will be counted if a GyrA mutation is encountered regardless of whether or not an acquired fluoroquinolone resistance is also present. The exception is Omp mutations, which do not contribute to the drug class count as their effect depends on the strain background and the presence of acquired beta-lactamase enzymes; hence this information is provided in a separate column, and interpretation is left to the user.
  * Note that since a drug class can have multiple resistance determinants, the gene count is typically higher than the class count.




### Serotype prediction

#### Basic capsule prediction with _wzi_ allele typing
By default, Kleborate will report the closest match amongst the _wzi_ alleles in the [BIGSdb](http://bigsdb.pasteur.fr/klebsiella/klebsiella.html). This is a marker of capsule locus (KL) type, which is highly predictive of capsule (K) serotype. Although there is not a 1-1 relationship between wzi allele and KL/K type, there is a strong correlation (see [Wyres et al, MGen 2016](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000102) and [Brisse et al, J Clin Micro 2013](https://jcm.asm.org/content/51/12/4073.long)). Note the _wzi database_ is populated with alleles from the _Klebsiella pneumoniae_ species complex and is not reliable for other species.

The _wzi_ allele can provide a handy way of spotting the virulence-associated types (wzi=K1, wzi2=K2, wzi5=K5); or spotting capsule switching within clones, e.g. you can tell which ST258 lineage you have from the wzi type (wzi154: the main lineage II; wzi29: recombinant lineage I; others: probably other recombinant lineages).

#### Capsule (K) and O antigen (LPS) serotype prediction using Kaptive
You can optionally turn on capsule and O antigen typing using the dedicated capsule typing tool [Kaptive](https://github.com/katholt/Kaptive). Note that the Kaptive database comprises O and K loci characterised in the _Klebsiella pneumoniae_ species complex (see [Wyres et al, MGen 2016](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000102) for K loci, [Wick et al, J Clin Micro 2018](https://jcm.asm.org/content/56/6/e00197-18) for O loci). These loci are sometimes also found in other _Klebsiella_ species but you should expect many novel loci outside the KpSC that will not be detected here.

* `--kaptive_k` turns on Kaptive screening of the K locus
* `--kaptive_o` turns on Kaptive screening of the O locus
* `--kaptive` turns on both (is equivalent to `--kaptive_k --kaptive_o`)

Note that running Kaptive will significantly increase the runtime of Kleborate, but provide much more detailed information about the K and/or O loci and their genes.

If Kaptive is switched on the Kleborate report will include a column to indicate the Kaptive confidence match for the reported best-matching K or O locus (see [here](https://github.com/katholt/Kaptive/blob/master/README.md#output-files) for a description of the logic). We recommend reporting only K and O loci with a confidence level of "Good" or better. Calls with confidence level "Low" or "None" should be considered carefully as they may result from assembly problems (fragmentation) or novel variation in the K/O locus. 
If you think you have found a novel K or O locus and would like us to add it to the Kaptive database please get in touch.

## Example output

### Test commands

Run these commands to test out Kleborate using some of the test data provided in the /test directory of this repository:

```
# 1) basic genotyping (no resistance typing; K serotype prediction using wzi allele only)
kleborate -o results.txt -a Kleborate/test/sequences/GCF_002248955.1.fna.gz Kleborate/test/sequences/GCF_003095495.1.fna.gz Kleborate/test/sequences/GCF_000009885.1.fna.gz Kleborate/test/sequences/GCF_900501255.1.fna.gz Kleborate/test/sequences/GCF_000019565.1.fna.gz Kleborate/test/sequences/GCF_000492415.1.fna.gz Kleborate/test/sequences/GCF_000492795.1.fna.gz

# 2) with resistance typing (K serotype prediction using wzi allele only)
kleborate -o results_res.txt --resistance -a Kleborate/test/sequences/GCF_002248955.1.fna.gz Kleborate/test/sequences/GCF_003095495.1.fna.gz Kleborate/test/sequences/GCF_000009885.1.fna.gz Kleborate/test/sequences/GCF_900501255.1.fna.gz Kleborate/test/sequences/GCF_000019565.1.fna.gz Kleborate/test/sequences/GCF_000492415.1.fna.gz Kleborate/test/sequences/GCF_000492795.1.fna.gz

# 3) with resistance typing & full K/O serotype prediction using Kaptive (slower)
kleborate -o results_res_kaptive.txt --all -a Kleborate/test/sequences/GCF_002248955.1.fna.gz Kleborate/test/sequences/GCF_003095495.1.fna.gz Kleborate/test/sequences/GCF_000009885.1.fna.gz Kleborate/test/sequences/GCF_900501255.1.fna.gz Kleborate/test/sequences/GCF_000019565.1.fna.gz Kleborate/test/sequences/GCF_000492415.1.fna.gz Kleborate/test/sequences/GCF_000492795.1.fna.gz
```


### Concise results (stdout)

These are the concise Kleborate results that are printed to the terminal, for example 1:

strain | species | ST | virulence_score | Yersiniabactin | YbST | Colibactin | CbST | Aerobactin | AbST | Salmochelin | SmST | rmpA | rmpA2 | wzi | K_locus
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---
GCF_002248955.1 | Klebsiella pneumoniae | ST15 | 0 | - | 0 | - | 0 | - | 0 | - | 0 | - | - | wzi29 | KL106
GCF_003095495.1 | Klebsiella pneumoniae | ST258 | 0 | - | 0 | - | 0 | - | 0 | - | 0 | - | - | wzi154 | KL107
GCF_000009885.1 | Klebsiella pneumoniae | ST23 | 4 | ybt 2; ICEKp1 | 326 | - | 0 | iuc 1 | 1 | iro 3 | 18-1LV | rmpA_11(ICEKp1),rmpA_2(KpVP-1) | rmpA2_3-47% | wzi1 | KL1
GCF_900501255.1 | Klebsiella pneumoniae | ST86 | 3 | - | 0 | - | 0 | iuc 1 | 1 | iro 1 | 1 | rmpA_2(KpVP-1) | rmpA2_4*-50% | wzi2 | KL2 (KL30)
GCF_000019565.1 | Klebsiella variicola subsp. variicola | ST146 | 0 | - | 0 | - | 0 | - | 0 | - | 0 | - | - | wzi159 | KL30
GCF_000492415.1 | Klebsiella quasipneumoniae subsp. quasipneumoniae | ST1437 | 0 | - | 0 | - | 0 | - | 0 | - | 0 | - | - | wzi185 | KL46
GCF_000492795.1 | Klebsiella quasipneumoniae subsp. similipneumoniae | ST1435 | 0 | - | 0 | - | 0 | - | 0 | - | 0 | - | - | wzi183 | KL21

For example 2 (ie with resistance typing turned on):

strain | species | ST | virulence_score | resistance_score | Yersiniabactin | YbST | Colibactin | CbST | Aerobactin | AbST | Salmochelin | SmST | rmpA | rmpA2 | wzi | K_locus | AGly | Col | Fcyn | Flq | Gly | MLS | Ntmdz | Phe | Rif | Sul | Tet | Tgc | Tmt | Omp | Bla | Bla_Carb | Bla_ESBL | Bla_ESBL_inhR | Bla_broad | Bla_broad_inhR
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---
GCF_002248955.1 | Klebsiella pneumoniae | ST15 | 0 | 0 | - | 0 | - | 0 | - | 0 | - | 0 | - | - | wzi29 | KL106 | Aac3-IId^ | Mcr3-1* | - | GyrA-83F;GyrA-87A;ParC-80I | - | - | - | CatA1^ | - | - | TetA | - | - | - | SHV-28^ | - | - | - | - | -
GCF_003095495.1 | Klebsiella pneumoniae | ST258 | 0 | 3 | - | 0 | - | 0 | - | 0 | - | 0 | - | - | wzi154 | KL107 | Aac3-IId^;AadA2^;Aph3-Ia^;RmtB;Sat-2A;StrA^;StrB | MgrB-62%;PmrB-36% | - | GyrA-83I;ParC-80I | - | Erm42*;MphA | - | CatA1^ | - | SulI;SulII | TetG | - | DfrA12? | OmpK35-25%;OmpK36GD | TEM-1D^ | KPC-2 | CTX-M-14 | - | SHV-11 | -
GCF_000009885.1 | Klebsiella pneumoniae | ST23 | 4 | 0 | ybt 2; ICEKp1 | 326 | - | 0 | iuc 1 | 1 | iro 3 | 18-1LV | rmpA_11(ICEKp1),rmpA_2(KpVP-1) | rmpA2_3-47% | wzi1 | KL1 | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | SHV-11^ | -
GCF_900501255.1 | Klebsiella pneumoniae | ST86 | 3 | 0 | - | 0 | - | 0 | iuc 1 | 1 | iro 1 | 1 | rmpA_2(KpVP-1) | rmpA2_4*-50% | wzi2 | KL2 (KL30) | - | - | - | - | - | - | - | - | - | - | - | - | - | - | SHV-187* | - | - | - | - | -
GCF_000019565.1 | Klebsiella variicola subsp. variicola | ST146 | 0 | 0 | - | 0 | - | 0 | - | 0 | - | 0 | - | - | wzi159 | KL30 | - | - | - | - | - | - | - | - | - | - | - | - | - | - | LEN-24*;LEN-24* | - | - | - | - | -
GCF_000492415.1 | Klebsiella quasipneumoniae subsp. quasipneumoniae | ST1437 | 0 | 0 | - | 0 | - | 0 | - | 0 | - | 0 | - | - | wzi185 | KL46 | Aac6-Ib;StrA*;StrB* | - | - | - | - | - | - | CatA2* | - | SulII | - | - | DfrA14 | - | - | - | - | - | OKP-A-3* | -
GCF_000492795.1 | Klebsiella quasipneumoniae subsp. similipneumoniae | ST1435 | 0 | 0 | - | 0 | - | 0 | - | 0 | - | 0 | - | - | wzi183 | KL21 | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | OKP-B-7* | -

For example 3 (ie with resistance typing & Kaptive serotyping turned on):

strain | species | ST | virulence_score | resistance_score | Yersiniabactin | YbST | Colibactin | CbST | Aerobactin | AbST | Salmochelin | SmST | rmpA | rmpA2 | wzi | K_locus | K_locus_confidence | O_locus | O_locus_confidence | AGly | Col | Fcyn | Flq | Gly | MLS | Ntmdz | Phe | Rif | Sul | Tet | Tgc | Tmt | Omp | Bla | Bla_Carb | Bla_ESBL | Bla_ESBL_inhR | Bla_broad | Bla_broad_inhR
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---
GCF_002248955.1 | Klebsiella pneumoniae | ST15 | 0 | 0 | - | 0 | - | 0 | - | 0 | - | 0 | - | - | wzi29 | KL107 | None | O1/O2v2 | Very high | Aac3-IId^ | Mcr3-1* | - | GyrA-83F;GyrA-87A;ParC-80I | - | - | - | CatA1^ | - | - | TetA | - | - | - | SHV-28^ | - | - | - | - | -
GCF_003095495.1 | Klebsiella pneumoniae | ST258 | 0 | 3 | - | 0 | - | 0 | - | 0 | - | 0 | - | - | wzi154 | KL107 | Good | O2v2 | Good | Aac3-IId^;AadA2^;Aph3-Ia^;RmtB;Sat-2A;StrA^;StrB | MgrB-62%;PmrB-36% | - | GyrA-83I;ParC-80I | - | Erm42*;MphA | - | CatA1^ | - | SulI;SulII | TetG | - | DfrA12? | OmpK35-25%;OmpK36GD | TEM-1D^ | KPC-2 | CTX-M-14 | - | SHV-11 | -
GCF_000009885.1 | Klebsiella pneumoniae | ST23 | 4 | 0 | ybt 2; ICEKp1 | 326 | - | 0 | iuc 1 | 1 | iro 3 | 18-1LV | rmpA_11(ICEKp1),rmpA_2(KpVP-1) | rmpA2_3-47% | wzi1 | KL1 | Perfect | O1v2 | Very high | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | SHV-11^ | -
GCF_900501255.1 | Klebsiella pneumoniae | ST86 | 3 | 0 | - | 0 | - | 0 | iuc 1 | 1 | iro 1 | 1 | rmpA_2(KpVP-1) | rmpA2_4*-50% | wzi2 | KL2 | Very high | O1v1 | Very high | - | - | - | - | - | - | - | - | - | - | - | - | - | - | SHV-187* | - | - | - | - | -
GCF_000019565.1 | Klebsiella variicola subsp. variicola | ST146 | 0 | 0 | - | 0 | - | 0 | - | 0 | - | 0 | - | - | wzi159 | KL30 | Very high | O3/O3a | Very high | - | - | - | - | - | - | - | - | - | - | - | - | - | - | LEN-24*;LEN-24* | - | - | - | - | -
GCF_000492415.1 | Klebsiella quasipneumoniae subsp. quasipneumoniae | ST1437 | 0 | 0 | - | 0 | - | 0 | - | 0 | - | 0 | - | - | wzi185 | KL46 | Low | O3/O3a | Very high | Aac6-Ib;StrA*;StrB* | - | - | - | - | - | - | CatA2* | - | SulII | - | - | DfrA14 | - | - | - | - | - | OKP-A-3* | -
GCF_000492795.1 | Klebsiella quasipneumoniae subsp. similipneumoniae | ST1435 | 0 | 0 | - | 0 | - | 0 | - | 0 | - | 0 | - | - | wzi183 | KL21 | Very high | O12 | Very high | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | OKP-B-7* | -

### Full results (file)

Here are the full Kleborate results (including assembly quality metrics, and allele calls for all genes in the five MLST schemes) for example 3, written to `results_res_kaptive.txt`:

strain | species | species_match | contig_count | N50 | largest_contig | ambiguous_bases | ST | virulence_score | resistance_score | num_resistance_classes | num_resistance_genes | Yersiniabactin | YbST | Colibactin | CbST | Aerobactin | AbST | Salmochelin | SmST | rmpA | rmpA2 | wzi | K_locus | K_locus_problems | K_locus_confidence | K_locus_identity | K_locus_missing_genes | O_locus | O_locus_problems | O_locus_confidence | O_locus_identity | O_locus_missing_genes | Chr_ST | gapA | infB | mdh | pgi | phoE | rpoB | tonB | ybtS | ybtX | ybtQ | ybtP | ybtA | irp2 | irp1 | ybtU | ybtT | ybtE | fyuA | clbA | clbB | clbC | clbD | clbE | clbF | clbG | clbH | clbI | clbL | clbM | clbN | clbO | clbP | clbQ | iucA | iucB | iucC | iucD | iutA | iroB | iroC | iroD | iroN | AGly | Col | Fcyn | Flq | Gly | MLS | Ntmdz | Phe | Rif | Sul | Tet | Tgc | Tmt | Omp | Bla | Bla_Carb | Bla_ESBL | Bla_ESBL_inhR | Bla_broad | Bla_broad_inhR
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---
GCF_002248955.1 | Klebsiella pneumoniae | strong | 73 | 194261 | 362142 | no | ST15 | 0 | 0 | 5 | 4 | - | 0 | - | 0 | - | 0 | - | 0 | - | - | wzi29 | KL107 | ?-+ | None | 87.89% | KL107_05_wzb,KL107_06_wzc,KL107_07_wbaP,KL107_08,KL107_09,KL107_10,KL107_12,KL107_13,KL107_14,KL107_15 | O1/O2v2 | none | Very high | 98.52% |  | ST15 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | Aac3-IId^ | Mcr3-1* | - | GyrA-83F;GyrA-87A;ParC-80I | - | - | - | CatA1^ | - | - | TetA | - | - | - | SHV-28^ | - | - | - | - | -
GCF_003095495.1 | Klebsiella pneumoniae | strong | 676 | 16918 | 71716 | no | ST258 | 0 | 3 | 11 | 17 | - | 0 | - | 0 | - | 0 | - | 0 | - | - | wzi154 | KL107 | ? | Good | 99.94% |  | O2v2 | ? | Good | 98.38% |  | ST258 | 3 | 3 | 1 | 1 | 1 | 1 | 79 | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | Aac3-IId^;AadA2^;Aph3-Ia^;RmtB;Sat-2A;StrA^;StrB | MgrB-62%;PmrB-36% | - | GyrA-83I;ParC-80I | - | Erm42*;MphA | - | CatA1^ | - | SulI;SulII | TetG | - | DfrA12? | OmpK35-25%;OmpK36GD | TEM-1D^ | KPC-2 | CTX-M-14 | - | SHV-11 | -
GCF_000009885.1 | Klebsiella pneumoniae | strong | 2 | 5248520 | 5248520 | no | ST23 | 4 | 0 | 1 | 0 | ybt 2; ICEKp1 | 326 | - | 0 | iuc 1 | 1 | iro 3 | 18-1LV | rmpA_11(ICEKp1),rmpA_2(KpVP-1) | rmpA2_3-47% | wzi1 | KL1 | none | Perfect | 100.00% |  | O1v2 | none | Very high | 99.13% |  | ST23 | 2 | 1 | 1 | 1 | 9 | 4 | 12 | 9 | 7 | 9 | 6 | 5 | 1 | 1 | 6 | 7 | 7 | 6 | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | 1 | 1 | 1 | 1 | 1 | 21 | 2 | 19 | 5 | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | SHV-11^ | -
GCF_900501255.1 | Klebsiella pneumoniae | strong | 134 | 303226 | 623663 | no | ST86 | 3 | 0 | 0 | 0 | - | 0 | - | 0 | iuc 1 | 1 | iro 1 | 1 | rmpA_2(KpVP-1) | rmpA2_4*-50% | wzi2 | KL2 | none | Very high | 99.94% |  | O1v1 | none | Very high | 98.46% |  | ST86 | 9 | 4 | 2 | 1 | 1 | 1 | 27 | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | - | - | - | - | - | - | - | - | - | - | - | - | - | - | SHV-187* | - | - | - | - | -
GCF_000019565.1 | Klebsiella variicola subsp. variicola | strong | 3 | 5641239 | 5641239 | no | ST146 | 0 | 0 | 0 | 0 | - | 0 | - | 0 | - | 0 | - | 0 | - | - | wzi159 | KL30 | * | Very high | 95.40% |  | O3/O3a | none | Very high | 98.72% |  | ST146 | 16 | 24 | 30 | 27 | 36 | 22 | 55 | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | LEN-24*;LEN-24* | - | - | - | - | -
GCF_000492415.1 | Klebsiella quasipneumoniae subsp. quasipneumoniae | strong | 10 | 5263297 | 5263297 | yes | ST1437 | 0 | 0 | 5 | 6 | - | 0 | - | 0 | - | 0 | - | 0 | - | - | wzi185 | KL46 | ?+* | Low | 98.02% |  | O3/O3a | none | Very high | 95.96% |  | ST1437 | 17 | 19 | 69 | 39 | 185 | 21 | 238 | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | Aac6-Ib;StrA*;StrB* | - | - | - | - | - | - | CatA2* | - | SulII | - | - | DfrA14 | - | - | - | - | - | OKP-A-3* | -
GCF_000492795.1 | Klebsiella quasipneumoniae subsp. similipneumoniae | strong | 2 | 5142035 | 5142035 | yes | ST1435 | 0 | 0 | 1 | 0 | - | 0 | - | 0 | - | 0 | - | 0 | - | - | wzi183 | KL21 | * | Very high | 96.81% |  | O12 | none | Very high | 99.01% |  | ST1435 | 18 | 88 | 128 | 116 | 11 | 99 | 237 | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | OKP-B-7* | -


#### Code testing

Unit tests are available in the /test directory of this repository


#### Visualising outputs in microreact

A helper script to assist users in formating their Kleborate results for viewing in Microreact is provided in the /scripts directory of this repository


## Typing from Illumina reads

If you don't have good quality assemblies, MLST assignment for the chromosomal & virulence locus schemes can also be achieved direct from _K. pneumoniae_ reads using [SRST2](https://github.com/katholt/srst2):

* Download the YbST, CbST, AbST, SmST allele sequences and profile tables from the [data directory](https://github.com/katholt/Kleborate/tree/master/kleborate/data) in this repository.
* Install [SRST2](https://github.com/katholt/srst2) if you don't already have it (`git clone https://github.com/katholt/srst2`).
* Run SRST2, setting the `--mlst_scheme` and `--mlst_definitions` to point to the YbST or CbST allele sequences and profile tables.

Note that currently you can only run SRST2 with one MLST scheme at a time, so in order to type MLST, YbST and CbST you will need to run three separate commands:
```
srst2 --input_pe reads_1.fastq.gz reads_2.fastq.gz --output YbST --log --mlst_db ybt_alleles.fasta --mlst_definitions YbST_profiles.txt
srst2 --input_pe reads_1.fastq.gz reads_2.fastq.gz --output CbST --log --mlst_db clb_alleles.fasta --mlst_definitions CbST_profiles.txt
srst2 --input_pe reads_1.fastq.gz reads_2.fastq.gz --output Klebs --log --mlst_db Klebsiella_pneumoniae.fasta --mlst_definitions kpnuemoniae.txt
```



## Contact us

Kleborate is under active development with many other Klebs genomic analysis tools and projects in progress. 

Please get in touch via the GitHub [issues tracker](https://github.com/katholt/Kleborate/issues) if you have any issues, questions or ideas.

For more on our lab, including other software, see [http://holtlab.net](http://holtlab.net)



## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)



<br>

-------------

Stop! Kleborate and listen<br>
ICE<i>Kp</i> is back with my brand-new invention<br>
If there was a problem, Klebs'll solve it<br>
Check out the hook while Klebs evolves it
