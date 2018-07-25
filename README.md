<p align="center"><img src="logo.png" alt="Kleborate" width="400"></p>

Kleborate is a tool to screen _Klebsiella_ genome assemblies for:
 * MLST sequence type
 * species (e.g. _K. pneumoniae_, _K. quasipneumoniae_, _K. variicola_, etc.)
 * ICEKp associated virulence loci: yersiniabactin (_ybt_), colibactin (_clb_)
 * virulence plasmid associated loci: salmochelin (_iro_), aerobactin (_iuc_), hypermucoidy (_rmpA_, _rmpA2_)
 * antimicrobial resistance genes, including quinolone resistance SNPs and colistin resistance truncations
 * K and O capsule types, via [Kaptive](https://github.com/katholt/Kaptive)

A manuscript describing the Kleborate software in full is currently in preparation. 

In the meantime, if you use Kleborate, please cite the component schemes that you report:<br>
> Yersiniabactin and colibactin (ICEKp) [__Lam MMC, et al.__ Genetic diversity, mobilisation and spread of the yersiniabactin-encoding mobile element ICEKp in _Klebsiella_ pneumoniae populations _MGen_  2018.](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000196)

> Kaptive for capsule (K) serotyping:
[Wyres, K. et al. Identification of _Klebsiella_ capsule synthesis loci from whole genome data. Microbial Genomics (2016).](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000102)

> Kaptive for O antigen serotyping:
[Kaptive Web: user-friendly capsule and lipopolysaccharide serotype prediction for _Klebsiella_ genomes. Journal of Clinical Microbiology (2018).](http://jcm.asm.org/content/56/6/e00197-18)


## Table of Contents

* [Background](#background)
* [Requirements](#requirements)
* [Installation](#installation)
* [Basic usage](#basic-usage)
* [Full usage](#full-usage)
* [Screening details](#screening-details)
   * [MLST](#mlst)
   * [Virulence loci](#virulence-loci)
   * [Other virulence genes](#other-virulence-genes)
   * [Resistance gene detection](#resistance-gene-detection)
   * [Scores and counts](#scores-and-counts)
   * [<em>Klebsiella</em> species](#klebsiella-species)
   * [Capsule typing](#capsule-typing)
* [Example output](#example-output)
   * [Test data](#test-data)
   * [Concise results (stdout)](#concise-results-stdout)
   * [Full results (file)](#full-results-file)
* [Typing from Illumina reads](#typing-from-illumina-reads)
* [Contact us](#contact-us)
* [License](#license)



## Background

_Klebsiella pneumoniae_ (_Kp_) is a commensal bacterium that causes opportunistic infections, with a handful of hypervirulent lineages recognised as true human pathogens. Evidence is now mounting that other _Kp_ strains carrying acquired siderophores (yersiniabactin, salmochelin and aerobactin) and/or the genotoxin colibactin are also highly pathogenic and can cause invasive disease.

Our goal is to help identify emerging pathogenic _Kp_ lineages, and to make it easy for people who are using genomic surveillance to monitor for antibiotic resistance to also look out for the convergence of antibiotic resistance and virulence. To help facilitate that, in this repo we share code for genotyping virulence and resistance genes in _K. pneumoniae_. A table of pre-computed results for 2500 public Klebs genomes is also provided in the [data directory](https://github.com/katholt/Kleborate/tree/master/data).



## Requirements

Software requirements:
* Python (either 2.7 or 3)
* [setuptools](https://pypi.python.org/pypi/setuptools) (required to install Kleborate)
  * To install: `pip install setuptools`
* BLAST+ command line tools (`makeblastdb`, `blastn`, etc.)
  * Version 2.2.30 or later is needed, as earlier versions have a bug with the `culling_limit` parameter.
* [Mash](https://github.com/marbl/Mash) is required to use the `--species` option

As input, Kleborate takes _Klebsiella_ genome assemblies (either completed or draft). If you have unassembled reads, try assembling them with our [Unicycler](https://github.com/rrwick/Unicycler) assembler which works great on Illumina or hybrid Illumina + Nanopore/PacBio reads).



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



## Basic usage

__Screen some genomes for MLST and virulence loci:__<br>
`kleborate -o results.txt -a *.fasta`

__Also screen for resistance genes:__<br>
`kleborate --resistance -o results.txt -a *.fasta`

__Turn on all of Kleborate's optional screens (resistance genes, species check and both K and O loci):__<br>
`kleborate --all -o results.txt -a *.fasta`



## Full usage

```
usage: kleborate -a ASSEMBLIES [ASSEMBLIES ...] [-r] [-s] [--kaptive_k]
                 [--kaptive_o] [-k] [--all] [-o OUTFILE]
                 [--kaptive_k_outfile KAPTIVE_K_OUTFILE]
                 [--kaptive_o_outfile KAPTIVE_O_OUTFILE] [-h] [--version]

Kleborate: a tool for characterising virulence and resistance in Klebsiella

Required arguments:
  -a ASSEMBLIES [ASSEMBLIES ...], --assemblies ASSEMBLIES [ASSEMBLIES ...]
                        FASTA file(s) for assemblies

Screening options:
  -r, --resistance      Turn on resistance genes screening (default: no
                        resistance gene screening)
  -s, --species         Turn on Klebsiella species identification (requires
                        Mash, default: no species identification)
  --kaptive_k           Turn on Kaptive screening of K loci (default: do not
                        run Kaptive for K loci)
  --kaptive_o           Turn on Kaptive screening of O loci (default: do not
                        run Kaptive for O loci)
  -k, --kaptive         Equivalent to --kaptive_k --kaptive_o
  --all                 Equivalent to --resistance --species --kaptive

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

### MLST

Multilocus sequencing typing of _Klebsiella_ follows the schemes described at the [_Klebsiella pneumoniae_ BIGSdb hosted at the Pasteur Institute](http://bigsdb.pasteur.fr/klebsiella/klebsiella.html). The alleles and schemes are stored in the [data directory](https://github.com/katholt/Kleborate/tree/master/kleborate/data) of this repository.

Some notes on Kleborate's MLST calls:
* Kleborate makes an effort to report the closest matching ST / clonal group if a precise match is not found.
* Imprecise allele matches are indicated with a `*`.
* Imprecise ST calls are indicated with `-nLV`, where n indicates the number of loci that disagree with the ST reported. So `258-1LV` indicates a single-locus variant of (SLV) of ST258, i.e. 6/7 loci match ST258.



### Virulence loci

We recently explored the diversity of the _Kp_ integrative conjugative element (ICE<i>Kp</i>), which mobilises the yersiniabactin locus _ybt_, using genomic analysis of a diverse set of 2498 _Kp_ (see [this preprint](http://biorxiv.org/content/early/2017/01/04/098178)). Overall, we found _ybt_ in about a third of all _Kp_ genomes and _clb_ in about 14%. We identified 17 distinct lineages of _ybt_ embedded within 14 structural variants of ICE<i>Kp</i> (some of which include the colibactin _clb_ or salmochelin _iro_ synthesis loci, annotated reference sequences for each ICE<i>Kp</i> variant are included in the [data directory](https://github.com/katholt/Kleborate/tree/master/kleborate/data) of this repository) that can integrate at any of four tRNA-Asn sites in the chromosome. Three of these 17 lineages were associated with three distinct lineages of colibactin. Our analyses reveal hundreds of ICE<i>Kp</i> transmission events affecting hundreds of chromosomal _Kp_ lineages, including nearly two dozen transfers into the globally disseminated carbapenem-resistant clonal group 258. Additionally, we identify a lineage of _ybt_ that is plasmid-encoded, representing a new mechanism for _ybt_ dispersal in _Kp_ populations. Based on this analysis, we developed a MLST-style approach for identifying virulence loci from genome data. 

We further explored the genetic diversity of the aerobactin (_iuc_) and salmochelin (_iro_) loci, both of which are often mobilised together by plasmids, and occasionally found on the chromosome (_iro_ with _ybt_ mobilised by ICE<i>Kp1</i> and a chromosomal _iuc_ locus associated with ST67 _Kp_ subsp _rhinoscleromatis_). We developed a MLST-style typing scheme for these loci, and identified five _iro_ and six _iuc_ lineages. The most common represented plasmid lineages were those related to the virulence plasmid pK2044 from NTUH-K2044 (_iuc1_ and _iro1_) and Kp52.145 plasmid II from Kp52.145 (_iuc2_ and _iro2_). Some lineages (_iuc2A_, _iuc3_ and _iro4_) were associated with novel plasmids that have not yet been previously described in _Kp_, and while others appears to be closely related to to _E. coli_ antimicrobial resistance plasmids with _iuc5_ and _iro5_.  ICE<i>Kp1</i> was associated with _iro3_ and _Kp_ subsp _rhinoscleromatis_ with _iuc4_.

Please note that the aerobactin iuc and salmochelin iro lineage names have been updated between Kleborate version 0.2.0 and 0.3.0 in accordance to the lineage names and associated mobile genetic elements (as listed in the preprint), with the update as follows:

v0.2.0       v0.3.0
--------------------
iuc 1        iuc 5
iuc 2        iuc 1
iuc 3A       iuc 2A
iuc 3B       iuc 2
iuc 4        iuc 3
iuc 5        iuc 4

iro 1        iro 5
iro 2        iro 4
iro 3        iro 1
iro 4        iro 2
iro 5        iro 3

Kleborate examines four different virulence loci in _Klebsiella_: yersiniabactin (_ybt_), colibactin (_clb_), aerobactin (_iuc_) and salmochelin (_iro_).
* For each virulence locus, Kleborate will call the sequence type using the same logic as the MLST described above.
* If the locus is not detected, Kleborate reports the ST as `0` and the lineage as `-`.
* Kleborate will also report the lineage associated with the virulence sequence types (and for yersiniabactin, ICE<i>Kp</i> structure as well).



### Other virulence genes

In addition to the virulence loci, Kleborate screens for the _rmpA_ and _rmpA2_ genes which result in a hypermucoid phenotype by upregulating capsule production.

* Results in these columns should be interpreted as simply binary yes/no calls.
* By default, a gene is called as present if it is detected in a single sequence with >90% identity and >80% coverage of the allele sequence from the virulence plasmid pLVPK. Note that _rmpA_ and _rmpA2_ are ~85% homologous and are reported separately.
* If multiple hits to the same query sequence are found in a given assembly, we attempt to report these separately. The NTUH-K2044 genome carries _rmpA_ in two locations (one on the virulence plasmid, one in the ICE<i>Kp1</i>), and this should be reported as `rmpA;rmpA` as seen in the example above.



### Resistance gene detection

By using the `--resistance` option, Kleborate will screen for resistance genes against the ARG-Annot database of acquired resistance genes ([SRST2](https://github.com/katholt/srst2) version), which includes allelic variants. It attempts to report the best matching variant for each locus in the genome:
* Imprecise allele matches are indicated with `*`.
* If the length of match is less than the length of the reported allele (i.e. a partial match), this is indicated with `?`.
* Note that narrow spectrum beta-lactamases AmpH and SHV () are core genes in _K. pneumoniae_ and so should be detected in most genomes.
  * These genes include: SHV (_K. pneumoniae_), LEN (_K. variicola_), OKP (_K. quasipneumoniae_) and AmpH (all of the above species) 
  * See [this paper](http://www.pnas.org/content/112/27/E3574.long) for more information.
* Note that _oqxAB_ are also core genes  in _K. pneumoniae_, but have been removed from this version of the ARG-Annot DB as they don't actually confer resistance to fluoroquinolones

Using the `--resistance` option also turns on screening for resistance-conferring mutations:
* Fluoroquinolone resistance SNPs: GyrA 83 & 87 and ParC 80 & 84.
* Colistin resistance due to truncation or loss of MgrB or (less than 90% gene coverage counts as a truncation/loss).

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
* Phe (phenicols)
* Rif (rifampin)
* Sul (sulfonamides)
* Tet (tetracyclines)
* Tmt (trimethoprim)


### Scores and counts

Kleborate outputs a simple categorical virulence score, and if resistance screening is enabled, an antimicrobial resistance score as well. These scores provide a rough categorisation of the strains to facilitate monitoring resistance-virulence convergence:

* The virulence score ranges from 0 to 5:
  * 0 = no virulence loci
  * 1 = yersiniabactin only
  * 2 = yersiniabactin and colibactin, or colibactin only 
  * 3 = aerobactin and/or salmochelin only (without yersiniabactin or colibactin)
  * 4 = aerobactin and/or salmochelin with yersiniabactin (without colibactin)
  * 5 = yersiniabactin, colibactin and aerobactin and/or salmochelin
* The resistance score ranges from 0 to 3:
  * 0 = no ESBL, no carbapenemase (regardless of colistin resistance)
  * 1 = ESBL, no carbapenemase (regardless of colistin resistance)
  * 2 = Carbapenemase without colistin resistance (regardless of ESBL)
  * 3 = Carbapenemase with colistin resistance (regardless of ESBL)

When resistance screening is enabled, Kleborate also quantifies how many resistance genes are present and how many resistance classes have at least one gene. Since a resistance class can have multiple genes (as is often the case for the intrinsic genes in the Bla class), the gene count is typically higher than the class count.



### _Klebsiella_ species

By using the `--species` option, Kleborate will attempt to identify the species of _Klebsiella_. It does this by comparing the assembly using Mash to a curated set of _Klebsiella_ assemblies [from NCBI](https://www.ncbi.nlm.nih.gov/assembly) and reporting the species of the closest match. Kleborate considers a Mash distance of ≤ 0.01 to be a strong species match. A distance of > 0.01 and ≤ 0.03 is a weak match and might indicate that your sample is a novel lineage or a hybrid between multiple _Klebsiella_ species. 

Here is an annotated tree of the reference assemblies, made by [mashtree](https://github.com/lskatz/mashtree):
<p align="center"><img src="images/species_tree.png" alt="Klebsiella species tree" width="90%"></p>

Kleborate is designed for the well-studied group of species at the top right of the tree which includes the 'big three': _pneumoniae_, _quasipneumoniae_ (two subspecies) and _variicola_. _K. quasivariicola_ is more recently characterised and described here: [Long 2017](http://genomea.asm.org/content/5/42/e01057-17). The Kp5 group does not yet have a species name and was described in this paper: [Blin 2017](http://onlinelibrary.wiley.com/doi/10.1111/1462-2920.13689/abstract). More distant _Klebsiella_ species (_oxytoca_, _michiganensis_, _grimontii_ and _aerogenes_) are also included, but the virulence profiles of these are less well characterised and deserve further attention.

Kleborate will also call other species in Enterobacteriaceae, as different species sometimes end up in _Klebsiella_ collections. These names are again assigned based on the clades in a mashtree, but were not as carefully curated as the _Klebsiella_ species (so take them with a grain of salt).



### Capsule typing

By default, Kleborate will report the closest match amongst the _wzi_ alleles in the BIGSdb. This is a marker of capsule locus (KL) type, which is highly predictive of capsule (K) serotype. Although there is not a 1-1 relationship between wzi allele and KL/K type, there is a strong correlation (see [Wyres et al, MGen 2016](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000102)). The _wzi_ allele can provide a handy way of spotting the virulence-associated types (wzi=K1, wzi2=K2, wzi5=K5); or spotting capsule switching within clones, e.g. you can tell which ST258 lineage you have from the wzi type (wzi154: the main lineage II; wzi29: recombinant lineage I; others: probably other recombinant lineages).

You can optionally turn on capsule typing using the dedicated capsule typing tool [Kaptive](https://github.com/katholt/Kaptive):
* `--kaptive_k` turns on Kaptive screening of the K locus
* `--kaptive_o` turns on Kaptive screening of the O locus
* `--kaptive` turns on both (is equivalent to `--kaptive_k --kaptive_o`)

This will significantly increase the runtime of Kleborate, but provide much more detailed information about the K and/or O loci and their genes.


## Example output

### Test data

Run these commands to download some well-known _Klebsiella_ genomes and run Kleborate with all optional screens enabled:

```bash
wget -O NTUH-K2044.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/009/885/GCA_000009885.1_ASM988v1/GCA_000009885.1_ASM988v1_genomic.fna.gz
wget -O Klebs_HS11286.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/240/185/GCA_000240185.2_ASM24018v2/GCA_000240185.2_ASM24018v2_genomic.fna.gz
wget -O MGH78578.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/016/305/GCA_000016305.1_ASM1630v1/GCA_000016305.1_ASM1630v1_genomic.fna.gz

kleborate  --all -o results.txt -a *.fasta.gz
```


### Concise results (stdout)

These are the concise Kleborate results that it prints to the terminal:

| strain        | species               | ST   | virulence_score | resistance_score | Yersiniabactin | YbST | Colibactin | CbST | Aerobactin | AbST | Salmochelin | SmST   | hypermucoidy | wzi   | K_locus | K_locus_confidence | O_locus | O_locus_confidence | AGly                                       | Col | Fcyn | Flq               | Gly | MLS | Ntmdz | Phe          | Rif | Sul        | Tet  | Tmt     | Bla                   | Bla_Carb | Bla_ESBL          | Bla_ESBL_inhR | Bla_broad | Bla_broad_inhR          | 
|---------------|-----------------------|------|-----------------|------------------|----------------|------|------------|------|------------|------|-------------|--------|--------------|-------|---------|--------------------|---------|--------------------|--------------------------------------------|-----|------|-------------------|-----|-----|-------|--------------|-----|------------|------|---------|-----------------------|----------|-------------------|---------------|-----------|-------------------------| 
| Klebs_HS11286 | Klebsiella pneumoniae | ST11 | 1               | 2                | ybt 9; ICEKp3  | 15   | -          | 0    | -          | 0    | -           | 0      | -            | wzi74 | KL103   | Very high          | O2v1    | Very high          | StrB;StrA*;AadA2*;RmtB;Aac3-IId*?          | -   | -    | ParC-80I;GyrA-83I | -   | -   | -     | -            | -   | SulII      | TetG | DfrA12? | AmpH*                 | KPC-2    | CTX-M-14;CTX-M-14 | -             | SHV-11    | TEM-30*;TEM-30*;TEM-30* | 
| MGH78578      | Klebsiella pneumoniae | ST38 | 0               | 1                | -              | 0    | -          | 0    | -          | 0    | -           | 0      | -            | wzi50 | KL52    | Perfect            | OL101   | High               | AadA1-pm*?;Aac6-Ib;StrB;Aph3''Ia;StrA;AadB | -   | -    | GyrA-83Y          | -   | -   | -     | CmlA5;CatA1* | -   | SulI;SulII | TetD | -       | AmpH*;SHV-187*;OXA-9* | -        | SHV-12            | -             | -         | TEM-54*;TEM-30*         | 
| NTUH-K2044    | Klebsiella pneumoniae | ST23 | 3               | 0                | ybt 2; ICEKp1  | 326  | -          | 0    | iuc 2      | 1    | iro 5       | 18-1LV | rmpA;rmpA    | wzi1  | KL1     | Perfect            | O1v2    | Very high          | -                                          | -   | -    | -                 | -   | -   | -     | -            | -   | -          | -    | -       | AmpH;SHV-190*         | -        | -                 | -             | -         | -                       | 



### Full results (file)

Here are the full Kleborate results, written to `results.txt`:

| strain        | species               | species_match | contig_count | N50     | largest_contig | ST   | virulence_score | resistance_score | num_resistance_classes | num_resistance_genes | Yersiniabactin | YbST | Colibactin | CbST | Aerobactin | AbST | Salmochelin | SmST   | hypermucoidy | wzi   | K_locus | K_locus_problems | K_locus_confidence | K_locus_identity | K_locus_missing_genes | O_locus | O_locus_problems | O_locus_confidence | O_locus_identity | O_locus_missing_genes | Chr_ST | gapA | infB | mdh | pgi | phoE | rpoB | tonB | ybtS | ybtX | ybtQ | ybtP | ybtA | irp2 | irp1 | ybtU | ybtT | ybtE | fyuA | clbA | clbB | clbC | clbD | clbE | clbF | clbG | clbH | clbI | clbL | clbM | clbN | clbO | clbP | clbQ | AGly                                       | Col | Fcyn | Flq               | Gly | MLS | Ntmdz | Phe          | Rif | Sul        | Tet  | Tmt     | Bla                   | Bla_Carb | Bla_ESBL          | Bla_ESBL_inhR | Bla_broad | Bla_broad_inhR          | 
|---------------|-----------------------|---------------|--------------|---------|----------------|------|-----------------|------------------|------------------------|----------------------|----------------|------|------------|------|------------|------|-------------|--------|--------------|-------|---------|------------------|--------------------|------------------|-----------------------|---------|------------------|--------------------|------------------|-----------------------|--------|------|------|-----|-----|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|--------------------------------------------|-----|------|-------------------|-----|-----|-------|--------------|-----|------------|------|---------|-----------------------|----------|-------------------|---------------|-----------|-------------------------| 
| Klebs_HS11286 | Klebsiella pneumoniae | strong        | 7            | 5333942 | 5333942        | ST11 | 1               | 2                | 9                      | 17                   | ybt 9; ICEKp3  | 15   | -          | 0    | -          | 0    | -           | 0      | -            | wzi74 | KL103   | *                | Very high          | 96.69%           |                       | O2v1    | none             | Very high          | 97.72%           |                       | ST11   | 3    | 3    | 1   | 1   | 1    | 1    | 4    | 14   | 11   | 14   | 5    | 9    | 22   | 19   | 10   | 5    | 11   | 11   | -    | -    | -    | -    | -    | -    | -    | -    | -    | -    | -    | -    | -    | -    | -    | StrB;StrA*;AadA2*;RmtB;Aac3-IId*?          | -   | -    | ParC-80I;GyrA-83I | -   | -   | -     | -            | -   | SulII      | TetG | DfrA12? | AmpH*                 | KPC-2    | CTX-M-14;CTX-M-14 | -             | SHV-11    | TEM-30*;TEM-30*;TEM-30* | 
| MGH78578      | Klebsiella pneumoniae | strong        | 6            | 5315120 | 5315120        | ST38 | 0               | 1                | 7                      | 15                   | -              | 0    | -          | 0    | -          | 0    | -           | 0      | -            | wzi50 | KL52    | none             | Perfect            | 100.00%          |                       | OL101   | *                | High               | 94.91%           |                       | ST38   | 2    | 1    | 2   | 1   | 2    | 2    | 2    | -    | -    | -    | -    | -    | -    | -    | -    | -    | -    | -    | -    | -    | -    | -    | -    | -    | -    | -    | -    | -    | -    | -    | -    | -    | -    | AadA1-pm*?;Aac6-Ib;StrB;Aph3''Ia;StrA;AadB | -   | -    | GyrA-83Y          | -   | -   | -     | CmlA5;CatA1* | -   | SulI;SulII | TetD | -       | AmpH*;SHV-187*;OXA-9* | -        | SHV-12            | -             | -         | TEM-54*;TEM-30*         | 
| NTUH-K2044    | Klebsiella pneumoniae | strong        | 2            | 5248520 | 5248520        | ST23 | 3               | 0                | 0                      | 0                    | ybt 2; ICEKp1  | 326  | -          | 0    | iuc 2      | 1    | iro 5       | 18-1LV | rmpA;rmpA    | wzi1  | KL1     | none             | Perfect            | 100.00%          |                       | O1v2    | none             | Very high          | 99.13%           |                       | ST23   | 2    | 1    | 1   | 1   | 9    | 4    | 12   | 9    | 7    | 9    | 6    | 5    | 1    | 1    | 6    | 7    | 7    | 6    | -    | -    | -    | -    | -    | -    | -    | -    | -    | -    | -    | -    | -    | -    | -    | -                                          | -   | -    | -                 | -   | -   | -     | -            | -   | -          | -    | -       | AmpH;SHV-190*         | -        | -                 | -             | -         | -                       | 



## Typing from Illumina reads

MLST assignment can also be achieved direct from reads using [SRST2](https://github.com/katholt/srst2):

* Download the YbST and CbST allele sequences and profile tables from the [data directory](https://github.com/katholt/Kleborate/tree/master/kleborate/data) in this repository.
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
