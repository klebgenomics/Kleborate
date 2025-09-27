<p align="center"><picture><source srcset="images/logo-dark.png" media="(prefers-color-scheme: dark)"><img src="images/logo.png" alt="Kleborate logo" width="400"></picture></p>

Kleborate was primarily developed to screen genome assemblies of _Klebsiella pneumoniae_ and the _Klebsiella pneumoniae_ species complex (KpSC) for:

* Species (e.g. _K. pneumoniae_, _K. quasipneumoniae_, _K. variicola_, etc.)
* _K. pneumoniae_ species complex MLST
* ICEKp associated virulence loci: yersiniabactin (_ybt_), colibactin (_clb_), salmochelin (_iro_), hypermucoidy (_rmp_)
* Virulence plasmid associated loci: salmochelin (_iro_), aerobactin (_iuc_), hypermucoidy (_rmp_, _rmpA2_)
* Antimicrobial resistance determinants: acquired genes, SNPs, gene truncations and intrinsic β-lactamases
* Prediction of phenotypic resistance to ciprofloxacin
* K (capsule) and O antigen (LPS) serotype prediction, via _wzi_ alleles and [Kaptive](https://github.com/klebgenomics/Kaptive)

Kleborate v3 includes a rewrite of the code to (i) replace the use of BLAST with [minimap2](https://lh3.github.io/minimap2/minimap2.html) (faster and less buggy); and (ii) introduce a modular structure making it easy to add new typing modules, including for other species. Currently, functionality for other species is limited to MLST for _Klebsiella oxytoca_ species complex and _Escherichia coli_ but more is in development.

A list of changes from v2 is available in the [documentation](https://kleborate.readthedocs.io/en/latest/index.html#changes-from-v2).

**Documentation**

For information on how to install and run Kleborate v3, please visit the [Docs](https://kleborate.readthedocs.io/en/latest/).

**Citation**

If you use Kleborate, please cite the paper: Lam, MMC. et al. A genomic surveillance framework and genotyping tool for _Klebsiella pneumoniae_ and its related species complex, Nature Communications (2021). [https://doi.org/10.1038/s41467-021-24448-3](https://doi.org/10.1038/s41467-021-24448-3)

If you use the Kaptive calls for K and O locus typing please also cite: Wyres, KL. et al. Identification of _Klebsiella_ capsule synthesis loci from whole genome data. Microbial Genomics (2016). [http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000102](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000102)

If you use the ciprofloxacin resistance predictions, please cite the paper: Tsang et al. Ciprofloxacin resistance in _Klebsiella pneumoniae_: phenotype prediction from genotype and global distribution of resistance determinants. BioRxiv (2025). [https://doi.org/10.1101/2025.09.24.678318](https://doi.org/10.1101/2025.09.24.678318)
