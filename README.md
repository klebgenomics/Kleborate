<p align="center"><img src="logo.png" alt="Kleborate" width="400"></p>

Kleborate is a tool to screen genome assemblies of _Klebsiella pneumoniae_ and the _Klebsiella pneumoniae_ species complex (KpSC) for:
 * MLST sequence type
 * species (e.g. _K. pneumoniae_, _K. quasipneumoniae_, _K. variicola_, etc.)
 * ICE<i>Kp</i> associated virulence loci: yersiniabactin (_ybt_), colibactin (_clb_), salmochelin (_iro_), hypermucoidy (_rmpA_)
 * virulence plasmid associated loci: salmochelin (_iro_), aerobactin (_iuc_), hypermucoidy (_rmpA_, _rmpA2_)
 * antimicrobial resistance determinants: acquired genes, SNPs, gene truncations and intrinsic β-lactamases
 * K (capsule) and O antigen (LPS) serotype prediction, via _wzi_ alleles and [Kaptive](https://github.com/katholt/Kaptive)


### Installation, usage and documentaton
All the info you need to install Kleborate, run it, and interpret results is in the [Kleborate wiki](https://github.com/katholt/Kleborate/wiki)!

A step-by-step tutorial, illustrating how to use Kleborate and interpret the data, is available at [bit.ly/kleborate-workshop](bit.ly/kleborate-workshop)

### Citations
If you use Kleborate, please cite the paper: [Lam, MMC. et al. A genomic surveillance framework and genotyping tool for Klebsiella pneumoniae and its related species complex. Nature Communications (2021)](https://www.nature.com/articles/s41467-021-24448-3). 

If you turn on the [Kaptive](https://github.com/katholt/Kaptive) option for full K and O typing, please also cite Kaptive directly:
[Identification of _Klebsiella_ capsule synthesis loci from whole genome data. Microbial Genomics (2016).](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000102)

### Related tools
You may also be interested in [Kleborate-viz](https://kleborate.erc.monash.edu/), a ShinyR app for visualising Kleborate output.

Users not comfortable runnnig Kleborate locally via the commandline may like to try uploading your genomes (reads or assemblies) to the drag-and-drop genomic epi online platform [PathogenWatch](https://pathogen.watch/) - all genomes identified as _Klebsiella pneumoniae_ will automatically be passed to Kleborate and Kaptive to generate a genotyping report, which can be downloaded in tabular format compatible with Kleborate-viz. See this preprint for more about [Klebsiella Pathogenwatch](https://www.biorxiv.org/content/10.1101/2021.06.22.448967v2).

