.. KleborateModular documentation master file, created by
   sphinx-quickstart on Thu Apr 25 06:02:56 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. toctree::
   :maxdepth: 1
   :hidden:

   Installation
   Usage
   modules
   Creating-New-Modules


########################
Introducing Kleborate v3
########################

Kleborate was primarily developed to screen genome assemblies of *Klebsiella pneumoniae* and the *Klebsiella pneumoniae* species complex (KpSC) for:

* Species (e.g. *K. pneumoniae*\ , *K. quasipneumoniae*\ , *K. variicola*\ , etc.)
* *K. pneumoniae* MLST sequence type
* *ICEKp* associated virulence loci: yersiniabactin (*ybt*), colibactin (*clb*), salmochelin (*iro*), hypermucoidy (*rmpA*)
* Virulence plasmid associated loci: salmochelin (\ *iro*\ ), aerobactin (\ *iuc*\ ), hypermucoidy (\ *rmpA*\ , *rmpA2*\ )
* Antimicrobial resistance determinants: acquired genes, SNPs, gene truncations and intrinsic β-lactamases
* K (capsule) and O antigen (LPS) serotype prediction, via *wzi* alleles and `Kaptive <https://github.com/klebgenomics/Kaptive>`_


`Kleborate v3 <https://github.com/klebgenomics/Kleborate>`_ includes a rewrite of the code to (i) replace the use of BLAST with minimap (faster and less buggy); and (ii) introduce a modular structure making it easy to add new typing modules, including for other species.


For *K. pneumoniae* species complex, **Kleborate v3 can reproduce the outputs of Kleborate v2 by running the preset modules for KpSC via:**
 

.. code-block:: Python

   kleborate -a *.fasta -o kleborate_results -p kpsc --trim_headers

(Note the command has changed from Kleborate v2, the above is equivalent to running ``kleborate --all -o results.txt -a *.fasta``  with Kleborate v2 and includes all resistance and Kaptive-based typing)

**New modules for other species are in development,** for now these include MLST schemes for *Klebsiella oxytoca* species complex and *Escherichia coli* (see the Modules page).


Citations
----------

If you use Kleborate, please cite the paper: Lam, MMC. et al. A genomic surveillance framework and genotyping tool for *Klebsiella pneumoniae* and its related species complex, *Nature Communications* (2021). `<https://www.nature.com/articles/s41467-021-24448-3>`_


If you use the Kaptive calls for K and O locus typing please also cite Wyres, KL. et al. Identification of *Klebsiella* capsule synthesis loci from whole genome data. *Microbial Genomics* (2016). `<http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000102>`_

The following papers provide more information on the component schemes and genotyping incorporated in Kleborate:

..
   
   Yersiniabactin and colibactin (*ICEKp*):
   Lam, MMC. et al. Genetic diversity, mobilisation and spread of the yersiniabactin-encoding mobile element *ICEKp* in *Klebsiella pneumoniae* populations. *Microbial Genomics* (2018). `Microbial Genomics <http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000196>`_

   Aerobactin and salmochelin:
   Lam, MMC. et al. Tracking key virulence loci encoding aerobactin and salmochelin siderophore synthesis in *Klebsiella pneumoniae*. *Genome Medicine* (2018). `Genome Medicine <https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-018-0587-5>`_

   Kaptive for capsule (K) serotyping:
   Wyres, KL. et al. Identification of *Klebsiella* capsule synthesis loci from whole genome data. *Microbial Genomics* (2016). `Microbial Genomics 2 <http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000102>`_

   Kaptive for O antigen (LPS) serotyping:
   Wick, RR et. al. Kaptive Web: user-friendly capsule and lipopolysaccharide serotype prediction for *Klebsiella* genomes. *Journal of Clinical Microbiology* (2018). `Journal of Clinical Microbiology <http://jcm.asm.org/content/56/6/e00197-18>`_


Changes from v2
----------------

When Kleborate v3 is run using the ``-p kpsc`` option to run preset modules for *K. pneumoniae* the same logic is implemented as Kleborate v2, plus the following changes/updates:

* MLST & virulence databases updated (April 2024)
* Column ``Chr_ST``  has been removed in v3, as it is redundant with ``ST`` 
* AMR database updated based on CARD v3.2.9 (June 2024)
* Added ``$`` to indicate when PmrB or MgrB have a mutation in the start codon that may disrupt translation (in ``Col_mutations`` column)
* Added check for synonymous mutation in ompK36 (25 C > T) associated with increased resistance to carbapenems (in ``Omp_mutations`` column)
* Updated to use Kaptive v3, which has some changes to the names of output variables:
  * ``K_locus_missing_genes``  has been renamed ``K_Missing_expected_genes`` 
  * ``O_locus_missing_genes``  has been renamed ``O_Missing_expected_genes`` 
* Provided new AMR genotyping report compatible with the `hAMRonization <https://github.com/pha4ge/hAMRonization/blob/master/schema/PHA4GE%20AMR%20Gene%20%26%20Variant%20Specification.csv>`_. standard developed by the Public Health Alliance for Genomic Epidemiology (PHA4GE).



Tutorial
--------

A step-by-step tutorial for Kleborate v3 is available at `kleborate-workshop <https://docs.google.com/document/d/1R61bQbBngpiDB2Gl_eXigePBVakYZEjy/edit>`_, covering:


* Kleborate's features and their scientific rationale
* How to run Kleborate 
* Examples, illustrating how to run and interpret results


Public reports
----------------

The `Kleborate paper <https://www.nature.com/articles/s41467-021-24448-3>`_ reports results of genotyping ~10,000 public genomes that have been filtered to remove redundant sequences (e.g. outbreak clusters, identified as small genome-wide mash distance with same year, location and genotypes), with Kleborate v2. The results can be explored in `Microreact <https://bit.ly/klebMR>`_ (which shows the mash tree, Kleborate output & curated metadata) or `Kleborate-Viz <https://kleborate.erc.monash.edu/>`_ (R shiny app). Kleborate-Viz also has the EuSCAPE dataset preloaded, or you can view your own Kleborate results.

Kleborate is also included in `Klebsiella Pathogenwatch <https://pathogen.watch/>`_ which shows interactive trees, maps and line lists for *Klebsiella pneumoniae* and allows you to analyse your own data in context of the public collections. See `this paper <https://doi.org/10.1093/cid/ciab784>`_ for an example of how to use it.

Contact us
----------

Kleborate is under active development with many other Klebs genomic analysis tools and projects in progress (see `github.com/klebgenomics <https://github.com/klebgenomics>`_). 

Please get in touch via the GitHub `issues tracker <https://github.com/klebgenomics/Kleborate/issues>`_ if you have any issues, questions or ideas.

For more on our lab, including other software, see `http://holtlab.net <http://holtlab.net>`_

License
-------

`GNU General Public License, version 3 <https://www.gnu.org/licenses/gpl-3.0.html>`_



