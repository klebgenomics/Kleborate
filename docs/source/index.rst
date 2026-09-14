.. Kleborate documentation master file, created by
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

**For kpsc typing please cite**:

1. Lam, MMC. et al. A genomic surveillance framework and genotyping tool for *Klebsiella pneumoniae* and its related species complex, *Nature Communications* (2021). `<https://www.nature.com/articles/s41467-021-24448-3>`_

2. Stanton, TD. et al. Fast and accurate in silico antigen typing with Kaptive 3, *Microbial Genomics* (2025). `<https://doi.org/10.1099/mgen.0.001428>`_

3. Bogaerts, B. et al. MiST: rapid, accurate and flexible (core-genome) multi-locus sequence typing (MLST) allele calling from draft genomes, *BMC Genomics* (2025). `<https://doi.org/10.1186/s12864-025-12324-z>`_

4. Brisse, S. et al. Virulent clones of *Klebsiella pneumoniae*: Identification and evolutionary scenario based on genomic and phenotypic characterization, *PLoS ONE* (2009). `<https://doi.org/10.1371/journal.pone.0004982>`_

5. Diancourt, L. et al. Multilocus sequence typing of *Klebsiella pneumoniae* nosocomial isolates, *Journal of Clinical Microbiology* (2005). `<https://doi.org/10.1128/JCM.43.8.4178-4182.2005>`_

6. Hennart, M. et al. A Dual Barcoding Approach to Bacterial Strain Nomenclature: Genomic Taxonomy of *Klebsiella pneumoniae* Strains, *Molecular Biology and Evolution* (2022). `<https://doi.org/10.1093/molbev/msac135>`_

**For Kosc typing please cite**:

1. Ashcroft, MM. et al. A capsule polysaccharide synthesis locus database for the *Klebsiella oxytoca* Species Complex, *bioRxiv* (2026). `<https://doi.org/10.64898/2026.07.16.739023>`_

**For Escherichia typing please cite**:

1. Wirth, T. et al. Sex and virulence in *Escherichia coli*: an evolutionary perspective, *Molecular Microbiology* (2006). `<https://doi.org/10.1111/j.1365-2958.2006.05172.x>`_

2. Ingle, D. et al. Evolution of atypical enteropathogenic *E. coli* by repeated acquisition of LEE pathogenicity island variants, *Nature Microbiology* (2016). `<https://doi.org/10.1038/nmicrobiol.2015.10>`_

3. Feldgarden, M. et al. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence, *Scientific Reports* (2021). `<https://doi.org/10.1038/s41598-021-91456-0>`_

4. Yassine, I. et al. ShigaPass: an in silico tool predicting *Shigella* serotypes from whole-genome sequencing assemblies, *Microbial Genomics* (2023). `<https://doi.org/10.1099/mgen.0.000961>`_

5. Bessonov, K. et al. ECTyper: in silico *Escherichia coli* serotype and species prediction from raw and assembled whole-genome sequence data, *Microbial Genomics* (2021). `<https://doi.org/10.1099/mgen.0.000728>`_



The following papers provide more information on the component schemes and genotyping incorporated in Kleborate:

..
   
   Yersiniabactin and colibactin (*ICEKp*):
   Lam, MMC. et al. Genetic diversity, mobilisation and spread of the yersiniabactin-encoding mobile element *ICEKp* in *Klebsiella pneumoniae* populations. *Microbial Genomics* (2018). `Microbial Genomics <http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000196>`_

   Aerobactin and salmochelin:
   Lam, MMC. et al. Tracking key virulence loci encoding aerobactin and salmochelin siderophore synthesis in *Klebsiella pneumoniae*. *Genome Medicine* (2018). `Genome Medicine <https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-018-0587-5>`_

   Kaptive for capsule (K) serotyping:
   Stanton, TD. et al. Fast and accurate in silico antigen typing with Kaptive 3, *Microbial Genomics* (2025). `<https://doi.org/10.1099/mgen.0.001428>`_

   Kaptive for capsule (K) serotyping:
   Wyres, KL. et al. Identification of *Klebsiella* capsule synthesis loci from whole genome data. *Microbial Genomics* (2016). `Microbial Genomics 2 <http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000102>`_

   Kaptive for O antigen (LPS) serotyping:
   Wick, RR et. al. Kaptive Web: user-friendly capsule and lipopolysaccharide serotype prediction for *Klebsiella* genomes. *Journal of Clinical Microbiology* (2018). `Journal of Clinical Microbiology <http://jcm.asm.org/content/56/6/e00197-18>`_


Changes from v2
----------------

When Kleborate v3 is run using the ``-p kpsc`` option to run preset modules for *K. pneumoniae* the same logic is implemented as Kleborate v2, plus the following changes/updates:

* Updated MLST & virulence databases 
* Column ``Chr_ST``  has been removed in v3, as it is redundant with ``ST`` 
* Updated AMR database to Kleborate_AMRdb_v3.3
* Mutations are reported using  `HGVS nomenclature <https://hgvs-nomenclature.org/stable/recommendations/protein/substitution/>`_
* Added ``p.(Met1?)`` to indicate when pmrB or mgrB have a mutation in the start codon that may disrupt translation (in ``Col_mutations`` column)
* Added check for synonymous mutation in ompK36 (25 C > T) associated with increased resistance to carbapenems (in ``Omp_mutations`` column)
* Added checks for ``D`` ompK36 loop 3 (L3) insertions
* Added new AMR genotyping report compatible with the `hAMRonization <https://github.com/pha4ge/hAMRonization/blob/master/schema/PHA4GE%20AMR%20Gene%20%26%20Variant%20Specification.csv>`_ standard developed by the Public Health Alliance for Genomic Epidemiology (PHA4GE).
* Added a module for Ciprofloxacin resistance prediction
* Updated *rmp* typing to detect expression of the *rmp* locus
* Added new module for typing *peg-344* gene
* Added new module for cgMLST and Lin codes
* Added new module for typing *mrk* operon
* Updated to use Kaptive v3, which has some changes to the names of output variables:
   ``K_locus_missing_genes``  has been renamed ``K_Missing_expected_genes`` 
   ``O_locus_missing_genes``  has been renamed ``O_Missing_expected_genes`` 
* Updated assembly statistics module to use qualibact-v1.0 curated thresholds
* Added new genome specification report compatible with the `PHA4GE Microbial Genotyping Data Specification <https://github.com/pha4ge/genotyping-specification>`_
* Added new module for KoSC K and O locus typing 

* Added new modules for *Escherichia* species: pathotyping, typing of the LEE pathogenicity island, ClermonTyping, typing of stx types using StxTyper, O:H serotyping using ECTyper, pks typing, AMR typing using AMRFinderPlus, typing of Group 2 and 3 CPS using Kaptive



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



