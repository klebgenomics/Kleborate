
****************************************************
Modules for *Escherichia* species
****************************************************


.. code-block:: Python

   --preset escherichia

These modules will be run if the ``enterobacterales__species``\   module confirms the input assembly as a member of the *Escherichia* genus. 

.. _escherichia__mlst_achtman:

.. _escherichia__mlst_pasteur:

*E. coli* MLST
------------

.. code-block:: Python

   -m escherichia__mlst_achtman, escherichia__mlst_pasteur

Genomes identified as belonging to the *Escherichia* genus are subjected to MLST using two 7-locus schemes.

The Pasteur scheme is described at the *Escherichia coli* `Database at Pasteur Institute <https://bigsdb.pasteur.fr/ecoli/>`_. Please see `here <https://bigsdb.pasteur.fr/ecoli/references/>`_ for more details.

The Achtman scheme is hosted at `\Enterobase <https://enterobase.warwick.ac.uk/>`_.

The genes included in each scheme are noted in the Outputs table below.

A copy of the MLST alleles and ST definitions used in each module is stored in the ``/data``  directory of the module.


*E. coli* MLST parameters
+++++++++++++++++++++++++++

``--escherichia_mlst_achtman_min_identity`` 

Minimum alignment percent identity for *Escherchia-Achtman* MLST (default: 90.0)

``--escherichia_mlst_achtman_min_coverage`` 

Minimum alignment percent coverage for Escherchia-Achtman MLST (default: 80.0)

``--escherichia_mlst_achtman_required_exact_matches`` 

At least this many exact matches are required to call an ST (default: 3)

``--escherichia_mlst_pasteur_min_identity`` 

Minimum alignment percent identity for Escherchia-Pasteur MLST (default: 90.0)

``--escherichia_mlst_pasteur_min_coverage`` 

Minimum alignment percent coverage for Escherchia-Pasteur MLST (default: 80.0)

``--escherichia_mlst_pasteur_required_exact_matches`` 

At least this many exact matches are required to call an ST (default: 4)

*E. coli* MLST outputs
++++++++++++++++++++++

Output of the Pasteur *E. coli* MLST module is the following columns:

.. list-table::

   * - ST
     - sequence type

   * - dinB, icdA, pabB, polB, putP, trpA, trpB, uidA
     - allele number

Output of the Achtman *E. coli* MLST module is the following columns:

.. list-table::

   * - ST
     - sequence type

   * - adk, fumC, gyrB, icd, mdh, purA, recA
     - allele number

* Kleborate makes an effort to report the closest matching ST if a precise match is not found.
* Imprecise allele matches are indicated with a ``*``.
* Imprecise ST calls are indicated with ``-nLV``\ , where n indicates the number of loci that disagree with the ST reported. So ``258-1LV`` indicates a single-locus variant (SLV) of ST258, i.e. 6/7 loci match ST258.


.. _escherichia__pathovar:

*E. coli* Pathotyping
---------------------

.. code-block:: Python

   -m escherichia__pathovar

*Escherichia coli* is broadly split into 2 groups: intestinal diarrheagenic *E. coli* (DEC), and extra-intestinal E. coli (ExPEC) `see paper <https://pmc.ncbi.nlm.nih.gov/articles/PMC5156508/>`_. Diarrheagenic Escherichia coli (DEC) encompasses enteropathogenic *E. coli* (EPEC), enterotoxigenic *E. coli* (ETEC), enterohaemorrhagic *E. coli* (EHEC), Shiga toxin-producing *E. coli* (STEC), enteroaggregative *E. coli* (EAEC), enteroinvasive *E. coli* (EIEC), and diffusely adherent *E. coli* (DAEC) `see <https://doi.org/10.3389/fcimb.2016.00141/>`_. Additionally, *Shigella* is also regarded as a DEC pathotype as it closely resembles EIEC in terms of virulence attributes and pathogenecity

The majority of DEC pathotypes are defined by the possession of one or more pathotype-specific virulence markers. However, in EAEC, DAEC and AIEC, the role of these markers in virulence is not proven. 

Virulence markers of diarrheagenic *E. coli* 
++++++++++++++++++++++++++++++++++++++++++++++

.. list-table:: 
   :header-rows: 1

   * - **E. coli Pathotype**
     - **Defining marker**
     - **Virulence determinants**
     - **Location of determinants**
     - **Diagnostic targets for PCR**
     - **Other diagnostic targets**

   * - EPEC
     - LEE pathogenicity island
     - LEE pathogenicity island
     - Pathogenicity island
     - eae
     - bfpA

   * - EIEC/Shigella
     - pINV
     - pINV
     - Plasmid
     - ipaH
     - Other ipa genes

   * - ETEC
     - ST or LT
     - ST or LT\nPlus colonisation factors
     - Plasmid; transposon
     - Elt, est
     - 

   * - EHEC
     - Shiga toxin
     - Shiga toxin 1 and 2
     - Prophages
     - Stx1, stx2
     - Eae, ehxA

   * - EAEC
     - pAA; aggregative adhesion
     - Not known
     - Plasmid
     - aggR, aatA, aaiC
     - `-`

   * - DAEC
     - Afa/ Dr adhesins
     - Not known
     - Not known
     - afa/Dr adhesins
     - `-`

   * - AIEC
     - Adherent-invasive phenotype
     - Not known
     - Not known
     - none
     - `-`



This module classifies *E. coli* pathotypes based on the presence or absence of virulence marker genes using a curated database `VirulenceFinder <http://www.genomicepidemiology.org/>`_ gene database.  Input genomes are aligned to the database using Minimap2, Kleborate then applies `pathotype-calling logic described in EnteroBase <https://enterobase.readthedocs.io/en/latest/pipelines/backend-pipeline-phylotypes.html?highlight=pathovar/>`_ to classify pathotypes

Additionally, Kleborate distinguish *Shigella* species based on the serotype-specific O-antigen biosynthetic gene cluster. The module aligns *E. coli* genomes against a curated reference sequence derived from the *Shigella* serotyping pipeline, `shigatyper <https://github.com/CFSAN-Biostatistics/shigatyper>`_ using Minimap2.

The VirulenceFinder DB, shigatyper reference sequence and marker definition are found in the **/data**  directory of this module.


*E. coli* Pathovar parameters
++++++++++++++++++++++++++++++++++

 
``--escherichia__pathovar_min_identity``

Minimum alignment percent identity for pathotype (default: 90.0)

``--escherichia__pathovar_min_coverage``

Minimum alignment percent coverage for pathotype (default: 80.0)


*E. coli* Pathovar outputs
++++++++++++++++++++++++++++

.. list-table:: 
   :header-rows: 1

   * - Pathotype
     - Predicted pathotype

   * - Stx1, Stx2, ST, LT, eae, ipaH
     - Virulence markers


.. _escherichia__mlst_lee:

Typing the LEE pathogenicity island of *E. coli*
----------------------------------------------

.. code-block:: Python

   -m escherichia__mlst_lee

Locus of enterocyte effacement (LEE) is a ~40 kb chromosomal pathogenicity island composed of 41 core genes organized into five operons  `Elliot et al., 1998 <https://onlinelibrary.wiley.com/doi/10.1046/j.1365-2958.1998.00783.x>`_. It encodes an (i) outer membrane adhesive protein, known as intimin protein that encodes eae gene (ii) type III secretion system (T3SS), and (iii) translocated receptor (Tir) as well as translocons, chaperones, regulators and secreted effector proteins that are linked to virulence.

Kleborate includes a module for subtyping of the LEE pathogenicity island. Details of the LEE subtypes and lineages can be found in this `Nature Microbiology paper <https://www.nature.com/articles/nmicrobiol201510>`_.

The LEE typing database is based on analysis of >250 LEE-containing *E. coli* genomes and includes 7 loci (eae (intimin), tir, espA, espB, espD, espH, espZ). The data is provided as a MLST-style database, in which combinations of alleles are assigned to a LEE subtype, to facilitate a common nomenclature for LEE subtypes. Each sequence in the database represents a cluster of closely related alleles that have been assigned to the same locus type. The LEE scheme includes three distinct lineages: Lineage 1 consists of LEE subtypes 1-2; Lineage 2 consists of LEE subtypes 3-8; Lineage 3 consists of LEE subtypes 9-30.

The sequences and  MLST-style profile definitions are stored in the **/data**  directory of this module.


Parameters
++++++++++

``--escherichia__mlst_LEE_min_identity``

Minimum alignment percent identity for escherichia__mlst_LEE (default: 90.0)

``--escherichia__mlst_LEE_min_coverage``

Minimum alignment percent coverage for escherichia_mlst_LEE (default: 80.0)

``escherichia__mlst_LEE_mlst_required_exact_matches``

At least this many exact matches are required to call an ST (default: 3)


*E. coli*  LEE MLST outputs
++++++++++++++++++++++++++++

Output of the *E. coli* LEE MLST module is the following columns:


.. list-table::

   * - LEE_ST
     - Sequence type

   * - LEE_lineage
     - Lineage

   * - LEE_eae, LEE_tir, LEE_espA, LEE_espB, LEE_espD, LEE_espH, LEE_espZ

     - allele number (LEE locus)


* Kleborate makes an effort to report the closest matching ST if a precise match is not found.
* Imprecise allele matches are indicated with a ``*``.
* Imprecise ST calls are indicated with ``-nLV``\ , where n indicates the number of loci that disagree with the ST reported. So ``ST10-3LV`` indicates a three-locus variant (SLV) of ST10, i.e. 4/7 loci match ST10.


.. _escherichia__stxtyper:


Stxtyper
-----------

.. code-block:: Python

   -m escherichia__stxtyper

Shiga toxins (Stxs) are a key virulence factor of Stx-producing *Escherichia coli* (STEC). They are also found in *Shigella dysenteriae 1*. Stxs belong to the AB-type toxin family and are divided into two antigenically distinct groups, Stx1 and Stx2. Each group contains several variants/subtypes six Stx1 (a, b, c, d, e, f) and seven Stx2 (a, b, c, d, e, f, and g) `Yano et al., 2023 <https://www.nature.com/articles/s41598-023-32111-8>`_ and `Melton-Celsa 2014 <https://pmc.ncbi.nlm.nih.gov/articles/PMC4270005/>`_. These toxins are encoded by bacteriophages (lysogenic Stx phage) and STEC strains can produce either single Stx subtype or a combination Stx subtypes.

This module will run StxTyper to determine the stx type. See the `StxTyper documentation <https://github.com/ncbi/stxtyper>`_ for more details of how it works.


StxTyper Outputs
+++++++++++++++++++++

StxTyper results are output in the following columns:

.. list-table::
   :header-rows: 1

   * - Column Name
     - Description
   * - Stx_type
     - The stx type, if operon is complete, stx will be reported as stx1a, for other values of operon stx_type will be stx1, stx2, or just stx if the algorithm can't resolve at all.
   * - operon
     - Status the operon found. (COMPLETE, PARTIAL, PARTIAL_CONTIG_END-for partial operons that could be split by contig boundaries due to sequencing or assembly artifacts, EXTENDED- coding sequence extends beyond the reference stop codon for one or both of the reference proteins, INTERNAL_STOP-for Stx operons where one of the subunits has a nonsense mutation, FRAMESHIFT-where StxTyper detected an indel in the coding sequence, AMBIGUOUS-StxTyper found an ambiguous base in the query sequence, COMPLETE_NOVEL-a full-length stx operon that is not typeable)
   * - identity
     - Percent identity for both A and B subunits
   * - target_start
     - start of the alignments
   * - target_stop
     - End of the alignments
   * - target_strand
     - strand the target is on
   * - A_reference
     - Closest reference protein for the A subunit
   * - A_identity
     - percent identity to the reference for the A subunit
   * - A_reference_subtype
     - Subtype assigned to the reference sequence for the A subunit
   * - A_coverage
     - Percentage of the reference for the A subunit that is covered by the alignment
   * - B_reference
     - Closest reference protein for the B subunit
   * - B_reference_subtype
     - Subtype assigned to the reference sequence for the B subunit
   * - B_identity
     - Percent identity to the reference for the B subunit
   * - B_coverage
     - Percentage of the reference for the B subunit that is covered by the alignment


.. _escherichia__ectyper:

*E. coli* O:H serotyping
----------------------

.. code-block:: Python

   -m escherichia__ectyper

*E. coli* serotypes are defined by combinations of O (lipopolysaccharide) and H (flagellar) antigens. Currently there are ~183 O-groups and 53 H-types that have been defined serologically `Ørskov and Ørskov 1984 <https://www.sciencedirect.com/science/article/abs/pii/S0580951708704471/>`_.


O-antigen 
++++++++++

The O-antigen is an integral component of the Lipopolysaccharide (LPS) found in the outer membrane of the bacteria. LPS comprises three components: lipid A, a core oligosaccharide, and the O-specific polysaccharide chain (O antigen).  The O-antigen domain exhibits significant variability consisting of 10 to 25 repeating oligosaccharide units, with each unit containing two to seven sugar residues `Liu et al., 2020 <https://pmc.ncbi.nlm.nih.gov/articles/PMC7685785/>`_. The genes responsible for synthesis of O-antigens are usually present as a gene cluster and are located between the two chromosomal housekeeping genes galF and gnd/ugd `Iguchi et al 2014 <https://pmc.ncbi.nlm.nih.gov/articles/PMC4379981/>`_. Major pathways involved in the assembly, synthesis and transport of O-antigen include, the Wzy pathway the Wzx/Wzy-dependent pathway, encoded by the wzx (O-antigen flippase) and wzy (O-antigen polymerase) genes, and the ABC transporter pathway, encoded by wzm and wzt. These genes are ideal biomarkers for predicting O antigen types.  


H antigens 
++++++++++

H antigens (flagellar) are surface proteins composed of repeated molecules of the protein flagellin, which facilitate bacterial motility. These antigens are numbered from H1 to H56 (H13, H22, and H50 are not used) and are distinct from the O and K antigens. Flagellin is encoded by the fliC gene on the chromosomal locus or its homologues (non-fliC flagellin-coding genes such as flkA, fllA, and flmA). Of the 53 well know H antigen types, 44 are conferred by expression of the fliC gene,  the remaining 9 H types are  encoded by non-fliC flagellin genes. Specifically H3, H35, H36, H47,and H53 are encoded by flkA, H44 and H55 by fllA, H54 by flmA, and H17 by flnA.


Kleborate uses ECTyper for in silico serotyping. See `ECTyper paper <https://pmc.ncbi.nlm.nih.gov/articles/PMC8767331/>`_. for more details 

Outputs
+++++++

Outputs of the ECTyper module is the following columns:

.. list-table:: 
   :header-rows: 0

   * - O-type
     - O antigen

   * - H-type
     - H antigen

   * - Serotype
     - Predicted O and H antigen(s)

   * - QC
     - The Quality Control value summarising the overall quality of prediction

   * - Evidence
     - How many alleles in total used to both call O and H antigens

   * - GeneScores
     - ECTyper O and H antigen gene scores in 0 to 1 range

   * - AllelesKeys
     - Best matching ECTyper database allele keys used to call the serotype

   * - GeneIdentities(%)
     - %identity values of the query alleles

   * - GeneCoverages(%)
     - %coverage values of the query alleles

   * - GeneLengths
     - allele lengths of the query alleles

   * - Warnings
     - An additional warning linked to the quality control status or any other error message(s).


.. _ClermonTyping:


ClermonTyping
----------------------

.. code-block:: Python

   -m escherichia__ezclermont


*Escherichia* genera is composed of the following clades *Escherichia albertii*, *E. fergusonii*, five cryptic *Escherichia* clades (I–V) and *E. coli* sensu stricto. Additionally, the *E. coli* species can be divided into seven main phylogroups termed A, B1, B2, C, D, E and F. 

Kleborate assigns genomes to different clades and phylogroups using `EzClermont tool <https://pmc.ncbi.nlm.nih.gov/articles/PMC7656184/>`_.The ClermonTyping is based on the concept of in vitro PCR assays.


Parameters
++++++++++

``--escherichia__ezclermont_min_length``

minimum contig length to consider.default: 500


Outputs
+++++++

.. list-table:: 
   :header-rows: 1

   * - **Clermont_type**
     - **Phylotype**

   * - **Clermont_profile**
     - Presence or absence of the PCR product


.. _Escherichia AMR:


Escherichia AMR
------------------------

.. code-block:: Python

   -m escherichia__amr


This module screens input genomes against Antimicrobial Resistance Reference Gene Database of acquired antimicrobial resistance genes using the `AMRFinderPlus tool <https://www.nature.com/articles/s41598-021-91456-0/>`_ . Identified AMR determinants and point mutations are grouped by drug class.


AMR parameters
++++++++++++++++++

``--organism`` 

To screen for point mutations

``-t , --threads`` 

Number of threads for alignment


AMR outputs
++++++++++++++++++

Results of the *Escherichia* AMR module are grouped by drug class 

.. list-table::
   :header-rows: 0

   * - Aminoglycoside
     - Aminoglycoside resistance genes

   * - Fluoroquinolone
     - Fluoroquinolone resistance genes

   * - Fosfomycin
     - Fosfomycin resistance genes

   * - Sulfonamide
     - Sulfonamide resistance genes

   * - Tetracycline
     - Tetracycline resistance genes

   * - Glycopeptide
     - Glycopeptide resistance genes

   * - Colistin
     - colistin resistance genes

   * - Phenicol
     - phenicol resistance genes

   * - Macrolide
     - Macrolide resistance genes

   * - Rifamycin
     - Rifampin resistance genes

   * - Trimethoprim
     - Trimethoprim resistance genes

   * - BetaLactam
     - Beta-lactamases

   * - Carbapenem
     - Carbapenemases

   * - Cephalosporin
     - Third-generation Cephalosporin

   * - Methicillin
     - Methicillin resistance genes

   * - Other Classes
     - Other categories
