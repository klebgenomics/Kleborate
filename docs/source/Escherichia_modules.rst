
****************************************************
Modules for *Escherichia* species
****************************************************


.. code-block:: Python

   --preset escherichia

These modules will be deployed if the ``enterobacterales__species``\   module confirms the input assembly as a member of the *Escherichia* genus. 

.. _escherichia__mlst_achtman:

.. _escherichia__mlst_pasteur:

MLST
------------

.. code-block:: Python

   -m escherichia__mlst_achtman, escherichia__mlst_pasteur

The Achtman scheme is hosted on `EnteroBase <https://enterobase.warwick.ac.uk/>`_ and described in `this paper <https://doi.org/10.1111/j.1365-2958.2006.05172.x>`_.

We also provide an option for users to run MLST using the Pasteur scheme by running:
``-m escherichia__mlst_pasteur``

The Pasteur scheme is hosted in the BIGSdb-Pasteur *Escherichia coli* `database <https://bigsdb.pasteur.fr/ecoli/>`_ and described in `this paper <https://doi.org/10.1186/1471-2164-9-560>`_.

The genes included in each scheme are noted in the Outputs table below.

A copy of the MLST alleles and ST definitions used in each module is stored in the ``data/``  directory of the module.


Parameters
+++++++++++++++++++++++++++

``--escherichia_mlst_achtman_min_identity`` 

Minimum alignment percent identity for *Escherichia-Achtman* MLST (default: 90.0)

``--escherichia_mlst_achtman_min_coverage`` 

Minimum alignment percent coverage for Escherichia-Achtman MLST (default: 80.0)

``--escherichia_mlst_achtman_required_exact_matches`` 

Minimum number of exact allele matches required to call an ST (default: 3).

``--escherichia_mlst_pasteur_min_identity`` 

Minimum alignment percent identity for Escherichia-Pasteur MLST (default: 90.0)

``--escherichia_mlst_pasteur_min_coverage`` 

Minimum alignment percent coverage for Escherichia-Pasteur MLST (default: 80.0)

``--escherichia_mlst_pasteur_required_exact_matches`` 

Minimum number of exact allele matches required to call an ST (default: 4).

Outputs
++++++++++++++++++++++

Output of the Achtman *E. coli* MLST module include the following columns:

.. list-table::
   :header-rows: 0

   * - **Column name**
     - **Content**

   * - ``ST_Achtman``
     - Sequence type

   * - ``adk``, ``fumC``, ``gyrB``, ``icd``, ``mdh``, ``purA``, ``recA``
     - Allele numbers for the Achtman scheme loci


Output of the Pasteur *E. coli* MLST module includes the following columns:

.. list-table::
   :header-rows: 0

   * - **Column name**
     - **Content**

   * - ``ST_Pasteur``
     - Sequence type.

   * - ``dinB``, ``icdA``, ``pabB``, ``polB``, ``putP``, ``trpA``, ``trpB``, ``uidA``
     - Allele numbers for the Pasteur scheme loci.

Notes
+++++

* Kleborate attempts to report the closest matching ST if a precise match is not found.
* Imprecise allele matches are indicated with a ``*``.
* Imprecise ST calls are indicated with ``-nLV``\ , where n indicates the number of loci that differ from the ST reported. For example, ``131-1LV`` indicates a single-locus variant (SLV) of ST131, i.e. 6/7 loci match ST131.


.. _escherichia__cgMLST:

cgMLST and LIN codes
------------------

.. code-block:: Python

   -m escherichia__cgMLST

This module performs cgMLST and  LIN code typing using `MiST <https://github.com/BioinformaticsPlatformWIV-ISP/MiST>`_  tool. 


Outputs
++++++++++++++++++
cgMLST results are output in the following columns:

.. list-table::
   :header-rows: 1

   * - **Column name**
     - **Content**

   * - cgST
     - Best matching scgST

   * - LIN code
     - LIN code / Partial LINcode for input strain

   * - Sublineage
     - Sublineage for input strain

   * - Clonal group
     - Clonal group for input strain


.. _escherichia__ezclermont:


Phylgroups
----------------------

.. code-block:: Python

   -m escherichia__ezclermont


The *Escherichia* genus comprises several clades, including *Escherichia albertii*, *E. fergusonii*, five cryptic *Escherichia* clades (I–V) and *E. coli* sensu stricto. Within *E. coli*, strains can be further divided into seven main phylogroups: A, B1, B2, C, D, E and F. 

Kleborate assigns Escherichia genomes to these clades and phylogroups using `EzClermont tool <https://doi.org/10.1099/acmi.0.000143>`_, which is based on *in vitro* PCR assay logic.


Parameters
++++++++++

``--escherichia__ezclermont_min_length``

Minimum contig length to consider. *Default:* ``500``


Outputs
+++++++

.. list-table:: 
   :header-rows: 0

   * - **Column name**
     - **Content**

   * - ``Clermont_type``
     - Assigned phylogroup or clade.

   * - ``Clermont_profile``
     - Presence or absence pattern of PCR products.



.. _escherichia__pathovar:

Pathotyping
---------------------

.. code-block:: Python

   -m escherichia__pathovar

*Escherichia coli* is broadly divided into 2 groups: intestinal diarrheagenic *E. coli* (DEC), and extra-intestinal *E. coli* (ExPEC). DEC encompasses several clinically relevant pathotypes: enteropathogenic *E. coli* (EPEC), enterotoxigenic *E. coli* (ETEC), enterohaemorrhagic *E. coli* (EHEC), Shiga toxin-producing *E. coli* (STEC), enteroaggregative *E. coli* (EAEC), enteroinvasive *E. coli* (EIEC), and diffusely adherent *E. coli* (DAEC) see `paper <https://doi.org/10.3389/fcimb.2016.00141>`_. Additionally, *Shigella* is considered a DEC pathotype due to its genetic and pathogenetic similarity to EIEC, see `paper <https://doi.org/10.1128/mbio.00882-23>`_.

The majority of DEC pathotypes are defined by specific virulence markers. However, for EAEC, DAEC and AIEC, the pathogenic role of proposed markers is not well established. 

This module classifies *E. coli* genomes into DEC pathotypes based on the presence or absence of key virulence marker genes, determined by alignment with Minimap2 to reference sequences for key marker genes, extracted from the `VirulenceFinder <https://cge.food.dtu.dk/services/VirulenceFinder/>`_ database.

Marker genes considered are:

.. list-table:: 
   :header-rows: 1

   * - **Gene**
     - **Marker of**

   * - *ipaH*
     - pINV plasmid shared by Shigella and EIEC

   * - *sta1* and *ltcA*
     - ST, heat-stable enterotoxin and LT, heat-labile enterotoxin, characteristic of ETEC
   
   * - *stx1A*, *stx1B*, *stx2A*, *stx2B*
     - shiga-toxins Stx1 and Stx2

   * - *bfpA*
     - *bfp* adhesin characteristic of EPEC

   * - *eae*
     - intimin, characteristic of EPEC and EHEC

The module also calls `ShigaPass <https://github.com/imanyass/ShigaPass>`_ to predict Shigella serotypes, and to differentiate *Shigella* from enteroinvasive *Escherichia coli* (EIEC).

The combination of marker genes detected are then used to classify genomes into pathotypes, using logic adapted from `EnteroBase <https://enterobase.readthedocs.io/en/latest/pipelines/backend-pipeline-phylotypes.html?highlight=pathovar/>`_ as follows:

.. list-table:: 
   :header-rows: 1

   * - **Pathotype**
     - **Marker/s Present**
     - **Marker/s Absent**

   * - EIEC
     - pINV
     - ST and LT, Shigella serotypes

   * - ETEC
     - ST or LT
     - pINV

   * - STEC
     - Stx1 or Stx2
     - intimin

   * - EHEC
     - Stx1 or Stx2, and intimin
     - 

   * - EPEC
     - intimin and bfp
     - Stx1 and Stx2

   * - aEPEC
     - intimin
     - Stx1 and Stx2 and bfp

   * - Shigella
     - Shigella serotypes
     - 


Parameters
++++++++++++++++++++++++++++++++++

``--escherichia__pathovar_min_identity``

Minimum alignment percent identity for pathotype (default: 90.0).

``--escherichia__pathovar_min_coverage``

Minimum alignment percent coverage for pathotype (default: 80.0).


Outputs
++++++++++++++++++++++++++++

.. list-table:: 
   :header-rows: 0

   * - **Column name**
     - **Content**

   * - ``Pathotype``
     - Predicted pathotype

   * - ``Stx1``, ``Stx2``, ``ST``, ``LT``, ``eae``, ``ipaH``
     - Virulence markers


Additionally, Kleborate includes a separate ``escherichia__vfdb`` module that screens *E. coli* genome assemblies against all markers in the `VirulenceFinder <https://cge.food.dtu.dk/services/VirulenceFinder/>`_ database (VFDB).


.. _escherichia__mlst_lee:

LEE pathogenicity island typing
----------------------------------------------

.. code-block:: Python

   -m escherichia__mlst_lee

Locus of enterocyte effacement (LEE) is a ~40 kb chromosomal pathogenicity island composed of 41 core genes organized into five operons  `Elliot et al., 1998 <https://onlinelibrary.wiley.com/doi/10.1046/j.1365-2958.1998.00783.x>`_. It encodes an (i) outer membrane adhesive protein, known as intimin protein and encoded by the _eae_ gene; (ii) type III secretion system (T3SS), and (iii) translocated receptor (Tir) as well as translocons, chaperones, regulators and secreted effector proteins that are linked to virulence.

Kleborate includes a module for subtyping of the LEE pathogenicity island. Details of the the scheme, and the LEE subtypes and lineages it identifies, can be found in this `paper <https://doi.org/10.1038/nmicrobiol.2015.10>`_.

The LEE typing database is based on analysis of >250 LEE-containing *E. coli* genomes and includes 7 loci (eae (intimin), tir, espA, espB, espD, espH, espZ). The data is provided as a MLST-style database, in which combinations of alleles are assigned to a LEE subtype, to facilitate a common nomenclature for LEE subtypes. Each sequence in the database represents a cluster of closely related alleles that have been assigned to the same locus type. The LEE scheme includes three distinct lineages: Lineage 1 consists of LEE subtypes 1-2; Lineage 2 consists of LEE subtypes 3-8; Lineage 3 consists of LEE subtypes 9-30.

The reference sequences and MLST-style profile definitions are included in the ``data/``  directory of this module.


Parameters
++++++++++

``--escherichia__mlst_LEE_min_identity``

Minimum alignment percent identity for ``escherichia__mlst_LEE``. *Default:* ``90.0``

``--escherichia__mlst_LEE_min_coverage``

Minimum alignment percent coverage for ``escherichia_mlst_LEE``. *Default:* ``80.0``

``escherichia__mlst_LEE_mlst_required_exact_matches``

Minimum number of exact allele matches required to assign an ST. *Default:* ``3``


Outputs
++++++++++++++++++++++++++++

The output of the *E. coli* LEE MLST module includes the following columns:

.. list-table::

   * - **Column name**
     - **Content**

   * - ``LEE_ST``
     - Assigned LEE sequence type.

   * - ``LEE_lineage``
     - Lineage associated with the LEE ST.

   * - ``LEE_eae``, ``LEE_tir``, ``LEE_espA``, ``LEE_espB``, ``LEE_espD``, ``LEE_espH``, ``LEE_espZ``

     - Allele numbers for each LEE locus.

Additional Notes
++++++++++++++++

* Kleborate attempts to report the closest matching ST if an exact match is not found.
* Imprecise allele matches are indicated with a ``*``.
* Imprecise ST calls are indicated with ``-nLV``\ , where n indicates the number of loci that disagree with the ST reported. For example, ``ST10-3LV`` indicates a three-locus variant (SLV) of ST10 (i.e. 4/7 loci match ST10).


.. _escherichia__stxtyper:

Stx typing
-----------

.. code-block:: Python

   -m escherichia__stxtyper

Shiga toxins (Stxs) are key virulence factors of Stx-producing *Escherichia coli* (STEC). They are also found in *Shigella dysenteriae 1*. Stxs belong to the AB-type toxin family and are divided into two antigenically distinct groups: Stx1 and Stx2. Each group contains several variants/subtypes—six for Stx1 (a, b, c, d, e, f) and seven for Stx2 (a, b, c, d, e, f, and g) (`Yano et al., 2023 <https://doi.org/10.1038/s41598-023-32111-8>`_, `Melton-Celsa, 2014 <https://doi.org/10.1128/microbiolspec.EHEC-0024-2013>`_). These toxins are encoded by lysogenic bacteriophages (Stx phage) and STEC strains may produce either single Stx subtype or a combination of subtypes.

This module will run `StxTyper <https://doi.org/10.3390/microorganisms14081607>` to determine the _stx_ type. See the `StxTyper documentation <https://github.com/ncbi/stxtyper>`_ for more details of how it works.


Outputs
+++++++++++++++++++++

StxTyper results are output in the following columns:

.. list-table::
   :header-rows: 1

   * - **Column name**
     - **Content**

   * - ``Stx_type``
     - The Shiga toxin type. If the operon is complete, the subtype will be reported (e.g., ``stx1a``). If the operon is incomplete or ambiguous, a broader designation is used: ``stx1``, ``stx2``, or simply ``stx`` if the algorithm cannot resolve at further.
   * - ``operon``
     - Status the operon detected. Possible values:
       ``COMPLETE`` – Full operon found.  
       ``PARTIAL`` – Operon incomplete.  
       ``PARTIAL_CONTIG_END`` – Partial operon likely truncated at contig boundary.  
       ``EXTENDED`` – Coding sequence extends beyond the reference stop codon for one or both subunits.  
       ``INTERNAL_STOP`` – A subunit contains a nonsense mutation.  
       ``FRAMESHIFT`` – Indel detected in coding sequence.  
       ``AMBIGUOUS`` – Ambiguous base(s) found in the sequence.  
       ``COMPLETE_NOVEL`` – Full-length operon that cannot be typed.
   * - ``identity``
     - Percent identity for both A and B subunits.
   * - ``target_start``
     - Start position of the alignment.
   * - ``target_stop``
     - End position of the alignment.
   * - ``target_strand``
     - Strand orientation of the target sequence.
   * - ``A_reference``
     - Closest reference protein for the A subunit.
   * - ``A_identity``
     - Percent identity to the reference for the A subunit.
   * - ``A_reference_subtype``
     - Subtype assigned to the reference sequence for the A subunit.
   * - ``A_coverage``
     - Percentage of the A subunit reference sequence covered by the alignment.
   * - ``B_reference``
     - Closest reference protein for the B subunit.
   * - ``B_reference_subtype``
     - Subtype assigned to the reference sequence for the B subunit.
   * - ``B_identity``
     - Percent identity to the reference for the B subunit.
   * - ``B_coverage``
     - Percentage of the B subunit reference sequence covered by the alignment.

.. _escherichia__pks:

Polyketide synthetase (colibactin) detection 
----------------------------------------------

.. code-block:: bash

  -m escherichia__pks

This module screens for the presence or absence of *clbB*, gene in the colibactin polyketide synthase (PKS) biosynthesis gene cluster, by aligning genome assemblies against the *clbB* `reference sequence <https://www.ncbi.nlm.nih.gov/nuccore/AM229678.1?from=41762&to=51382>`_.

Parameters
++++++++++

``--escherichia__pks_min_identity``

Minimum alignment percent identity for detecting *clbB* (default: 90.0)

``--escherichia__pks_min_coverage``

Minimum alignment percent coverage for detecting *clbB* (default: 80.0)


Outputs
+++++++++++

.. list-table::

   * - **Column name**
     - **Content**

   * - clbB
     - presence or absence of the gene


.. _escherichia__ectyper:

*E. coli* O:H serotyping
----------------------

.. code-block:: Python

   -m escherichia__ectyper

*E. coli* serotypes are defined by combinations of O (lipopolysaccharide) and H (flagellar) antigens. Currently there are ~183 O-groups and 53 H-types that have been defined serologically `Ørskov and Ørskov 1984 <https://www.sciencedirect.com/science/article/abs/pii/S0580951708704471/>`_.

**O-antigen**

The O-antigen is an integral component of the Lipopolysaccharide (LPS) found in the outer membrane of the bacteria. LPS comprises three components: lipid A, a core oligosaccharide, and the O-specific polysaccharide chain (O antigen).  The O-antigen domain exhibits significant variability consisting of 10 to 25 repeating oligosaccharide units, with each unit containing two to seven sugar residues `Liu et al., 2020 <https://pmc.ncbi.nlm.nih.gov/articles/PMC7685785/>`_. The genes responsible for synthesis of O-antigens are usually present as a gene cluster and are located between the two chromosomal housekeeping genes galF and gnd/ugd `Iguchi et al 2014 <https://pmc.ncbi.nlm.nih.gov/articles/PMC4379981/>`_. Major pathways involved in the assembly, synthesis and transport of O-antigen include, the Wzy pathway the Wzx/Wzy-dependent pathway, encoded by the wzx (O-antigen flippase) and wzy (O-antigen polymerase) genes, and the ABC transporter pathway, encoded by wzm and wzt. These genes are ideal biomarkers for predicting O antigen types.  


**H antigens**

H antigens (flagellar) are surface proteins composed of repeated molecules of the protein flagellin, which facilitate bacterial motility. These antigens are numbered from H1 to H56 (H13, H22, and H50 are not used) and are distinct from the O and K antigens. Flagellin is encoded by the fliC gene on the chromosomal locus or its homologues (non-fliC flagellin-coding genes such as flkA, fllA, and flmA). Of the 53 well known H antigen types, 44 are conferred by expression of the fliC gene,  the remaining 9 H types are  encoded by non-fliC flagellin genes. Specifically H3, H35, H36, H47,and H53 are encoded by flkA, H44 and H55 by fllA, H54 by flmA, and H17 by flnA.

Kleborate uses ECTyper for in silico serotyping. See `ECTyper paper <https://pmc.ncbi.nlm.nih.gov/articles/PMC8767331/>`_ for more details.

Outputs
+++++++

Outputs of the ECTyper module is the following columns:

.. list-table:: 
   :header-rows: 0

   * - **Column name**
     - **Content**

   * - ``O-type``
     - Predicted O antigen.

   * - ``H-type``
     - Predicted H antigen.

   * - ``Serotype``
     - Combined prediction of O and H antigens.

   * - ``QC``
     - Quality control values summarising the overall confidence of the serotype prediction.

   * - ``Evidence``
     - Total number of alleles used to call both O and H antigens.

   * - ``GeneScores``
     - ECTyper gene scores for O and H antigens, ranging from 0 to 1.

   * - ``AllelesKeys``
     - Best-matching allele keys from the ECTyper database used for serotype assignment.

   * - ``GeneIdentities(%)``
     - Percent identity values of the query alleles.

   * - ``GeneCoverages(%)``
     - Percent coverage values for the query alleles.

   * - ``GeneLengths``
     - Gene lengths ( in base pairs) of the query alleles.

   * - ``Warnings``
     - Additional messages related to QC status or other issues affecting serotype prediction.


.. _escherichia__kaptive:

Capsule (K) typing
--------------

.. code-block:: Python

   -m escherichia__kaptive

This module will run the `Kaptive <https://github.com/klebgenomics/kaptive>`_ v3 tool to type *E. coli* group 2 and group 3 capsular (K) loci. See the `capsular K-typing database <https://github.com/rgladstone/EC-K-typing>`_ for more details.


Parameters
+++++++++++++++++++

``--ecoli_kps``

Group 2 + 3 CPS database 


Outputs
+++++++++++++++++

Kaptive results are output in the following columns:

.. list-table::
   :header-rows: 1

   * - **Column name**
     - **Content**
   * - K_locus
     - The locus type which most closely matches the assembly.
   * - K_type
     - The predicted serotype/phenotype of the assembly.
   * - K_locus_confidence
     - Typeable or Untypeable.
   * - K_locus_problems
     - Characters indicating issues with the locus match (see problems).
   * - K_locus_identity
     - Weighted percent identity of the best matching locus to the assembly.
   * - K_Missing_expected_genes
     - Missing expected genes in the locus



.. _escherichia__amr:


Antimicrobial resistance
------------------------

.. code-block:: Python

   -m escherichia__amr

This module screens input genomes for acquired antimicrobial resistance genes and known resistance-associated point mutations using the `AMRFinderPlus tool <https://doi.org/10.1038/s41598-021-91456-0>`_. Identified determinants are grouped by drug class.


Parameters
++++++++++++++++++

``--organism``

Used to screen for point mutations in species-specific resistance markers (default Escherichia coli).

``-t , --threads`` 

Number of threads to use for alignment.


Outputs
++++++++++++++++++

Results of the *Escherichia* AMR module are grouped by drug class:

.. list-table::
   :header-rows: 0

   * - **Column name**
     - **Markers assigned as NCBI class/subclass**

   * - ``Aminoglycoside``
     - AMINOGLYCOSIDE

   * - ``Fluoroquinolone``
     - FLUOROQUINOLONE, QUINOLONE

   * - ``Fosfomycin``
     - FOSFOMYCIN

   * - ``Sulfonamide``
     - SULFONAMIDE

   * - ``Tetracycline``
     - TETRACYCLINE

   * - ``Glycopeptide``
     - GLYCOPEPTIDE

   * - ``Colistin``
     - COLISTIN

   * - ``Phenicol``
     - PHENICOL

   * - ``Macrolide``
     - MACROLIDE

   * - ``Rifamycin``
     - RIFAMYCIN

   * - ``Trimethoprim``
     - TRIMETHOPRIM

   * - ``BetaLactam``
     - BETA-LACTAM, CEPHALOTHIN

   * - ``Carbapenem``
     - CARBAPENEM

   * - ``Cephalosporin``
     - CEPHALOSPORIN (note these can be assumed to be 3rd generation cephalosporins)

   * - ``Methicillin``
     - METHICILLIN

   * - ``Other Classes``
     - other CLASS values
