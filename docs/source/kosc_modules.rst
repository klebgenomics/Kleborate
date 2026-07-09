****************************************************
Modules for *Klebsiella oxytoca* species complex
****************************************************

.. _klebsiella_oxytoca_complex__mlst:

.. code-block:: Python

   --preset kosc

These modules will be run if the ``enterobacterales__species``\   module confirms the input assembly as a member of the *K. oxytoca* species complex (KpSC) labelled in the tree above, and listed in the table below. 

.. list-table::

   * - Klebsiella oxytoca

   * - Klebsiella grimontii

   * - Klebsiella michiganensis

   * - Klebsiella pasteurii

   * - Klebsiella huaxiensis

   * - Klebsiella spallanzanii



KoSC MLST
---------

.. code-block:: Python

   -m klebsiella_oxytoca_complex__mlst

Genomes identified as belonging to the *K. oxytoca* species complex are subjected to MLST using the 7-locus scheme described at the  *K. oxytoca* `\at PubMLST <https://pubmlst.org/organisms/klebsiella-oxytoca>`_.

A copy of the MLST alleles and ST definitions is stored in the /data directory of this module.


KoSC MLST parameters
+++++++++++++++++++++++

``--klebsiella_oxytoca_complex__mlst_min_identity`` 

Minimum alignment percent identity for klebsiella_oxytoca_complex MLST (default: 90.0)

``--klebsiella_oxytoca_complex__mlst_min_coverage`` 

Minimum alignment percent coverage for klebsiella_oxytoca_complex MLST (default: 80.0)

``--klebsiella_oxytoca_complex__mlst_required_exact_matches``

At least this many exact matches are required to call an ST (default: 3)


KoSC MLST outputs
++++++++++++++++++++

Output of the KoSC MLST module is the following columns:

.. list-table::

   * - ST
     - sequence type

   * - gapA, infB, mdh, pgi, phoE, rpoB, tonB
     - allele number

* Kleborate makes an effort to report the closest matching ST if a precise match is not found.
* Imprecise allele matches are indicated with a ``*``.
* Imprecise ST calls are indicated with ``-nLV``\ , where n indicates the number of loci that disagree with the ST reported. So ``258-1LV`` indicates a single-locus variant (SLV) of ST258, i.e. 6/7 loci match ST258.


KoSC virulence typing
---------------------
Genomes identified as belonging to the *K. oxytoca* species complex are subjected to typing for yersiniabactin, colibactin, aerobactin, salmochelin and rmp loci, using the following modules:

.. code-block:: Python

   -m klebsiella__ybst, klebsiella__cbst, klebsiella__abst, klebsiella__smst, klebsiella__rmst

These modules were primarily designed for typing of *K. pneumoniae* species complex, and the databases are populated from variation detected in KpSC genomes. However they can appear in KoSC genomes and so typing is included in the KoSC preset.



KoSC K and O locus typing with Kaptive
-----------------------------------------

.. code-block:: Python

   -m kosc_kaptive

This module will run the `Kaptive <https://github.com/klebgenomics/kaptive>`_ v3 tool to identify capsule (K) and O antigen loci. See the Kaptive `documentation <https://kaptive.readthedocs.io/en/latest/>`_ for more details of how Kaptive works, tutorials, and citations.


Kaptive parameters
+++++++++++++++++++

``-t , --threads``

Number of threads for alignment (default: 1)

``--k-db, KoSC_K_locus_database.gbk``

Kaptive database for K-locus typing

``--o-db, KoSC_O_locus_database.gbk``

Kaptive database for o-locus typing



Kaptive outputs
+++++++++++++++++

Kaptive results are output in the following columns:

.. list-table::
   :header-rows: 1

   * - Column Name
     - Description
   * - Best match locus
     - The locus type which most closely matches the assembly.
   * - Best match type
     - The predicted serotype/phenotype of the assembly.
   * - Match confidence
     - Typeable or Untypeable.
   * - Problems
     - Characters indicating issues with the locus match (see problems).
   * - Identity
     - Weighted percent identity of the best matching locus to the assembly.
   * - Coverage
     - Weighted percent coverage of the best matching locus in the assembly.
   * - Length discrepancy
     - If the locus was found in a single piece, this is the difference between the locus length and the assembly length.
   * - Expected genes in locus
     - A fraction indicating how many of the genes in the best matching locus were found in the locus part of the assembly.
   * - Expected genes in locus, details
     - Gene names for the expected genes found in the locus part of the assembly.
   * - Missing expected genes
     - A string listing the gene names of expected genes that were not found.
