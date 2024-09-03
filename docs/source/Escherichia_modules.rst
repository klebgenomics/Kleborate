
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


E. coli MLST parameters
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

E. coli MLST outputs
++++++++++++++++++++++

Output of the Pasteur E. coli MLST module is the following columns:

.. list-table::

   * - ST
     - sequence type

   * - dinB, icdA, pabB, polB, putP, trpA, trpB, uidA
     - allele number

Output of the Achtman E. coli MLST module is the following columns:

.. list-table::

   * - ST
     - sequence type

   * - adk, fumC, gyrB, icd, mdh, purA, recA
     - allele number

* Kleborate makes an effort to report the closest matching ST if a precise match is not found.
* Imprecise allele matches are indicated with a ``*``.
* Imprecise ST calls are indicated with ``-nLV``\ , where n indicates the number of loci that disagree with the ST reported. So ``258-1LV`` indicates a single-locus variant (SLV) of ST258, i.e. 6/7 loci match ST258.
