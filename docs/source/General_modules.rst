
########################
General modules
########################

.. _species_detection:

Species detection
-----------------

``-m enterobacterales__species``


This module will attempt to identify the species of each input assembly. It does this by comparing the assembly using `Mash <https://mash.readthedocs.io/>`_ to a curated set of *Klebsiella* and other *Enterobacteriaceae* assemblies from NCBI, and reporting the species of the closest match. 

Parameters
++++++++++++++++++

``--enterobacterales__species_strong``

Mash distance threshold for a strong species match (default: 0.02)

``--enterobacterales__species_weak``

Mash distance threshold for a weak species match (default: 0.04)

Outputs
+++++++

Output of the species typing module is the following columns:

.. list-table::

   * - species
     - Species name (scientific name)

   * - species_match
     - Strength of the species call indicated as ``strong``\  (Mash distance ≤ 0.02) or ``weak``\  (Mash distance of > 0.02 and ≤ 0.04, may be novel or hybrid species)

The quality and completeness of Kleborate results depends on the quality of the input genome assemblies. In general, you can expect good results from draft genomes assembled with tools like SPAdes from high-depth (>50x) Illumina data, however it is always possible that key genes subject to genotyping may be split across contigs, which can create problems for detecting and typing them accurately.


Contig stats
------------

.. _contig_stats:

.. code-block:: Python

   -m general__contig_stats

This module takes ``enterobacterales__species`` as a prerequisite and  generates some basic assembly statistics to help users understand their typing results in the context of assembly quality, although we recommend users conduct more comprehensive QC themselves before typing genomes (e.g. screen for contamination, etc).

The module reports a standard set of assembly quality metrics (see Outputs below).


It will also flag in the ``QC_warnings``\  column if an assembly size falls outside those specified in the ``species_specification.txt``\  in the module directory, or if N50 <10 kbp or ambiguous bases (Ns) are detected in the sequence.

Outputs
+++++++

Output of the contig stats module is the following columns:

.. list-table::

   * - contig_count
     - Number of contigs in the input assembly

   * - N50
     - `N50 <https://en.wikipedia.org/wiki/N50,_L50,_and_related_statistics>`_ calculated from the contig sizes

   * - largest_contig
     - Size of largest contig (in bp)

   * - total_size
     - Total assembly size (in bp)

   * - ambiguous_bases
     - Detection of ambiguous bases (yes or no). If yes, the number of ambiguous bases is also provided in brackets.

   * - QC_warnings
     - List of QC issues detected, including: ``ambiguous_bases``\ (ambiguous bases detected) ``N50``\ (N50 < 10 kbp), ``total_size`` (genome size falls outside expected range).
