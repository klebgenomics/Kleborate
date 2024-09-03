
########################
Usage
########################

Input files
-----------

KleborateModular takes genome assemblies in FASTA format (can be gzipped). It will work on either draft or completed assemblies, though completed is better because it reduces the risk of fragmented genes/loci.

If you have unassembled reads, you will need to assemble them before running Kleborate. You can try our `Unicycler <https://github.com/rrwick/Unicycler>`_ assembler which works great on Illumina or hybrid Illumina + Nanopore/PacBio reads. Or check out this `guide to bacterial genome assembly <https://github.com/rrwick/Trycycler/wiki/Guide-to-bacterial-genome-assembly>`_.

Basic usage
------------------

Run with preset modules for *K. pneumoniae* species complex (KpSC):

.. code-block:: Python

   ./kleborate-runner.py -a *.fasta -o results.txt -p kpsc

Run with specified modules only:

.. code-block:: Python

   ./kleborate-runner.py -a *.fasta -o results.txt -m 

Run with preset modules for *E. coli* or other *Escherichia*, on gzipped assemblies:

.. code-block:: Python

   ./kleborate-runner.py  -a *.fasta.gz -o results.txt -p escherichia

Check available modules, check version, print help:

.. code-block:: Python

   ./kleborate-runner.py [--list_modules] [--version] [-h]


Input/output:
  -a ASSEMBLIES [ASSEMBLIES ...], --assemblies ASSEMBLIES [ASSEMBLIES ...]
                                        FASTA file(s) for assemblies

  -o OUTFILE, --outfile OUTFILE         
                                        File for detailed output (default: Kleborate_results.txt)

Modules:
  --list_modules         
                                        Print a list of all available modules and then quit (default: False)

  -p PRESET, --preset PRESET         
                                        Module presets, choose from 

                                        - kpsc (*Klebsiella pnuemoniae* species complex)
                                        - kosc (*Klebsiella oxytoca* species complex)
                                        - escherichia  (*Escherichia*)


  -m MODULES, --modules MODULES         
                                        Comma-delimited list of Kleborate modules to use


Help:
     -h, --help         
                                        Show a help message and exit
     --help_all         
                                        Show a help message with all module options
     --version         
                                        Show program's version number and exit


Available modules and parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Generic modules
+++++++++++++++

enterobacterales__species module:
  --enterobacterales__species_strong ENTEROBACTERALES__SPECIES_STRONG
                                        Mash distance threshold for a strong species match (default: 0.02)
  --enterobacterales__species_weak ENTEROBACTERALES__SPECIES_WEAK
                                        Mash distance threshold for a weak species match (default: 0.04)


Klebsiella pneumoniae SC modules
++++++++++++++++++++++++++++++++

klebsiella__abst module:
  --klebsiella__abst_min_identity KLEBSIELLA__ABST_MIN_IDENTITY
                                        Minimum alignment percent identity for aerobactin MLST (default: 90.0)
  --klebsiella__abst_min_coverage KLEBSIELLA__ABST_MIN_COVERAGE
                                        Minimum alignment percent coverage for aerobactin MLST (default: 80.0)
  --klebsiella__abst_required_exact_matches KLEBSIELLA__ABST_REQUIRED_EXACT_MATCHES
                                        At least this many exact matches are required to call an ST (default: 3)

klebsiella__cbst module:
  --klebsiella__cbst_min_identity KLEBSIELLA__CBST_MIN_IDENTITY
                                        Minimum alignment percent identity for colibactin MLST (default: 90.0)
  --klebsiella__cbst_min_coverage KLEBSIELLA__CBST_MIN_COVERAGE
                                        Minimum alignment percent coverage for colibactin MLST (default: 80.0)
  --klebsiella__cbst_required_exact_matches KLEBSIELLA__CBST_REQUIRED_EXACT_MATCHES
                                        At least this many exact matches are required to call an ST (default: 8)

klebsiella__rmst module:
  --klebsiella__rmst_min_identity KLEBSIELLA__RMST_MIN_IDENTITY
                                        Minimum alignment percent identity for Rmp MLST (default: 90.0)
  --klebsiella__rmst_min_coverage KLEBSIELLA__RMST_MIN_COVERAGE
                                        Minimum alignment percent coverage for Rmp MLST (default: 80.0)
  --klebsiella__rmst_required_exact_matches KLEBSIELLA__RMST_REQUIRED_EXACT_MATCHES
                                        At least this many exact matches are required to call an ST (default: 2)

klebsiella__smst module:
  --klebsiella__smst_min_identity KLEBSIELLA__SMST_MIN_IDENTITY
                                        Minimum alignment percent identity for salmochelin MLST (default: 90.0)
  --klebsiella__smst_min_coverage KLEBSIELLA__SMST_MIN_COVERAGE
                                        Minimum alignment percent coverage for salmochelin MLST (default: 80.0)
  --klebsiella__smst_required_exact_matches KLEBSIELLA__SMST_REQUIRED_EXACT_MATCHES
                                        At least this many exact matches are required to call an ST (default: 2)

klebsiella__ybst module:
  --klebsiella__ybst_min_identity KLEBSIELLA__YBST_MIN_IDENTITY
                                        Minimum alignment percent identity for yersiniabactin MLST (default: 90.0)
  --klebsiella__ybst_min_coverage KLEBSIELLA__YBST_MIN_COVERAGE
                                        Minimum alignment percent coverage for yersiniabactin MLST (default: 80.0)
  --klebsiella__ybst_required_exact_matches KLEBSIELLA__YBST_REQUIRED_EXACT_MATCHES
                                        At least this many exact matches are required to call an ST (default: 2)

klebsiella_pneumo_complex__mlst module:
  --klebsiella_pneumo_complex__mlst_min_identity KLEBSIELLA_PNEUMO_COMPLEX__MLST_MIN_IDENTITY
                                        Minimum alignment percent identity for klebsiella_pneumo_complex_MLST (default: 90.0)
  --klebsiella_pneumo_complex__mlst_min_coverage KLEBSIELLA_PNEUMO_COMPLEX__MLST_MIN_COVERAGE
                                        Minimum alignment percent coverage for klebsiella_pneumo_complex_MLST (default: 80.0)
  --klebsiella_pneumo_complex__mlst_required_exact_matches KLEBSIELLA_PNEUMO_COMPLEX__MLST_REQUIRED_EXACT_MATCHES
                                        At least this many exact matches are required to call an ST (default: 3)

klebsiella_pneumo_complex__amr module:
  --klebsiella_pneumo_complex__amr_min_identity KLEBSIELLA_PNEUMO_COMPLEX__AMR_MIN_IDENTITY
                                        Minimum alignment percent identity for klebsiella_pneumo_complex Amr results (default: 90.0)
  --klebsiella_pneumo_complex__amr_min_coverage KLEBSIELLA_PNEUMO_COMPLEX__AMR_MIN_COVERAGE
                                        Minimum alignment percent coverage for klebsiella_pneumo_complex Amr results (default: 80.0)
  --klebsiella_pneumo_complex__amr_min_spurious_identity KLEBSIELLA_PNEUMO_COMPLEX__AMR_MIN_SPURIOUS_IDENTITY
                                        Minimum alignment percent identity for klebsiella_pneumo_complex Amr spurious results (default: 80.0)
  --klebsiella_pneumo_complex__amr_min_spurious_coverage KLEBSIELLA_PNEUMO_COMPLEX__AMR_MIN_SPURIOUS_COVERAGE
                                        Minimum alignment percent coverage for klebsiella_pneumo_complex Amr spurious results (default: 40.0)

klebsiella_pneumo_complex__kaptive module:
  -t , --threads                        Number of threads for alignment (default: 1)



Klebsiella oxytoca SC modules
+++++++++++++++++++++++++++++

klebsiella_oxytoca_complex__mlst module:
  --klebsiella_oxytoca_complex__mlst_min_identity KLEBSIELLA_OXYTOCA_COMPLEX__MLST_MIN_IDENTITY
                                        Minimum alignment percent identity for klebsiella_oxytoca_complex MLST (default: 90.0)
  --klebsiella_oxytoca_complex__mlst_min_coverage KLEBSIELLA_OXYTOCA_COMPLEX__MLST_MIN_COVERAGE
                                        Minimum alignment percent coverage for klebsiella_oxytoca_complex MLST (default: 80.0)
  --klebsiella_oxytoca_complex__mlst_required_exact_matches KLEBSIELLA_OXYTOCA_COMPLEX__MLST_REQUIRED_EXACT_MATCHES
                                        At least this many exact matches are required to call an ST (default: 3)


Eschericia SC modules
+++++++++++++++++++++

escherichia__mlst_achtman module:
  --escherichia_mlst_achtman_min_identity ESCHERICHIA_MLST_ACHTMAN_MIN_IDENTITY
                                        Minimum alignment percent identity for *Escherchia-Achtman* MLST (default: 90.0)
  --escherichia_mlst_achtman_min_coverage ESCHERICHIA_MLST_ACHTMAN_MIN_COVERAGE
                                        Minimum alignment percent coverage for Escherchia-Achtman MLST (default: 80.0)
  --escherichia_mlst_achtman_required_exact_matches ESCHERICHIA_MLST_ACHTMAN_REQUIRED_EXACT_MATCHES
                                        At least this many exact matches are required to call an ST (default: 3)

escherichia__mlst_pasteur module:
  --escherichia_mlst_pasteur_min_identity ESCHERICHIA_MLST_PASTEUR_MIN_IDENTITY
                                        Minimum alignment percent identity for Escherchia-Pasteur MLST (default: 90.0)
  --escherichia_mlst_pasteur_min_coverage ESCHERICHIA_MLST_PASTEUR_MIN_COVERAGE
                                        Minimum alignment percent coverage for Escherchia-Pasteur MLST (default: 80.0)
  --escherichia_mlst_pasteur_required_exact_matches ESCHERICHIA_MLST_PASTEUR_REQUIRED_EXACT_MATCHES
                                        At least this many exact matches are required to call an ST (default: 4)
