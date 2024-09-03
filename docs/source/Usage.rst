
########################
Usage
########################

Input files
-----------

Genome assemblies in FASTA format (can be gzipped). 

Can be either draft or completed assemblies (completed is better because it reduces the risk of fragmented genes/loci).


Basic usage
-----------

Run with preset modules for *K. pneumoniae* species complex (KpSC), which reproduces the behaviour of Kleborate v2:

.. code-block:: Python

   kleborate -a *.fasta.gz -o kleborate_results -p kpsc --trim_headers

- ``-a *.fasta.gz``: Specifies the input files (assemblies) to be analyzed (.fasta or .fasta.gz).
- ``-o``: Specifies the directory where the output files will be saved (one output file per species/complex detected).
- ``-p``: Specifies the preset modules to run (kpsc, kosc, escherichia).
- ``--trim_headers``: Trim module names from column headers in the output.


Run with specified modules only, e.g. AMR typing for *K. pneumoniae* species complex:

.. code-block:: Python

   kleborate -a *.fasta -o kleborate_results -m klebsiella_pneumo_complex__amr

(A list of modules is available via ```kleborate --list_modules``` or `here <https://kleboratemodular.readthedocs.io/en/latest/modules.html>`_)


Run with preset modules for *K. oxytoca* species complex

.. code-block:: Python

   kleborate -a *.fasta -o kleborate_results -p kosc


Run with preset modules for *E. coli* or other *Escherichia*, on gzipped assemblies:

.. code-block:: Python

   kleborate  -a *.fasta.gz -o kleborate_results -p escherichia


Check available modules, check version, print help:

.. code-block:: Python

   kleborate [--list_modules] [--version] [-h]



Output files
--------------------

Output files are tab-delimited (.txt) files, one per species/complex detected, named in the format: klebsiella_pneumo_complex_output.txt.
Columns included in each output file will depend on the modules that are run; essentially each module creates a set of results columns that are added to the output file for the relevant species/complex. By default, each column name is preprended with the name of the module that generated it. This can be turned off using --trim_headers when running kleborate, or these column headers can be stripped off later using the trim_headers.py script.


Parameters
----------

**Input/output:**

``-a ASSEMBLIES [ASSEMBLIES ...], --assemblies ASSEMBLIES [ASSEMBLIES ...]``
    FASTA file(s) for assemblies, optionally gzipped (.gz)

``-o OUTDIR, --outdir OUTDIR``
    Directory for storing output files (default: Kleborate_results)

``--trim_headers``
    Trim headers in the output files (switch on to remove module names from the column headers in the output files). Alternatively, users can trim the headers off later using this script: `trim_headers.py <https://github.com/klebgenomics/KleborateModular/blob/main/kleborate/shared/trim_headers.py>`_

**Modules:**

``-p PRESET, --preset PRESET``         
    Module presets, choose from:

.. list-table::

   * - kpsc
     - *Klebsiella pnuemoniae* species complex

   * - kosc
     - *Klebsiella oxytoca* species complex
                                        
   * - escherichia 
     - *Escherichia* genus

``--list_modules``         
    Print a list of all available modules and then quit (default: False)

``-m MODULES, --modules MODULES``         
    Comma-delimited list of Kleborate modules to use

Module-specific parameters can be found `here <https://kleboratemodular.readthedocs.io/en/latest/modules.html>`_


**Help:**
     
``-h, --help``       
    Show a help message and exit

``--help_all``         
    Show a help message with all module options

``--version``         
    Show program's version number and exit





