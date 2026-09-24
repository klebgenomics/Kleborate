########################
Installation
########################

Dependencies
=============
Kleborate requires the following software and libraries to be installed and available in your path:


* `Python <https://www.python.org/>`_ 
* `rammappy <https://github.com/tomdstanton/rammappy>`_ 
* `Biopython <https://biopython.org/>`_ 
* `Mash <https://github.com/marbl/Mash>`_ 
* `Minimap2 <https://github.com/lh3/minimap2>`_ 
* `ectyper <https://github.com/phac-nml/ecoli_serotyping>`_ 
* `stxtyper <https://github.com/ncbi/stxtyper>`_
* `ncbi-amrfinderplus <https://github.com/ncbi/amr>`_
* `EzClermont <https://github.com/nickp60/EzClermont>`_
* `Kaptive <https://github.com/klebgenomics/Kaptive>`_  
* `MiST <https://github.com/BioinformaticsPlatformWIV-ISP/MiST>`_
* `ShigaPass <https://github.com/imanyass/ShigaPass>`_


Install Kleborate 
~~~~~~~~~~~~~~~~~~~~~~~~~~

Create a conda environment containing Kleborate dependancies::


conda create -n klebsiella_analysis \
    -c conda-forge \
    -c bioconda \
    --strict-channel-priority \
    python=3.11 \
    minimap2 \
    mash \
    mist_typing \
    ezclermont \
    ectyper \
    ncbi-stxtyper \
    shigapass \
    ncbi-amrfinderplus -y

   

Activate the environment and install kleborate (and Kaptive) using pip::

    conda activate klebsiella_analysis
    pip install kleborate

Or Bioconda::

    conda install -c bioconda kleborate
    pip install rammappy

Database set up
============================
Before running Kleborate, set up the reference databases required

AMRFinderPlus Database
~~~~~~~~~~~~~~~~~~~~~~~~~~

NOTE: AMRFinderPlus is only used by the module ``-m escherichia__amr`` and preset ``-p escherichia``, so if you are not analysing *E. coli/Shigella* you don't need to do this step.

Download the latest AMRFinderPlus database:

.. code-block:: bash

   amrfinder -u


KpSC cgMLST Database (MiST)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

NOTE: KpSC cgMLST is only used by the module ``-m kpsc__cgmlst`` and preset ``-p kpsc``, so if you are not analysing *K. pneumoniae* species complex you don't need to do this step.

The ``kpsc__cgmlst`` module performs core-genome MLST (cgMLST) allele calling using `MiST <https://github.com/BioinformaticsPlatformWIV-ISP/MiST>`_ and uses the ``scgMLST629_S`` scheme hosted on the `Institut Pasteur BIGSdb instance <https://bigsdb.pasteur.fr/klebsiella/>`_.

This requires a local, indexed copy of the cgMLST scheme to be stored in ``kleborate/modules/kpsc__cgmlst/data``. As the scheme is under license, we cannot distribute it with Kleborate, and you will need to download the latest version yourself. To simplify this, we have supplied the Python script `setup_cgmlst.py <https://github.com/klebgenomics/Kleborate/blob/development/kleborate/shared/setup_cgmlst.py>`_ download and index the scheme database:

.. code-block:: bash

   python setup_cgmlst.py

Prerequisites
--------------
Before running the ``setup_cgmlst.py`` script, ensure the following tools are installed and available in your ``PATH``:

* **MiST** (``mist``)
* **bigsdb-downloader** (Required for authenticated Pasteur downloads: ``pip install bigsdb-downloader``)


What the script does
--------------
1. **Verifies Dependencies:** Confirms ``mist`` is accessible in your environment
2. Installs ``bigsdb-downloader``
3. Locates Kleborate data path (``kleborate/modules/kpsc__cgmlst/data``)
4. Downloads the ``scgMLST629_S`` scheme from the Institut Pasteur BIGSdb instance
   
  
Download Modes
--------------
When running ``setup_cgmlst.py``, you will be prompted to select one of two download modes:

1. **Standard download (Public)** Pulls public scheme data (frozen on 2024-12-31) without requiring login credentials.
2. **Latest Pasteur database (Authenticated)**
   Pulls the most up-to-date scheme data directly from Institut Pasteur and requires OAuth authentication. Uses the ``bigsdb_auth`` downloader with credentials stored in ``.bigsdb_tokens/``.

Pasteur Credential Setup (Mode 2)
--------------------------------------
If you select **Mode 2** and valid tokens are not found in ``.bigsdb_tokens/access_tokens``, the script initiates the OAuth setup:

1. **Obtain API Client Credentials:**
   Register for database access via the `Institut Pasteur BIGSdb Portal <https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl>`_.
   * Request an OAuth Client Key and Secret by emailing ``bigsdb@pasteur.fr``.

2. **Run Authentication via the Setup Script and input API Keys:**
   Enter your ``Client ID`` and ``Client Secret`` at the terminal prompts.

3. **Authorize in Browser:**
   Open the generated URL in your browser, log in to your Pasteur account, and copy the verification code.

4. **Complete Verification:**
   Paste the verification code into the terminal. Access tokens will be saved to ``.bigsdb_tokens/``, and subsequent runs will skip re-authentication.

Output Files
------------
Upon completion, the following files will be downloaded inside the target folder:

* ``kleb_scgmlst_s/`` — Raw scheme FASTA alleles and ``profiles.tsv``.
* ``kleb_scgmlst_s-index/`` — Indexed binary database used by ``mist`` during Kleborate runs.



*E. coli* cgMLST Database (MiST)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

NOTE: E. coli cgMLST is only used by the module ``-m escherichia__cgmlst`` and preset ``-p escherichia``, so if you are not analysing *E. coli/Shigella* you don't need to do this step.

The ``escherichia__cgmlst`` module performs core-genome MLST (cgMLST) allele uses the **Escherichia.cgMLSTv1** scheme hosted on `EnteroBase <https://enterobase.warwick.ac.uk/species/index/ecoli>`_.

Run the `setup_ecoli_cgmlst.py <https://github.com/klebgenomics/Kleborate/blob/development/kleborate/shared/setup_ecoli_cgmlst.py>`_ script to download and index the scheme:

.. code-block:: bash

   python setup_ecoli_cgmlst.py

Output Files
------------
The following files will be created inside the ``kleborate/modules/ecoli__cgmlst/data`` folder:

* ``ecoli_cgmlst_v1/`` — Raw scheme FASTA alleles and ``profiles.tsv``.
* ``ecoli_cgmlst_v1-index/`` — Indexed binary database used by ``mist`` during allele calling.


EnteroBase API Token (required for LIN codes)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

EnteroBase schemes do not include LIN codes in the downloaded profiles file. To retrieve the LIN code for a matched cgST, this module queries EnteroBase's API token

Please request API access to the EnteroBase E. coli database at https://enterobase.warwick.ac.uk/
Once, you have obtained the token run; 

.. code-block:: bash

   echo 'YOUR_TOKEN' > ./enterobase_token
   chmod 600 ./enterobase_token

When running the module pass the ``--ecoli_entero_token`` path


.. code-block:: bash

     kleborate \
     -a *.fasta \
     -o kleborate_results \
     -m escherichia__cgmlst \
     --trim_headers \
     --ecoli_entero_token ./enterobase_token


See also
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* `BIGSdb_downloader Documentation <https://github.com/kjolley/BIGSdb_downloader>`_
* `MiST Repository <https://github.com/BioinformaticsPlatformWIV-ISP/MiST/wiki/lincodes>`_
* `Inferring LIN codes with MiST <https://github.com/BioinformaticsPlatformWIV-ISP/MiST/wiki/lincodes>`_
* `Inferring Klebsiella LIN codes with MiST <https://github.com/BioinformaticsPlatformWIV-ISP/MiST/wiki/Klebsiella-LINcodes-case-study>`_



Test installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To test that Kleborate is installed and working correctly, download the example genome assembly and run Kleborate using the  -p kpsc::

   wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/813/595/GCF_002813595.1_ASM281359v1/GCF_002813595.1_ASM281359v1_genomic.fna.gz
   kleborate -a GCF_002813595.1_ASM281359v1_genomic.fna.gz -o kleborate_test -p kpsc

If the installation is successful, the analysis should complete without errors and generate the expected output files.


Additional *K. pneumoniae* test datasets are provided to further validate your installation::

   kleborate -a test/kpsc_test/data/ -o kleborate_kpsc_test -p kpsc

The generated output should match the corresponding reference files located in
``test/kpsc_test/example_output/``:

``test/kpsc_test/example_output/klebsiella_pneumo_complex_output.txt``:

.. csv-table::
   :file: ../../test/kpsc_test/example_output/klebsiella_pneumo_complex_output.txt
   :delim: tab
   :header-rows: 1


``test/kpsc_test/example_output/klebsiella_pneumo_complex_hAMRonization_output.txt``:

.. csv-table::
   :file: ../../test/kpsc_test/example_output/klebsiella_pneumo_complex_hAMRonization_output.txt
   :delim: tab
   :header-rows: 1

``test/kpsc_test/example_output/klebsiella_pneumo_complex_genotype_spec.txt``:

.. csv-table::
   :file: ../../test/kpsc_test/example_output/klebsiella_pneumo_complex_genotype_spec.txt
   :delim: tab
   :header-rows: 1

   
