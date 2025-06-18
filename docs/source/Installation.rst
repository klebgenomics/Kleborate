########################
Installation
########################

Dependencies
=============
Kleborate requires the following software and libraries to be installed and available in your path:


* `Python <https://www.python.org/>`_ v3.9 or later
* `Biopython <https://biopython.org/>`_ v1.75 or later
* `Mash <https://github.com/marbl/Mash>`_ v2.0 or later
* `Minimap2 <https://github.com/lh3/minimap2>`_ 
* `Kaptive <https://github.com/klebgenomics/Kaptive>`_ 
* `DNA Features Viewer <https://edinburgh-genome-foundry.github.io/DnaFeaturesViewer/>`_
* ectyper <https://github.com/phac-nml/ecoli_serotyping>`_ 
* stxtyper <https://github.com/ncbi/stxtyper>`_
* ncbi-amrfinderplus <https://github.com/ncbi/amr>`_
* EzClermont <https://github.com/nickp60/EzClermont>`_


Install Kleborate 
~~~~~~~~~~~~~~~~~~~~~~~~~~

Create a conda environment and install Kleborate dependancies::

   conda create -n klebsiella_analysis -c bioconda python=3.9 minimap2 mash ezclermont ectyper stxtyper ncbi-amrfinderplus -y
   

Activate the conda environment and install kleborate using pip::
   
   conda activate klebsiella_analysis
   
   pip install kleborate

Or Bioconda:

   conda install -c bioconda kleborate


Download the AMRFinder database::

   amrfinder -u


Test installation::

   wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/813/595/GCF_002813595.1_ASM281359v1/GCF_002813595.1_ASM281359v1_genomic.fna.gz
   kleborate -a GCF_002813595.1_ASM281359v1_genomic.fna.gz -o kleborate_test -p kpsc
