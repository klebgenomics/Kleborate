########################
Installation
########################

Dépendances
=============
Kleborate exige que les logiciels et bibliothèques suivants soient installés et disponibles dans votre chemin :


* `Python <https://www.python.org/>`_ v3.9 or later
* `Biopython <https://biopython.org/>`_ v1.75 or later
* `Mash <https://github.com/marbl/Mash>`_ v2.0 or later
* `Minimap2 <https://github.com/lh3/minimap2>`_ 
* `Kaptive <https://github.com/klebgenomics/Kaptive>`_ 
* `DNA Features Viewer <https://edinburgh-genome-foundry.github.io/DnaFeaturesViewer/>`_
* `ectyper <https://github.com/phac-nml/ecoli_serotyping>`_ 
* `stxtyper <https://github.com/ncbi/stxtyper>`_
* `ncbi-amrfinderplus <https://github.com/ncbi/amr>`_
* `EzClermont <https://github.com/nickp60/EzClermont>`_


Installer Kleborate
~~~~~~~~~~~~~~~~~~~~~~~~~~

Créer un environnement conda et installer les dépendances de Kleborate::

   conda create -n klebsiella_analysis -c bioconda python=3.9 minimap2 mash ezclermont ectyper stxtyper ncbi-amrfinderplus -y
   

Activer l'environnement conda et installer kleborate en utilisant pip::
   
   conda activate klebsiella_analysis
   
   pip install kleborate

Ou Bioconda:

   conda install -c bioconda kleborate


Téléchargez la base de données AMRFinder::

   amrfinder -u


Tester l'installation::

   wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/813/595/GCF_002813595.1_ASM281359v1/GCF_002813595.1_ASM281359v1_genomic.fna.gz
   kleborate -a GCF_002813595.1_ASM281359v1_genomic.fna.gz -o kleborate_test -p kpsc
