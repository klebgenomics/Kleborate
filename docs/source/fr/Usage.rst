
########################
Utilisation
########################

Fichiers d'entrée
-----------

Les assemblages génomiques au format FASTA (peuvent être gzippés).

Peuvent être soit des assemblages incomplets ou complets (complet est mieux parce qu'il réduit le risque de fragmentation des gènes/loci).


Utilisation de base
-----------

Exécuter avec les modules prédéfinis pour le complexe d'espèces *K. pneumoniae* (KpSC), qui reproduit le comportement de Kleborate v2 :

.. code-block:: Python

   kleborate -a *.fasta.gz -o kleborate_results -p kpsc --trim_headers

- ``-a *.fasta.gz``: Spécifie les fichiers d'entrée (assemblages) à analyser (.fasta ou .fasta.gz).
- ``-o``: Spécifie le répertoire où les fichiers de sortie seront enregistrés (un fichier de sortie par espèce/complexe détecté).
- ``-p``: Spécifie les modules prédéfinis à exécuter (kpsc, kosc, escherichia).
- ``--trim_headers``: Elague les noms de module dans les en-têtes de colonnes dans la sortie.


Exécuter seulement avec des modules spécifiés, p.ex. typage AMR pour le complexe d'espèces *K. pneumoniae* :

.. code-block:: Python

   kleborate -a *.fasta -o kleborate_results -m klebsiella_pneumo_complex__amr

(Une liste de modules est disponible via ```kleborate --list_modules``` ou `ici <https://kleboratemodular.readthedocs.io/en/latest/modules.html>`_)


Exécuter avec des modules prédéfinis pour le complexe d'espèces *K. oxytoca*

.. code-block:: Python

   kleborate -a *.fasta -o kleborate_results -p kosc --trim_headers


Exécuter avec des modules préréglés pour *E. coli* ou d'autres *Escherichia*, sur des assemblages gzippés:

.. code-block:: Python

   kleborate  -a *.fasta.gz -o kleborate_results -p escherichia --trim_headers


Vérifier les modules disponibles, vérifier la version, imprimer l'aide:

.. code-block:: Python

   kleborate [--list_modules] [--version] [-h]



Fichiers de sortie
--------------------

Les fichiers de sortie sont des fichiers tabulés (.txt), un par espèce/complexe détecté, nommé dans le format: klebsiella_pneumo_complex_output.txt.
Les colonnes incluses dans chaque fichier de sortie dépendent des modules qui sont exécutés; essentiellement, chaque module crée un ensemble de colonnes de résultats qui sont ajoutées au fichier de sortie pour les espèces/complexes pertinents. Par défaut, chaque nom de colonne est précédé du nom du module qui l'a généré. Cela peut être désactivé en utilisant --trim_headers lors de l'exécution de kleborate, ou ces en-têtes de colonnes peuvent être supprimés plus tard en utilisant le script trim_headers.py.


Paramètres
----------

**Input/output:**

``-a ASSEMBLIES [ASSEMBLIES ...], --assemblies ASSEMBLIES [ASSEMBLIES ...]``
   Fichier(s) FASTA pour les assemblages, éventuellement gzippé (.gz)

``-o OUTDIR, --outdir OUTDIR``
   Répertoire pour stocker les fichiers de sortie (par défaut: Kleborate_results)

``--trim_headers``
   Elague les en-têtes dans les fichiers de sortie (activer pour supprimer les noms de module des en-têtes de colonne dans les fichiers de sortie). Alternativement, les utilisateurs peuvent élaguer les en-têtes plus tard en utilisant ce script: `trim_headers.py <https://github.com/klebgenomics/KleborateModular/blob/main/kleborate/shared/trim_headers.py>`_

**Modules:**

``-p PRESET, --preset PRESET`` 
Préréglages de modules, choisissez parmi :

.. list-table::

   * - kpsc
     - *Klebsiella pneumoniae* species complex

   * - kosc
     - *Klebsiella oxytoca* species complex
                                        
   * - escherichia 
     - *Escherichia* genus

``--list_modules``
   Imprime une liste de tous les modules disponibles puis quitte (par défaut: Faux)

``-m MODULES, --modules MODULES``
   Liste délimitée par des virgules des modules Kleborate à utiliser

Les paramètres spécifiques à chaque module peuvent être trouvés `ici <https://kleboratemodular.readthedocs.io/en/latest/modules.html>`_


**Aide:**
     
``-h, --help`` 
   Affiche un message d'aide et quitte

``--help_all``
   Affiche un message d'aide avec toutes les options du module

``--version``
   Affiche le numéro de version du programme et quitte





