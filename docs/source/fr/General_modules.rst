
########################
Modules généraux
########################

.. _species_detection:

Détection des espèces
-----------------

``-m enterobacterales__species``


Ce module tentera d'identifier l'espèce de chaque assemblage d'entrée. Il le fait en comparant l'assemblage à l'aide de `Mash <https://mash.readthedocs.io/>`_ à un ensemble de *Klebsiella* et d'autres *Enterobacteriaceae* de NCBI, et en déclarant l'espèce de la correspondance la plus proche.

Paramètres
++++++++++++++++++

``--enterobacterales__species_strong``

Seuil de distance de Mash pour une correspondance d'espèces fortes (par défaut : 0,02)

``--enterobactéries_species_faible``

Seuil de distance de Mash pour une correspondance d'espèce faible (par défaut : 0,04)


Sorties
+++++++

La sortie du module de détection des espèces consiste en les colonnes suivantes :

.. list-table::

   * - species
     - Nom de l'espèce (nom scientifique)

   * - species_match
     - La confiance en l'espèce, indiquée comme ``strong``\ (distance Mash ≤ 0,02) ou ``weak``\ (Distance Mash >0,02 et ≤0,04, peut être une espèce nouvelle ou hybride)

La qualité et l'exhaustivité des résultats de Kleborate dépendent de la qualité des assemblages génomiques d'entrée. En général, vous pouvez vous attendre à de bons résultats avec des génomes assemblés avec des outils comme SPAdes de haute profondeur (>50x). Cependant, il est toujours possible que les gènes clés soumis au génotypage puissent être interrompus entre 2 contigs, ce qui peut créer des problèmes pour les détecter et les typer avec précision.


Statistiques des contigs
------------

.. _contig_stats:

.. code-block:: Python

   -m general__contig_stats

Ce module prend ``enterobacterales_species`` comme un préalable et génère des statistiques de base pour aider les utilisateurs à comprendre leurs résultats de typage dans le contexte de la qualité de l`assemblage, bien que nous recommandons aux utilisateurs d`effectuer eux-mêmes un QC plus complet avant de génotyper des génomes (par exemple, un criblage de contamination, etc.).

Le module retourne un ensemble standard de paramètres de qualité d'assemblage (voir les sorties ci-dessous).


Il indiquera également dans la colonne ``QC_warnings``\ si une taille d'assemblage tombe en dehors de celles spécifiées dans ``species_specification.txt``\ dans le répertoire du module, ou si N50 <10 kbp ou si des bases ambiguës (Ns) sont détectées dans la séquence.

Sorties
+++++++

La sortie du module de statistiques des contigs consiste en les colonnes suivantes :


.. list-table::

   * - contig_count
     - Nombre de contigs dans l'assemblage d'entrée

   * - N50
     - `N50 <https://en.wikipedia.org/wiki/N50,_L50,_and_related_statistics>`_ calculé à partir des tailles de contig

   * - largest_contig
     - Taille du plus grand contig (en bp)

   * - total_size
     - Taille totale de l'assemblage (en bp)

   * - ambiguous_bases
     - Détection de bases ambiguës (yes or no). Si oui, le nombre de bases ambiguës est également indiqué entre parenthèses.

   * - QC_warnings
     - Liste des problèmes de QC détectés, y compris : ``ambiguous_bases``\ (bases ambiguës détectées) ``N50``\ (N50 < 10 kbp), ``total_size`` (la taille des génomes dépasse la plage prévue).
