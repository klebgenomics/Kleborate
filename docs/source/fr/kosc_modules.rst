****************************************************
Modules pour le complexe d'espèces *Klebsiella oxytoca*
****************************************************

.. _klebsiella_oxytoca_complex__mlst:

.. code-block:: Python

   --preset kosc

Ces modules seront exécutés si le module ``enterobacterales_species``\ confirme l'assemblage génomique d'entrée en tant que membre du complexe d'espèces *K. oxytoca* (KpSC) étiqueté dans l'arbre ci-dessus, et répertorié dans le tableau ci-dessous.

.. list-table::

   * - Klebsiella oxytoca

   * - Klebsiella grimontii

   * - Klebsiella michiganensis

   * - Klebsiella pasteurii

   * - Klebsiella huaxiensis

   * - Klebsiella spallanzanii



MLST du KoSC
---------

.. code-block:: Python

   -m klebsiella_oxytoca_complex__mlst

Les génomes identifiés comme appartenant au complexe d'espèces *K. oxytoca* sont soumis au MLST en utilisant le schéma de 7 gènes décrit pour *K. oxytoca* `\at PubMLST <https://pubmlst.org/org/organismes/klebsiella-oxytoca>`_.

Une copie des définitions des allèles MLST et ST est stockée dans le répertoire /data de ce module.


Paramètres du MLST KoSC
+++++++++++++++++++++++

``--klebsiella_oxytoca_complex__mlst_min_identity`` 

Pourcentage d'identité minimale d'alignement pour klebsiella_oxytoca_complex MLST (par défaut : 90.0)

``--klebsiella_oxytoca_complex__mlst_min_coverage``

Pourcentage minimal de couverture de l'alignement pour klebsiella_oxytoca_complex MLST (par défaut: 80.0)

``--klebsiella_oxytoca_complex__mlst_required_exact_matches``

Au moins ce nombre de correspondances exactes sont nécessaires pour appeler un ST (par défaut: 3)


Sorties du MLST KoSC
++++++++++++++++++++

La sortie du module KoSC MLST consiste en les colonnes suivantes :

.. list-table::

   * - ST
     - sequence type (séquençotype)

   * - gapA, infB, mdh, pgi, phoE, rpoB, tonB
     - numéro d'allèle

* Kleborate tente de signaler le ST correspondant le plus proche si une correspondance précise n'est pas trouvée.
* Les correspondances d'allèles imprécises sont indiquées avec un ``*``.
* Les appels de ST imprécis sont indiqués avec ``-nLV``\ , où n indique le nombre de loci qui ne sont pas en accord avec le ST signalé. Ainsi, ``258-1LV`` indique une variante monolocus (SLV) de ST258, c'est-à-dire 6 sur 7 loci correspondent à ST258.


Typage de virulence du KoSC
---------------------
Les génomes identifiés comme appartenant au complexe d'espèces *K. oxytoca* sont soumis au génotypage de la yersiniabactine, de la colibactine, de l'aérobactine, de la salmochéline et du locus rmp, en utilisant les modules suivants :

.. code-block:: Python

   -m klebsiella__ybst, klebsiella__cbst, klebsiella__abst, klebsiella__smst, klebsiella__rmst

Ces modules ont principalement été conçus pour le typage du complexe d'espèces *K. pneumoniae* et les bases de données sont peuplées de variants détectés dans les génomes du KpSC. Cependant, ils peuvent apparaître dans les génomes du KoSC et ainsi le typage est inclus dans le préréglage pour le KoSC.
