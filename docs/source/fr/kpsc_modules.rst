****************************************************
Modules pour complexe d'espèces *Klebsiella pneumoniae*
****************************************************


.. code-block:: Python

   --preset kpsc

Ces modules seront exécutés si le module ``enterobacterales__species``\ confirme l'assemblage génomomique d'entrée en tant que membre du complexe d'espèces *K. pneumoniae* (KpSC) étiqueté dans l'arbre ci-dessous.

Nous avons inclus les numéros de phylogroupes dans le tableau ci-dessous pour la compatibilité avec la littérature plus ancienne, mais ces noms ne sont pas utilisés dans la sortie Kleborate. Voir `cette revue <https://www.nature.com/articles/s41579-019-0315-1>`_ pour un aperçu du complexe d'espèces.


.. figure:: _static/kleborate_species_tree.png
   :align: center
   :width: 90%
   :alt: Klebsiella species tree

.. list-table::
   :header-rows: 1

   * - Species
     - Kp phylogroup\ :sup:`a`
     - Kp phylogroup (alternative)\ :sup:`b`
     - Reference
   * - *K. pneumoniae*
     - Kp1
     - KpI
     - `Brenner, D.J. 1979 Int J Syst Evol Microbiol 29: 38-41 <https://ijs.microbiologyresearch.org/content/journal/ijsem/10.1099/00207713-29-1-38>`_
   * - *K. quasipneumoniae* subsp *quasipneumoniae*
     - Kp2
     - KpIIa
     - `Brisse et al., 2014 Int J Syst Evol Microbiol 64:3146-52 <https://ijs.microbiologyresearch.org/content/journal/ijsem/10.1099/ijs.0.062737-0#tab2>`_
   * - *K. quasipneumoniae* subsp *similipneumoniae*
     - Kp4
     - KpIIb
     - `Brisse et al. 2014 Int J Syst Evol Microbiol 64:3146-52 <https://ijs.microbiologyresearch.org/content/journal/ijsem/10.1099/ijs.0.062737-0#tab2>`_
   * - *K. variicola* subsp *variicola*
     - Kp3
     - KpIII
     - `Rosenblueth et al. 2004 Syst Appl Microbiol 27:27-35 <https://www.sciencedirect.com/science/article/abs/pii/S0723202004702349?via%3Dihub>`_
   * - *K. variicola* subsp *tropica*
     - Kp5
     - `-`
     - `Rodrigues et al., 2019 Res Microbiol ﻿S0923-2508:﻿30019-1 <https://www.sciencedirect.com/science/article/pii/S0923250819300191?via%3Dihub>`_ (described as subsp *tropicalensis* in paper)
   * - *K. quasivariicola*
     - Kp6
     - `-`
     - `Long et al. 2017 Genome Announc 5: ﻿e01057-17 <https://mra.asm.org/content/5/42/e01057-17>`_
   * - *K. africana*
     - Kp7
     - `-`
     - `Rodrigues et al. 2019 Res Microbiol ﻿S0923-2508:﻿30019-1 <https://www.sciencedirect.com/science/article/pii/S0923250819300191?via%3Dihub>`_ (described as *africanensis* in this paper)


:sup:`a` Kp phylogroup numbers as described in `Rodrigues et al. 2019 <https://www.sciencedirect.com/science/article/pii/S0923250819300191?via%3Dihub>`_

:sup:`b` alternative (older) Kp phylogroup numbers as described in `Brisse et al. 2001 <https://ijs.microbiologyresearch.org/content/journal/ijsem/10.1099/00207713-51-3-915#tab2>`_ and `Fevre et al. 2005 <https://aac.asm.org/content/49/12/5149>`_ prior to the identification of *K. variicola* subsp *tropica*\ , *K. quasivariicola* and *K. africana*.

.. _klebsiella_pneumo_complex_mlst:

KpSC MLST
---------

.. code-block:: Python

   -m klebsiella_pneumo_complex__mlst

Les génomes identifiés par Kleborate comme appartenant au complexe d'espèces *K. pneumoniae* sont soumis au MLST à l'aide du schéma 7-gènes décrit dans la `\base de données BIGSdb-Pasteur *K. pneumoniae* hébergée à l'Institut Pasteur <https://bigsdb.pasteur.fr/klebsiella/>`_. Notez que ce schéma n'est pas restreint à *K. pneumoniae sensu stricto* mais couvre l'ensemble du complexe d'espèces.

Une copie des définitions des allèles MLST et des profiles des ST est stockée dans le répertoire **/data** de ce module.

Rhinoscleromatis et Ozaenae
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Le groupe clonal *K. pneumoniae* CG67 est connu sous le nom de *K. pneumoniae* subsp. *rhinoscleromatis* parce qu'il provoque le rhinosclerome (infection granulomateuse chronique du nez et des voies respiratoires supérieures), et le groupe clonal CG91 est connu sous le nom de *K. pneumoniae* subsp. *ozaenae* car il est associé à l'ozène (rhinite atrophique). Pour alerter les utilisateurs, lorsque des ST appartenant à ces groupes clonaux sont détectés par Kleborate, cela est signalé dans la colonne ST, par exemple 'ST67 (subsp. rhinoscleromatis)' ou 'ST97 (subsp. ozaenae)'.

Les ST pertinents sont:

.. list-table::

   * - **Species column**
     - **ST**
     - **MLST column**
   * - *K. pneumoniae*
     - 67, 68, 69, 3772, 3819
     - ST67 (subsp. rhinoscleromatis)
   * - *K. pneumoniae*
     - 90, 91, 92, 93, 95, 96, 97, 381, 777, 3193
       3766, 3768, 3771, 3781, 3782, 3784, 3802, 3803
     - ST91 (subsp. ozaenae)


Paramètres
++++++++++

``--klebsiella_pneumo_complex__mlst_min_identity``

Pourcentage d'identité minimale de l'alignement pour klebsiella_pneumo_complex_MLST (par défaut : 90.0)

``--klebsiella_pneumo_complex__mlst_min_coverage`` 

Pourcentage minimal de couverture d'alignement pour klebsiella_pneumo_complex_MLST (par défaut : 80.0)

``--klebsiella_pneumo_complex__mlst_required_exact_matches`` 

Au moins ce nombre de correspondances exactes sont nécessaires pour appeler un ST (par défaut: 3)


Sorties
+++++++

La sortie du module KpSC MLST comprends la colonne suivante:

.. list-table::

   * - ST
     - sequençotype

   * - gapA, infB, mdh, pgi, phoE, rpoB, tonB
     - numéro d'allele

* Kleborate tente de signaler le ST correspondant le plus proche si une correspondance exacte n'est pas trouvée.
* Les correspondances d'allèles imprécises sont indiquées avec un ``*``.
* Les appels de ST imprécis sont indiqués avec ``-nLV``\ , où n indique le nombre de loci qui ne sont pas en accord avec le ST signalé. Ainsi, ``258-1LV`` indique un variant monolocus (SLV, pour single-locus variant) de ST258, c'est-à-dire 6/7 loci correspondent à ST258.


Modules de virulence du KpSC
----------------------

Des modules de typage sont disponibles pour les cinq locus de virulence acquis qui sont associés à des infections invasives et qui se trouvent à une prévalence élevée chez les souches hypervirulentes *K. pneumoniae* : les sidérophores yersiniabactine (\ *ybt*\ ), aérobactine (\ *iuc*\ ) et salmochéline (\ *iro*\ ), la génotoxine colibactine (\ *clb*\ ) et le locus d'hypermucosité *rmpADC*. Chacun de ces loci comprend plusieurs gènes et ne sera signalé que si plus de 50% des gènes sont détectés.

Il y a aussi un module pour analyser le gène d'hypermucosité alternatif *rmpA2*.

Pour chaque module, si le locus cible est détecté, l'outil de génotypage fera:

* Appeler un séquençotype en utilisant la même logique que pour MLST 7-gènes
* Signaler la lignée phylogénétique associée à chaque ST, comme indiqué ci-dessous et détaillé dans les publications correspondantes.
* Signaler le variant structurel de l'élément génétique mobile habituellement associé à cette lignée phylogénétique (pour *ybt* et *rmpADC* seulement)

Les schémas de type MLST spécifiques aux locus *ybt*\ , *clb*\ , *iuc*\ , *iro* and *rmpADC*, et les allèles *rmpA2*, sont définis dans la base de données *K. pneumoniae* `BIGSdb-Pasteur <https://bigsdb.pasteur.fr/klebsiella/>`_.

Notes sur les rapports sur les allèles de virulence:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Les allèles de virulence sont traités de la même manière que les allèles MLST:

* Pour considérer une correspondance Minimap2, elle doit dépasser à la fois 80% d'identité et 40% de couverture (réglable via les options --min_spurious_identity et --min_spurious_coverage).
* Les correspondances qui ne répondent pas à 90% d'identité et 80% de couverture (réglables via les options ``--min_identity`` et ``--min_coverage``) sont signalées dans la colonne ``spurious_virulence_hits`` mais ne sont pas utilisées pour le séquençotypage.
* Les résultats imparfaits (soit <100 % d'identité ou <100 % de couverture) sont signalés avec un ``*``. Par exemple, ``15*` signifie qu'aucune correspondance parfaite n'a été trouvée, mais la correspondance la plus proche est l'allèle 15.
* Kleborate va ensuite traduire l'allèle en séquence d'acides aminés et chercher des troncations (exprimées en % de longueur d'acides aminés depuis le codon de départ). Si le résultat est inférieur à 90 %, il est ajouté au résultat (par exemple ``15*-42%``\ ).

Notes sur le rapport de séquençotype de virulence :
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Les ST des locus de virulence ne sont signalés que si >50% des gènes dans un locus sont détectés (p. ex. au moins 6 des 11 locus *ybt* sont requis pour déclarer un ST *ybt*).
* Si <50 % des gènes dans un locus sont détectés, Kleborate signale que le ST est ``0`` et la lignée ``-``.
* Si <100% mais >50% des gènes dans un locus sont détectés, Kleborate signalera le locus comme (incomplete), avec le ST correspondant le plus proche et sa lignée phylogénétique correspondante. Par exemple, si seulement 7 des 11 gènes *ybt* sont détectés, cela sera déclaré comme étant ``ybtX; ICEKpX (incomplete)``.
* Pour les génomes à copies multiples d'un locus de virulence (p. ex. une souche qui porte ICE *Kp1* et le plasmide KpVP-1 aura deux copies de *iro* et *rmp*\), Kleborate signalera et attribuera un ST ou un ST correspondant le plus proche à chacun de ces locus de virulence à condition que le locus soit relativement intact dans le génome (c'est-à-dire que plus de 50 % des gènes dans un locus soient présents sur un seul contig) et selon les critères ci-dessus.

.. _klebsiella__ybst:
.. _klebsiella__cbst:

Yersiniabactine et colibacine
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: Python

   -m klebsiella__ybst, klebsiella__cbst

Nous avons déjà exploré la diversité de l'élément conjugué intégratif de *K. pneumoniae* (ICE *Kp*), qui mobilise le locus yersiniabactine *ybt*, en utilisant l'analyse génomique d'un ensemble diversifié de 2498 *Klebsiella* (voir `cet article <http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000196>`_\ ). Dans l'ensemble, nous avons trouvé *ybt* dans environ un tiers des génomes *K. pneumoniae* (et *clb* dans environ 14%). Nous avons identifié 17 lignées distinctes de *ybt* (voir figure) intégrées dans 14 variantes structurales de ICE *Kp* qui peuvent s'intégrer à l'un des quatre sites de tRNA-Asn du chromosome. On a trouvé qu'un type était à diffusion plasmidique. Sur la base de cette analyse, nous avons développé une approche de type MLST pour l'attribution des séquençotypes de yersiniabactine (YbST) et de colibactine (CbST), qui est mise en œuvre dans Kleborate.

Notez que bien que ICE *Kp1* se trouve occasionnellement dans d'autres espèces du KpSC, et même dans d'autres genres d'Enterobacteriaceae (voir `publication originale <http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000196>`_\ ), la plupart des variantes connues incluses dans la base de données proviennent de *K. pneumoniae*.

Les bases de données d'allèles et de YbST ont été mis à jour pour la dernière fois en avril 2024. Le nombre de lignées ybt est maintenant de 28 et le nombre de variantes structurales ICE*Kp* est de 22.


Paramètres ybst
++++++++++++++++++

``--klebsiella__ybst_min_identity``

Pourcentage minimal d'identité d'alignement (par défaut: 90.0)

``--klebsiella__ybst_min_coverage``

Pourcentage minimal d'alignement pour la couverture (par défaut: 80.0)

``--klebsiella__ybst_required_exact_matches``

Au moins ce nombre de correspondances exactes sont nécessaires pour appeler un ST (par défaut: 6)


Sorties Ybst
++++++++++++++++++

La sortie du module ybst consiste en les colonnes suivantes:

.. list-table::

   * - Yersiniabactin
     - Lineage (ICEKp prediction)

   * - YbST
     - Yersiniabactin sequence type

   * - ybtS, ybtX, ybtQ, ybtP, ybtA, irp2, irp1, ybtU, ybtT, ybtE, fyuA
     - allele number (ybt locus)


Paramètres cbst
++++++++++++++++++

``--klebsiella__cbst_min_identity``

Pourcentage minimal d'identité (par défaut : 90,0)

``--klebsiella_cbst_min_coverage``

Pourcentage minimal d'alignement pour la couverture (par défaut : 80,0)

``--klebsiella__cbst_required_exact_matches``

Au moins ce nombre de correspondances exactes sont nécessaires pour appeler un ST (par défaut: 8)


Sorties cbst
++++++++++++++++++

La sortie du module cbst consiste en les colonnes suivantes:

.. list-table::

   * - Colibactin
     - Lineage

   * - CbST
     - Colibactin sequence type

   * - clbA, clbB, clbC, clbD, clbE, clbF, clbG, clbH, clbI, clbL, clbM, clbN, clbO, clbP, clbQ
     - allele number (clb / pks locus)

.. _klebsiella__abst:

.. _klebsiella__smst:

Aérobactine et salmochéline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: Python

   -m klebsiella__abst, klebsiella__smst

Nous avons également étudié la diversité génétique des loci aérobactine (\ *iuc*\ ) et salmochéline (\ *iro*\ ) parmi un ensemble de données de 2733 génomes de *Klebsiella* (voir `cette publication <https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-018-0587-5>`_\ ). Nous avons identifié cinq lignées *iro* et six lignées *iuc*, chacune d'elles étant associée à un emplacement spécifique dans les génomes de *K. pneumoniae* (principalement des plasmides de virulence). Sur la base de cette analyse, nous avons développé une approche de style MLST pour l'attribution des types de séquences d'aérobactine (AbST) et de salmochéline (SmST) qui est implémentée dans Kleborate.

* Les lignées les plus courantes sont *iuc1* et *iro1*\ , qui se trouvent ensemble sur le plasmide de virulence FIBk KpVP-1 (typifié par pK2044 ou pLVPK commun aux clones hypervirulents ST23, ST86, etc.).
* Les lignées *iuc2* et *iro2* ont été associées à l'autre plasmide de virulence FIBk KpVP-2 (typifié par Kp52.145 plasmide II de la souche de laboratoire K2 ST66 connue sous le nom de Kp52.145 ou CIP 52.145 ou B5055).
* *iuc5* et *iro5* proviennent de *E. coli* et sont portés (souvent ensemble) chez *E. coli* sur des plasmides FII qui peuvent être transférés à *K. pneumoniae*.
* Les lignées *iuc2A*\ , *iuc3* et *iro4* ont été associées à d'autres plasmides FIBk nouveaux qui n'avaient pas été décrits précédemment dans *K. pneumoniae*\ , mais dont les séquences sont incluses dans la `publication <https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-018-0587-5>`_.
* Le locus salmochéline présent dans ICE *Kp1* constitue sa propre lignée *iro3*\ , et le locus aérobactine présent dans le chromosome des souches ST67 *K. pneumoniae* subsp *rhinoscleromatis* constitue sa propre lignée *iuc4*.

Note sur la mise à jour de la séquence *iucA* :
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Dans la version 2.2.0 et antérieures de Kleborate, la majorité des allèles *iucA* avaient une longueur de séquence de 1791 pb, à l'exception de ceux associés à la lignée *iuc 5* qui ont une longueur de 1725 pb. A cet égard, *iucA* dans les génomes avec *iuc3* encode un codon stop prématuré entraînant une protéine IucA tronquée significativement et probablement non fonctionnelle (à 2% de la longueur de la séquence intacte en acides aminés), malgré des preuves expérimentales montrant une activité sidérophore chez les isolats *iuc3*\ positifs. À la lumière de cette évidence, les séquences des gènes *iucA* ayant une longueur de ~1791 pb ont été mises à jour à ~1725 pb en supprimant les 66 premiers pb. Ces changements sont capturés dans la version 2.3.0 de Kleborate et règlent le problème de troncation dans les génomes avec *iuc3*\. Les allèles *iucA* et les profils AbST suivants ont également été retirés en raison de la redondance de séquence après la mise à jour:

* allèles: iucA48, iucA49, iucA52
* profils: AbST 70, 82, 83

Les bases de données d'allèles et de ST ont été mis à jour pour la dernière fois en avril 2024.

Paramètres abst
++++++++++++++++++

``--klebsiella__abst_min_identity``

Pourcentage minimal d'identité (par défaut : 90,0)

``--klebsiella_abst_min_coverage``

Pourcentage minimal d'alignement (par défaut : 80,0)

``--klebsiella__abst_required_exact_matches``

Au moins ce nombre de correspondances exactes sont nécessaires pour appeler un ST (par défaut: 3)


Sorties abst
++++++++++++++++++

La sortie du module abst consiste en les colonnes suivantes:

.. list-table::

   * - Aerobactin
     - Lineage (plasmid prediction)

   * - AbST
     - Sequence type

   * - iucA, iucB, iucC, iucD, iutA
     - allele number (iuc locus)

Paramètres  smst
++++++++++++++++++

``--klebsiella__smst_min_identity``

Identité minimale en pourcentage d'alignement (par défaut: 90.0)

``--klebsiella_smst_min_coverage``

Pourcentage minimal d'alignement (par défaut : 80,0)

``--klebsiella__smst_required_exact_matches`` 

Au moins ce nombre de correspondances exactes sont nécessaires pour appeler un ST (par défaut: 2)

Sorties smst
++++++++++++++++++

La sortie du module smst consiste en les colonnes suivantes:

.. list-table::

   * - Salmochelin
     - Lineage (plasmid prediction)

   * - SmST
     - Sequence type

   * - iroB, iroC, iroD, iroN
     - allele number (iro locus)

.. _klebsiella__rmst:

.. _klebsiella__rmpa2:

Locus d'hypermucosité
^^^^^^^^^^^^^^^^^^

.. code-block:: Python

   -m klebsiella__rmst, klebsiella__rmpa2

Le locus *rmpA* est associé au phénotype d'hypermucosité qui est une caractéristique de virulence souvent observée chez les souches hypervirulentes de *K. pneumoniae*. Des travaux récents ont révélé que *rmpA* sert de régulateur de transcription pour les gènes *rmpD* et *rmpC*, et ensemble ces gènes constituent le locus *rmpADC* (ou *rmp*\ ). *rmpC* est impliqué dans la régulation de l'expression des capsules tandis que *rmpD* controle l'hypermucoviscosité (voir l'article sur `rmpC <https://mbio.asm.org/content/10/2/e00089-19>`_ et celui sur `rmpD <https://mbio.asm.org/content/11/5/e01750-20>`_ pour plus d'informations.)

À la lumière de cette information, nous avons examiné et extrait les séquences *rmpA*\ , *rmpD* et *rmpC* des 2733 génomes inclus dans l'étude de l'aérobactine et de la salmochéline, et avons généré un schéma de typage RmST. Nous avons observé quatre lignées *rmp* distinctes, qui étaient associées aux plasmides de virulence KpVP-1 (\ *rmp 1*\ ), KpVP-2 (\ *rmp 2*\ ), *iuc2A* ((\ *rmp 2A*\ ), ICE *Kp1* (rmp 3) et la lignée *rmp4* qui est associée à *K. pneumoniae* CG67 `Lam et al., 2024 BioRxiv <https://www.biorxiv.org/content/10.1101/2024.05.28.596137v1/>`_

Le module klebsiella_rmst crible les génomes pour *rmpADC* et signalera un séquençotype, ainsi que la lignée et l'élément génétique mobile associés.

Le gène *rmpA2* est homologue à *rmpA* et le module klebsiella__rmpa2 crible les génomes pour les allèles de *rmpA2*.

Remarque:
^^^^^^^^

* Les allèles pour chaque gène proviennent de `BIGSdb-pasteur <https://bigsdb.pasteur.fr/klebsiella/>`_\ , tandis que d'autres allèles *rmpA* ont également été ajoutés à Kleborate.
* Les gènes *rmpA* et *rmpA2* partagent ~83 % d'identité en nucléotides.
* Les correspondances Minimap2 nucléotidiques uniques (sans chevauchement) avec une identité >95% et une couverture >50% sont signalés. Notez que plusieurs correspondances du même gène sont signalées si elles sont trouvées. Par exemple, le génome NTUH-K2044 porte *rmpA* dans le plasmide de virulence et aussi dans ICE *Kp1* , qui est rapporté dans la colonne *rmpA* comme étant ``rmpA_11(ICEKp1),rmpA_2(KpVP-1)``.
* Comme pour les autres gènes de virulence, les troncations dans les gènes *rmpA* et *rmpA2* sont exprimées en pourcentage de la longueur des acides aminés depuis le codon de départ, par exemple ``rmpA_5-54%`` indique que la protéine RmpA est tronquée après 54 % de la longueur de la séquence intacte des acides aminés. Ces troncations semblent communes, en raison des insertions et des suppressions dans un poly-G, et entraînent presque certainement une perte de fonction protéique.


Paramètres rmst
++++++++++++++++++

``--klebsiella__rmst_min_identity``

Pourcentage minimal d'identité p(par défaut : 90.0)

``--klebsiella_rmst_min_coverage``

Pourcentage minimal d'alignement (par défaut : 80,0)

``--klebsiella_rmst_required_exact_matches «»

Au moins ce nombre de correspondances exactes sont nécessaires pour appeler un ST (par défaut: 2)


Sorties rmst
++++++++++++++++++

La sortie du module rmst consiste en les colonnes suivantes:

.. list-table::

   * - RmpADC
     - Lineage

   * - RmST
     - Sequence type

   * - rmpA, rmpD, rmpC
     - allele number (rmp locus)



Paramètres rmpA2
++++++++++++++++++

``--klebsiella__rmpa2_min_identity`` 

Pourcentage minimal d'identité pour les allèles rmpA2 (par défaut : 90.0)

``--klebsiella_rmpa2_min_coverage``

Pourcentage minimal de couverture d'alignement pour les allèles rmpA2 (par défaut : 80,0)


Sorties rmpA2
++++++++++++++++++

La sortie du module rmst consiste en les colonnes suivantes:

.. list-table::

   * - rmpA2
     - best matching allele



.. _klebsiella_pneumo_complex__virulence_score:

Score de virulence
^^^^^^^^^^^^^^^^^^

.. code-block:: Python

   -m klebsiella_pneumo_complex__virulence_score

Ce module prend ``klebsiella__abst``, ``klebsiella__cbst``, ``klebsiella__ybst`` comme préalables et calcule un score de virulence, qui va de 0 à 5 comme indiqué ci-dessous. Notez que ni le locus de salmochéline (iro) ni le rmpADC ne sont explicitement considérés dans le score de virulence, par simplicité. Les locus iro et rmpADC apparaissent généralement à côté du locus de l'aérobactine (iuc) sur les plasmides de virulence de Kp, et ainsi la présence de iuc (score de 3-5) implique généralement la présence d'iro et de rmpADC. Toutefois, nous privilégions iuc le calcul de la note, car l'aérobactine est spécifiquement associée à la croissance dans le sang et est un prédicteur plus fort du phénotype d'hypervirulence `voir cette revue <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6349525/>`_. Les loci iro et rmpADC sont également présents occasionnellement avec ybt, dans la variante ICEKp1, mais cela aura tout de même un score de 1.

.. list-table::

* - 0
- négatif pour l' ensemble des loci yersiniabactine (ybt), colibatine (clb), et l'aérobactine (iuc)

* - 1
- yersiniabactine seulement

* - 2
- yersiniabactine et colibactine (ou colibactine seulement)

* - 3
- aérobactine (sans yersiniabactine ou colibactine)

* - 4
- aérobactine avec yersiniabactine (sans colibactine)

* - 5
- yersiniabactine, colibactine et aérobactine



Sorties des scores de virulence
++++++++++++++++++++++++++++++++++++++

Le score de virulence est donné dans la colonne suivante:

.. list-table::

   * - virulence_score
     - Score de 0 à 5, comme décrit ci-dessus




.. _klebsiella_pneumo_complex__amr:

AMR de KpSC
--------

.. code-block:: Python

   -m klebsiella_pneumo_complex__amr

Les gènes de résistance acquis
^^^^^^^^^^^^^^^^^^

Ce module permet d'analyser les génome par rapport à une version curée de la base de données `CARD <https://card.mcmaster.ca/>`_ des allèles de gènes de résistance acquis (voir ci-dessous `spreadsheet <https://figshare.com/articles/dataset/CARD_v3_0_8_AMR_database_curation_for_Kleborate/13256759>`__ pour plus de détails sur la curation), et de les regrouper par classe d'antibiotique. Les gènes chromosomiques *fosA* et *oqxAB* qui sont intrinsèques à tous les KpSC ne sont pas rapportés et ne confèrent généralement pas de résistance à la fosfomycine ou aux fluoroquinolones chez ces espèces.

Kleborate suit la logique de choisir la meilleure correspondance pour allèle, l'annoter avec des informations supplémentaires et de le placer dans une colonne appropriée.

En bref:

* Les correspondance de nucléotides exactes sont préférées, suivies des correspondances d'acides aminés exactes, suivies des correspondances de nucléotides inexactes.
* Les annotations indiquent certains aspects de la correspondance : ``^`` (correspondance inexacte en nucléotides mais correspondance exacte en acides aminés), ``*`` (inexacte dans les deux cas), ``?`` (correspondance incomplète), ``-X%`` (séquence d'acide aminé tronquée), ``$`` (codon de départ muté, la traduction peut être perturbée).
* La colonne indique la confiance en la correspondance : les correspondances fortes vont dans la colonne pour leur classe d'antibiotique, les correspondances tronquées vont dans la colonne ``truncated_resistance_hits`` et les correspondances faibles d'identité/couverture vont dans la colonne ``spurious_resistance_hits``.

Et voici la logique plus en détail:

* Pour considérer une correspondance Minimap, elle doit dépasser à la fois 80% d'identité et 40% de couverture (ajustable via les options ``--min_spurious_identity`` et ``--min_spurious_coverage``).
* Si le résultat est 100% identité et 100% couverture, alors il sera déclaré sans autre annotation (par exemple ``TEM-15``\ ).
* Si aucune correspondance exacte avec les nucléotides n'est trouvée, Kleborate recherche une correspondance exacte avec les acides aminés, et le signalera avec un symbole ``^``. Par exemple, ``TEM-15^`` indique une correspondance exacte avec la séquence protéique TEM-15 mais avec une ou plusieurs différences de nucléotides.
* Si aucune correspondance exacte avec les acides aminés n'est trouvée, la correspondance nucléotidique la plus proche est signalée avec un symbole ``*``. Par exemple ``TEM-15*`` indique qu'il n'y a pas d'appariement précis de nucléotides ou d'acides aminés, mais l'appariement de nucléotides le plus proche est TEM-15.
Si le match est inférieur à 100% de couverture, une ?` est ajouté au résultat par exemple. ``TEM-15? indique une correspondance incomplète à 100% d’identité, et ``TEM-15*?` indique une correspondance incomplète à <100 % d’identité.
* Kleborate va ensuite traduire la correspondance en séquence d'acides aminés et chercher des troncations (exprimées en % de longueur d'acides aminés depuis le codon de départ). Si le résultat est inférieur à 90%, il est ajouté au résultat (par exemple ``TEM-15*-42%``\) et le résultat est indiqué dans la colonne ``truncated_resistance_hits``.
* Si la correspondance est inférieure à 90% d'identité ou à 80% de couverture nucléotidique (réglable par les options ``--min_identity`` et ``--min_coverage``), il est indiqué dans la colonne ``spurious_resistance_hits``. Autrement, il est indiqué dans la colonne pour sa classe (p. ex. ``Bla_ESBL_acquired``\ ).

Notez que Kleborate rapporte des résultats de résistance pour toutes les classes d'antimicrobiens avec des mécanismes de résistance attribuables avec confiance dans KpSC. Toutes ces substances ne sont pas utilisées cliniquement pour le traitement des infections par le KpSC (p. ex., MLS, Rif) mais elles sont toujours signalées car la présence de déterminants de résistance acquis à ces classes intéresse les chercheurs pour d'autres raisons (p. ex., ces gènes peuvent être des marqueurs utiles des MGE et de leur propagation; il est possible que ces agents soient utilisés contre d'autres organismes pour sélectionner le KpSC chez les patients co-infectés ou dans l'environnement). Pour une vue d'ensemble de la résistance aux antimicrobiens et des définitions consensus de la multirésistance, de la résistance étendue (XDR) et de la résistance complète (pan-resistance) chez les Enterobacteriaceae, voir `Magiorakos 2012 <https://www.cliniquemicrobiologieandinfection.com/article/S1198-743X(1461632-3/fulltext>`_\

Les bêta-lactamases SHV
^^^^^^^^^^^^^^^^^^^

Tous les KpSC portent un gène chromosomique de bêta-lactamase (SHV dans *K. pneumoniae*\ , LEN dans *K. varicola*\ , OKP dans *K. quasipneumoniae*\ ) qui confère une résistance cliniquement significative à l'ampicilline. Certains KpSC portent également des allèles SHV mobiles acquis, qui peuvent conférer une résistance supplémentaire aux inhibiteurs et/ou aux bêta-lactamines à spectre étendu.

Kleborate rapportera tous les allèles SHV qu'il détecte et les séparera en colonnes basées sur le phénotype de résistance qu'on prévoit :

* Les allèles SHV associés à la résistance à l'ampicilline seulement, seront rapportés dans la colonne ``Bla_chr`` parce qu'ils sont supposés représenter l'allèle chromosomique. Ces gènes ne sont pas inclus dans le nombre de gènes de résistance acquis.
* D'autres allèles SHV, par exemple ceux qui sont prédits pour coder les ESBL (bêta-lactamases à spectre étendu) ou les bêta-lactamases avec résistance aux inhibiteurs, seront signalés dans les colonnes ``Bla_ESBL_acquired`` ou ``Bla_inhR_acquired`` pertinentes (voir ci-dessous), car ces allèles SHV sont presque toujours portés sur des plasmides. (Cependant, il est possible d'avoir une mutation dans un gène SHV chromosomique qui donne une correspondance avec un allèle ESBL, qui serait également rapporté dans la colonne ``Bla_ESBL_acquired`` et compté comme un gène acquis parce qu'il est très difficile de faire la différence sans exploration manuelle du contexte génétique.)

Les mutations spécifiques et l'attribution des allèles à la classe sont détaillées dans cette publication de KlebNET-GSP: `Tsang et al, 2024 Microbial genomics <https://doi.org/10.1099/mgen.0.001294>`_.


Autres mutations chromosomiques associées à la RAM
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* mutations de résistance aux fluoroquinolones : GyrA 83 & 87 et ParC 80 & 84. Elles apparaissent dans la colonne ``Flq_mutations``.
* Résistance à la colistine due à la troncation ou à la perte des gènes core MgrB ou PmrB. Si ces gènes sont manquants ou tronqués, ces informations seront rapportées dans la colonne «Col_mutations» (les troncations sont exprimées en % de longueur en acides aminés depuis le codon de départ, s'il y a une mutation dans le codon de départ, ceci est indiqué comme «$» pour indiquer que le gène est présent mais ne peut pas être traduit correctement). Notez que si MgrB et PmrB sont présents et non tronqués, rien ne sera signalé à leur propos dans la colonne 'Col'.

* Les troncations OmpK35 et OmpK36 et les mutations ponctuelles réduisent la sensibilité aux bêta-lactamases (`insertions GD ou TD dans la troisième boucle <https://www.nature.com/articles/s41467-019-11756-y>`_ ou `C>T synonyme au nucléotide 25 <https://doi.org/10.1073/pnas.2203593119>`_ ``ompK36_c25t``). Ces informations seront déclarées dans la colonne ``Omp_mutations`` (les troncations sont exprimées en % de la longueur en acides aminés depuis le codon de départ). Notez que si un gène est fragmenté en plusieurs contigs, Kleborate tentera de prédire l'allèle correspondant le plus proche basé sur le fragment le plus long. Si ce fragment le plus long ne contient pas le début du gène, la troncation sera déclarée comme -0%. De plus, si ces gènes core sont présents et non tronqués, rien ne sera signalé dans la colonne 'Omp'. L'effet spécifique des mutations OmpK sur la sensibilité dépend de plusieurs facteurs, dont les combinaisons d'allèles OmpK35 et OmpK36 et les gènes de bêta-lactamase présents (c'est pourquoi nous les signalons dans leur propre colonne séparée des gènes Bla). Voir par exemple `puplication <https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1007218>` et `celui-ci <https://www.nature.com/articles/s41467-019-11756-y>`_ pour plus d'informations sur les gènes OmpK et la résistance.

Notez que ceux-ci ne comptent pas pour le nombre de gènes de résistance acquis, mais comptent pour les classes d'antibiotiques (à l'exception des mutations Omp, dont le spectre des effets dépend de la présence de bêta-lactamases acquises et donc leur impact sur des classes spécifiques de bêta-lactamines est difficile à prédire).


Paramètres AMR
++++++++++++++++++++++++++++++++++++++

``--klebsiella_pneumo_complex__amr_min_identity``

Pourcentage d'identité minimale d'alignement pour klebsiella_pneumo_complex (par défaut: 90.0)

``--klebsiella_pneumo_complex_amr_min_coverage``

Pourcentage minimal de couverture d'alignement pour klebsiella_pneumo_complex (par défaut: 80.0)

``--klebsiella_pneumo_complex__amr_min_spurious_identity``

Identité minimale d'alignement pour klebsiella_pneumo_complex pour résultats fallacieux (par défaut: 80.0)

``--klebsiella_pneumo_complex_amr_min_spurious_coverage``

Pourcentage minimal de couverture d'alignement pour klebsiella_pneumo_complex pour résultats fallacieux (par défaut: 40.0)

Sorties AMR
++++++++++++++++++++++++++++++++++++++

Les résultats du module KpSC AMR sont groupés par classe d'antibiotique (selon `ARG-Annot <https://www.ncbi.nlm.nih.gov/pubmed/24145532>`_ DB), les bêta-lactamases étant ensuite réparties en classes de Lahey (maintenant maintenues à `BLDB <http://www.bldb.eu/>`_\), comme suit:


.. list-table::

   * - AGly_acquired
     - aminoglycosides

   * - Col_acquired
     - colistine
     
   * - Fcyn_acquired
     - fosfomycine
     
   * - Flq_acquired
     - fluoroquinolones
     
   * - Gly_acquired
     - glycopeptides
     
   * - MLS_acquired
     - macrolides
     
   * - Phe_acquired
     - phenicols
     
   * - Rif_acquired
     - rifampine
     
   * - Sul_acquired
     - sulfonamides
     
   * - Tet_acquired
     - tetracyclines
     
   * - Tgc_acquired
     - tigecycline
     
   * - Tmt_acquired
     - trimethoprime
     
   * - Bla_acquired
     - beta-lactamases (autres que SHV) qui n'ont pas d'activité connue de type extended-spectrum, carbapenemase, ou inhibitor-resistance

   * - Bla_ESBL_acquired
     - extended-spectrum beta-lactamases, dont les allèles SHV avec activité BLSE connue
   
   * - Bla_ESBL_inhR_acquired
     - extended spectrum beta-lactamases avec resistance aux inhibiteurs, dont les allèles SHV connus 

   * - Bla_Carb_acquired
     - carbapenemases

   * - Bla_chr
     - allèles SHV associés seulement à la résistance à l'ampicilline (supposés chromosomiques)
   
   * - SHV_mutations
     - mutations dans la beta-lactamase SHV connues pour élargir le spectre
   
   * - Omp_mutations
     - mutations dans OmpK35 et OmpK36 (osmoporines) liées à la résistance
     
   * - Col_mutations
     - reporté si les gènes de MgrB ou PmrB ne sont pas intacts
     
   * - Flq_mutations
     - mutations dans la quinolone-resistance determining region de GyrA ou ParC
     
   * - truncated_resistance_hits
     - liste des gènes acquis dans lesquels la protéine est prédite tronquée (par ex. à cause d'un codon stop ou décalage du cadre de lecture)
     
   * - spurious_resistance_hits
     - liste des gènes acquis pour lesquels l'identité ou la couverture est sous le seuil (defaut <90% identity ou <80% couverture nucléotidique)



De plus, nous fournissons un nouveau rapport de génotypage AMR compatible avec le standard `hAMRonization <https://github.com/pha4ge/hAMRonization/blob/master/schema/PHA4GE%20AMR%20Gene%20%26%20Variant%20Specification.csv>`_. élaboré par `Alliance en santé publique pour l'épidémiologie génomique (PHA4GE) <https://www.biorxiv.org/content/10.1101/2024.03.07.583950v1>`_, améliorant ainsi l'interopérabilité des résultats de Kleborate concernant la résistance aux antibiotiques.


Rapport de hAMronisation pour Kleborate
++++++++++++++++++++++++++++++++++++++


.. list-table::

   * - Input_file_name
     - Le nom du fichier contenant les données de séquence à analyser

   * - Gene_symbol
     - Le nom abrégé d'un gène
     
   * - Mutation
     - Changements de séquence des acides aminés/nucléotides détectés dans la séquence analysée par rapport à une référence
     
   * - Genetic_variation_type
     - La classe de variation génétique détectée
     
   * - Drug_class
     - Ensemble de molécules d'antibiotique
     
   * - Input Sequence ID
     - Un identifiant de séquence(s) moléculaire(s) ou d'entrée d'une base de données de séquences
     
   * - Input_gene_length
     - La longueur (nombre de positions) d'une séquence de gènes cible soumise par un utilisateur
     
   * - Input_gene_start
     - La position du premier nucléotide dans une séquence génique analysée (séquence génique input)
     
   * - Input_gene_stop
     - La position du dernier nucléotide dans une séquence génique analysée (séquence génique input)
     
   * - Reference_gene_length
     - La longueur (nombre de positions) d'une séquence de référence génétique extraite d'une base de données
     
   * - Reference_gene_start
     - Position du premier nucléotide dans une séquence de gènes de référence (séquence utilisée pour la comparaison)
     
   * - Sequence_identity
     - L'identité des séquences est le nombre (%) de correspondances (caractères identiques) dans les positions d'un alignement de deux séquences moléculaires
     
   * - Coverage (percentage)
     - Le pourcentage de la séquence de référence couverte par la séquence d'intérêt

   * - Reference_accession
     - Un identifiant qui spécifie une entrée individuelle dans un entrepôt de séquences public

   * - Strand_orientation
     - L'orientation d'un élément génomique sur la molécule à double brin
     
   * - Software_name
     - Nom d'un outil informatique, application utilisée pour l'analyse des données
     
   * - Software_version
     - La version du logiciel utilisé pour analyser les données
     
   * - Reference_database_name
     - Identifiant d'une base de données biologique ou bioinformatique
     
   * - Reference_database_version
     - La version de la base de données contenant les séquences de référence utilisées pour l'analyse
     
   * - Input_protein_length
     - La longueur (nombre de positions) d'une séquence cible protéique soumise par un utilisateur
     
   * - Reference_protein_length
     - La longueur (nombre de positions) d'une séquence de référence protéique extraite d'une base de données


.. _Resistance scores and counts:

Score de résistance et nombre de résistances
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

L'exécution du module KpSC AMR exécute automatiquement des modules supplémentaires pour générer des comptes de gènes de résistance et de classes d'antibiotiques, et pour calculer un score de résistance. Ces modules prennent ``klebsiella_pneumo_complex_amr`` comme préalable et peuvent être spécifiés manuellement comme suit:

.. code-block:: Python

   -m klebsiella_pneumo_complex__resistance_score, klebsiella_pneumo_complex__resistance_gene_count, klebsiella_pneumo_complex__resistance_class_count


Score de résistance
++++++++++++++++++++++

Ce module calcule un score de résistance, qui va de 0 à 3 comme suit :

.. list-table::

   * - 0
     - pas d' ESBL, pas de carbapénemase (indépendamment de la résistance à la colstine)

   * - 1
     - ESBL, pas de carbapénémase (indépendamment de la résistance à la colistine)

   * - 2
     - Carbapénemase sans résistance à la colistine (indépendamment des gènes ESBL ou des mutations OmpK)

   * - 3
     - Carbapénemase avec résistance à la colistine (indépendamment des gènes ESBL ou des mutations OmpK)



Nombre de gènes de résistance et nombre de classes d'antibiotiques
+++++++++++++++++++++++++++++++++++++++++++++

Ce module mesure combien de gènes de résistance acquis sont présents et combien de classes de médicaments (en *plus* de l'ampicilline à laquelle KpSC est intrinsèquement résistant) ont au moins un déterminant de résistance détecté (c'est-à-dire ignorant les gènes enregistrés dans les colonnes Bla_chr et Bla_acquis).

A noter :

* La présence de *mutations*\ de résistance, et de formes non-ESBL de gènes core SHV/LEN/OKP, ne contribuent pas au nombre de *gènes* de résistance.
* Les mutations contribuent au nombre de classes, p.ex. la résistance à la fluoroquinolone sera comptée si une mutation GyrA est rencontrée, peu importe si un gène de résistance aux quinolones est acquis (\ *qnr*\ ). Les exceptions sont les mutations Omp, qui ne contribuent pas au nombre de classes de médicaments car leur effet dépend du fond de la souche et de la présence d'enzymes bêta-lactamase acquises; par conséquent, ces informations sont fournies dans une colonne séparée, et l'interprétation est laissée à l'utilisateur (voir la page `résistance aux antimicrobiens <https://github.com/katholt/Kleborate/wiki/Antimicrobial-résistance>`_).
* Les gènes rapportés dans les colonnes ``truncated_resistance_genes`` et ``spurious_resistance_genes`` ne contribuent pas aux dénombrements.
* Notez que, puisqu'une classe de médicaments peut avoir plusieurs déterminants de résistance, le nombre de gènes est généralement plus élevé que le nombre de classes.


Sorties des scores de résistance et nombres de résistances
++++++++++++++++++++++++++++++++++++++

Les scores et les nombres de résistance sont affichés dans les colonnes suivantes :

.. list-table::

   * - resistance_score
     - Score de 0-3, comme défini ci-dessus

   * - num_resistance_genes
     - Nombre de résistances acquises

   * - num_resistance_classes
     - Nombre de classes pour lesquelles un déterminant a été acquis (en plus de la résistance intrinsèque à l'ampicilline)


.. _klebsiella_pneumo_complex__cipro_prediction:


Prédiction de la résistance à la ciprofloxacine
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


.. code-block:: Python

   -m klebsiella_pneumo_complex__cipro_prediction

La prédiction de la résistance à la ciprofloxacine est effectuée sur la base de l'attribution du génome à l'un des dix profils génotypiques, basés sur:

i) Le nombre de mutations dans la région de détermination de la résistance aux quinolones (QRDR) de GyrA et de ParC;
(ii) Le nombre de gènes de résistance aux quinolones à médiation plasmidique (PMQR) (c.-à-d. gènes *qep* et *qnr*);
(iii) La présence ou l'absence de *aac(6')-Ib-cr*.
 

Chaque profil de génotype est associé à un phénotype ciprofloxacine, sous la forme d'une affectation catégorique (type sauvage S, type non sauvage I, type non sauvage R) et d'une concentration minimale inhibitrice (CMI).

L'association de chaque profil de génotype avec un phénotype est basée sur l'analyse de ~13000 génomes, par `KlebNET-GSP AMR Genotype-Phenotype Group <https://klebnet.org/amrgenopheno/>`_, et la force de l'évidence de cet ensemble de données est indiquée dans les colonnes **Positive predictive value** et **MIC**. La valeur prédictive positive du profil du génotype est exprimée par le nombre brut de génomes avec ce génotype et le nombre de ceux qui possèdent le phénotype associé. La colonne **MIC** indique la valeur médiane de la CMI et l'intervalle interquartile de toutes les valeurs de la CMI, pour les isolats ayant ce profil de génotype.

Le développement et la validation du classificateur de prédiction de la résistance à la ciprofloxacine sont décrits en détail dans le document suivant `<https://www.biorxiv.org/content/10.1101/2025.09.24.678318v1/>`_.


.. list-table::

   * - **Genotype profile**
     - **Phenotype prediction**
     - **Positive predictive value**
     - **MIC (mg/L), median [interquartile range]**
   * - 0^ QRDR, 0 PMQR, 0 aac(6`)-Ib-cr
     - wildtype S
     - 90.99% S (N=5168/5680)
     - 0.25 mg/L [0.25-0.25]
   * - 0 QRDR, 0 PMQR, 1 aac(6`)-Ib-cr
     - wildtype S
     - 65.22% S (N=105/161)
     - 0.25 mg/L [0.25-0.5]
   * - 0 QRDR, qnrB1, 0 aac(6`)-Ib-cr
     - nonwildtype I
     - 81.25% I/R (n=130/160)
     - 0.5 mg/L [0.5-1]
   * - 1 QRDR, 0 PMQR, 0 aac(6`)-Ib-cr
     - nonwildtype R
     - 77.67% R (N=80/103)
     - 1 mg/L [1-2]
   * - 1 QRDR, 0 PMQR, 1 aac(6`)-Ib-cr
     - nonwildtype R
     - 86.96% R (N=20/23)
     - 2 mg/L [1-2]
   * - >1 QRDR, 0 PMQR, * aac(6`)-Ib-cr
     - nonwildtype R
     - 99.22% R (N=2150/2167)
     - 2 mg/L [2-4]
   * - 0 QRDR, 1^ PMQR, 0 aac(6`)-Ib-cr
     - nonwildtype R
     - 77.47% R (N=423/546)
     - 1 mg/L [1-2]
   * - 0 QRDR, 1 PMQR, 1 aac(6`)-Ib-cr
     - nonwildtype R
     - 94.63% R (N=775/819)
     - 2 mg/L [1-2]
   * - 0 QRDR, >1 PMQR, * aac(6`)-Ib-cr
     - nonwildtype R
     - 97.06% R (N=66/68)
     - 2 mg/L [2-4]
   * - >0 QRDR, >0 PMQR, * aac(6`)-Ib-cr
     - nonwildtype R
     - 99.22% R (N=2421/2440)
     - 4 mg/L [4-4]

* ^ GyrA-87G et GyrA-87H ne sont pas inclus dans le compte QRDR, et qnrB1 est exclu du compte PMQR unique.
* \* indique que le gène peut être présent ou absent
* Notez que *aac(6`)-Ib-cr* est déclaré dans les colonnes AGly_acquired et Flq_acquired.



Les résultats de la prédiction de la résistance à la ciprofloxacine sont rapportés par Kleborate avec quatre colonnes supplémentaires:


.. list-table::

   * - Ciprofloxacin_prediction
     - Indique la prédiction catégorique du phénotype pour ce génome (type sauvage S, type non sauvage I, type non sauvage R)

   * - Ciprofloxacin_profile
     - Indique à lequel des dix profils (du tableau ci-dessus) ce génome a été attribué
     
   * - Ciprofloxacin_profile_support
     - Pourcentage indiquant la valeur prédictive positive du profil du génotype dans **Ciprofloxacin_profile** pour la catégorie S/I/R indiquée dans **Ciprofloxacin_prediction**, d'après les données du groupe Génotype-Phénotype KlebNET-GSP. La fraction entre parenthèses (N=n/x) indique le nombre total de génomes avec ce profil de génotype (n), et le nombre de ceux qui possèdent le phénotype associé (x), qui ont été utilisés pour calculer le pourcentage.
     
   * - Ciprofloxacin_MIC_prediction
     - Indique la distribution MIC observée pour le profil de génotype dans **Ciprofloxacin_profile**, sous forme de valeur médiane et d'intervalle interquartile, d'après les données KlebNET-GSP AMR énumérées dans la colonne **Ciprofloxacin_profile_support**.


.. _klebsiella_pneumo_complex__kaptive:


Typage des locus K et O de KpSC avec Kaptive
-----------------------------------------

.. code-block:: Python

   -m klebsiella_pneumo_complex__kaptive

Ce module exécutera l'outil `Kaptive <https://github.com/klebgenomics/kaptive>`_ v3 pour identifier les loci d'antigènes capsulaire (K) et somatique (O). Voir la documentation `Kaptive <https://kaptive.readthedocs.io/en/latest/>`_ pour plus de détails sur le fonctionnement de Kaptive, les tutoriels et les citations.

Paramètres Kaptive
+++++++++++++++++++

``-t , --threads``

Nombre de threads pour l'alignement (par défaut: 1)

``--k-db, kpsc_k``

Base de données Kaptive pour le typage K-locus

``--o-db, kpsc_o``

Base de données Kaptive pour le typage o-locus



Sorties Kaptive
+++++++++++++++++

Les résultats de Kaptive sont affichés dans les colonnes suivantes :

.. list-table::
   :header-rows: 1

   * - Column Name
     - Description
   * - Best match locus
     - Le type de locus qui correspond le plus à l'ensemblage.
   * - Best match type
     - Le sérotype/phénotype prédit de l'assemblage.
   * - Match confidence
     - Typable ou non.
   * - Problems
     - Caractères indiquant des problèmes avec le locus correspondant (voir problèmes).
   * - Identity
     - Pourcentage pondéré d'identité du locus le mieux assorti à l'assemblage.
   * - Coverage
     - Pourcentage pondéré de couverture du locus le mieux assorti de l'assemblage.
   * - Length discrepancy
     - Si le locus a été trouvé en une seule pièce, c'est la différence entre la longueur du locus et la longueur de l'assemblage.
   * - Expected genes in locus
     - Une fraction indiquant combien de gènes dans le locus le mieux assorti ont été trouvés dans la partie locus de l'assemblage.
   * - Expected genes in locus, details
     - Noms des gènes attendus trouvés dans la partie locus de l'assemblage.
   * - Missing expected genes
     - Une chaîne énumérant les noms de gènes attendus qui n'ont pas été trouvés.



.. _klebsiella_pneumo_complex__wzi:


Typage wzi pour la prédiction de l'antigène K
-----------------------------------------

.. code-block:: Python

   -m klebsiella_pneumo_complex__wzi

Ce module indique la correspondance la plus proche entre les allèles *wzi* dans `BIGSdb <http://bigsdb.pasteur.fr/klebsiella/klebsiella.html>`_. Il s'agit d'un marqueur du type de locus de capsule (KL), qui est prédictif du sérotype capsulaire (K). Bien qu'il n'existe pas de relation 1-1 entre l'allèle *wzi* et le type KL/K, il existe une forte corrélation (voir `Wyres et al., MGen 2016 <http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000102>`_ et `Brisse et al., J Clin Micro 2013 <https://jcm.asm.org/content/51/12/4073.long>`_). Notez que la base de données *wzi* est peuplée d'allèles du complexe des espèces *Klebsiella pneumoniae* et n'est pas fiable pour d'autres espèces.

L'allèle *wzi* peut fournir un moyen pratique de repérer les types de capsules associées à l'hypervirulence (wzi=K1, wzi2=K2, wzi5=K5); ou de repérer le changement de capsule à l'intérieur des clones, par exemple, vous pouvez indiquer quelle lignée ST258 vous avez à partir du type _wzi_ (wzi154: la lignée principale II; wzi29: lignée recombinante I; autres: probablement d'autres lignées recombinantes). Mais les prédictions des locus K du module Kaptive sont plus spécifiques et plus fiables.

Sorties Wzi
+++++++++++++++

Les résultats de typage wzi sont affichés dans les colonnes suivantes :

.. list-table::

   * - wzi
     - L'allèle wzi

   * - K_locus
     - Locus K généralement associé à cet allèle wzi

