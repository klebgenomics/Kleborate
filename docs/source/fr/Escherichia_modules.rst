
****************************************************
Modules pour les espèces *Escherichia*
****************************************************


.. code-block:: Python

   --preset escherichia

Ces modules seront déployés si le module ``enterobacterales__species``\ confirme l'assemblage d'entrée en tant que membre du genre *Escherichia*.

.. _escherichia__mlst_achtman:

.. _escherichia__mlst_pasteur:

*E. coli* MLST
------------

.. code-block:: Python

   -m escherichia__mlst_achtman, escherichia__mlst_pasteur

Les génomes identifiés comme appartenant au genre *Escherichia* sont soumis à l'analyse MLST à l'aide du schéma Achtman 7-locus.

Le schéma MLST Achtman est hébergé sur `EnteroBase <https://enterobase.warwick.ac.uk/>`_.

Nous offrons également une option aux utilisateurs pour exécuter l'analyse MLST en utilisant le schéma Pasteur en exécutant:
``-m escherichia__mlst_pasteur``

Le schéma Pasteur est décrit dans la base de données *Escherichia coli* gérée par `l'Institut du Pasteur <https://bigsdb.pasteur.fr/ecoli/>`_. Pour plus d'informations et de références, voir `BIGSdb <https://bigsdb.pasteur.fr/ecoli/references/>`_.

Les gènes inclus dans chaque schéma sont indiqués dans le tableau de sorties ci-dessous.

Une copie des allèles MLST et des définitions ST utilisées dans chaque module est stockée dans le répertoire ``/data`` du module.


Paramètres MLST *E. coli*
+++++++++++++++++++++++++++

``--escherichia_mlst_achtman_min_identity`` 

Pourcentage minimal d'identité d'alignement pour *Escherichia-Achtman* MLST (par défaut : 90.0)

``--escherichia_mlst_achtman_min_coverage``

Pourcentage minimal de couverture d'alignement pour Escherichia-Achtman MLST (par défaut: 80.0)

``--escherichia_mlst_achtman_required_exact_matches`` 

Nombre minimum d'allèles exacts requis pour appeler un ST (par défaut : 3).

``--escherichia_mlst_pasteur_min_identity`` 

Pourcentage minimal d'identité d'alignement pour Escherichia-Pasteur MLST (par défaut: 90.0)

``--escherichia_mlst_pasteur_min_coverage``

Pourcentage minimal de couverture d'alignement pour Escherichia-Pasteur MLST (par défaut: 80.0)

``--escherichia_mlst_pasteur_required_exact_matches``

Nombre minimum d'allèles exacts requis pour appeler un ST (par défaut : 4).

Sorties MLST *E. coli*
++++++++++++++++++++++

La sortie du module MLST Achtman *E. coli* comprend les colonnes suivantes:

.. list-table::
   :header-rows: 0

   * - ``ST``
     - Sequence type

   * - ``adk``, ``fumC``, ``gyrB``, ``icd``, ``mdh``, ``purA``, ``recA``
     - Numéros d'allèles pour les loci Achtman.


La sortie du module MLST Pasteur *E. coli* comprend les colonnes suivantes:

.. list-table::
   :header-rows: 0

   * - ``ST``
     - Sequence type.

   * - ``dinB``, ``icdA``, ``pabB``, ``polB``, ``putP``, ``trpA``, ``trpB``, ``uidA``
     - Numéros d'allèles pour les loci du schéma Pasteur.


Notes
-----

* Kleborate tente de signaler le ST correspondant le plus proche si une correspondance précise n'est pas trouvée.
* Les correspondances d'allèles imprécises sont indiquées avec un ``*``.
* Les appels de ST imprécis sont indiqués par ``-nLV``\ , où n indique le nombre de loci qui diffèrent du ST signalé. Par exemple, ``258-1LV`` indique un variant monolocus (SLV) de ST258, c'est-à-dire que 6 sur 7 loci correspondent à ST258.


.. _escherichia__pathovar:

*E. coli* Pathotypage
---------------------

.. code-block:: Python

   -m escherichia__pathovar

* Escherichia coli* est globalement divisé en 2 groupes: intestinal diarrheagenic *E. coli* (DEC), and extra-intestinal *E. coli* (ExPEC) `voir article <https://pmc.ncbi.nlm.nih.gov/articles/PMC5156508/>`_. Les DEC comprennent plusieurs pathotypes cliniquement pertinents : entéropathogène *E. coli* (EPEC), entérotoxinogène *E. coli* (ETEC), entérohémorragique *E. coli* (EHEC), productrice de shiga toxine *E. coli* (STEC), entéroaggrégatif *E. coli* (EAEC), entéroinvasif *E. coli* (EIEC) et diffusement adhérent *E. coli* (DAEC) `voir <https://pmc.ncbi.nlm.nih.gov/articles/PMC5114240/>`_ paper. De plus, *Shigella* est considéré comme un pathotype DEC en raison de sa similitude génétique et pathogénique avec les EIEC.

La majorité des pathotypes DEC sont définis par des marqueurs de virulence spécifiques. Toutefois, le rôle pathogène des marqueurs proposés n'est pas bien établi pour les EAEC, DAEC et AIEC.

Marqueurs de virulence des *E. coli* causant des diarrhées
++++++++++++++++++++++++++++++++++++++++++++++

.. list-table:: 
   :header-rows: 1

   * - **Pathotype**
     - **Defining marker**
     - **Virulence determinants**
     - **Location of determinants**
     - **PCR Diagnostic targets**
     - **Other diagnostic targets**

   * - EPEC
     - LEE pathogenicity island
     - LEE pathogenicity island
     - Pathogenicity island
     - ``eae``
     - ``bfpA``

   * - EIEC/*Shigella*
     - pINV
     - pINV
     - Plasmid
     - ``ipaH``
     - Other ``ipa`` genes

   * - ETEC
     - ST or LT
     - ST or LT\nPlus colonisation factors
     - Plasmid; transposon
     - ``elt``, ``est``
     - `-`

   * - EHEC
     - Shiga toxin
     - Stx1 and/or Stx2
     - Prophages
     - ``stx1``, ``stx2``
     - ``eae``, ``ehxA``

   * - EAEC
     - pAA; aggregative adhesion
     - Not known
     - Plasmid
     - ``aggR``, ``aatA``, ``aaiC``
     - `-`

   * - DAEC
     - Afa/ Dr adhesins
     - Not known
     - Not known
     - ``afa/Dr`` adhesins
     - `-`

   * - AIEC
     - Adherent-invasive phenotype
     - Not known
     - Not known
     - none
     - `-`

Comment ça marche
+++++++++++++

Ce module classe les génomes *E. coli* en pathotypes DEC en fonction de la présence ou de l'absence de gènes marqueurs de virulence à l'aide d'une base de données curée `VirulenceFinder <https://cge.food.dtu.dk/services/VirulenceFinder/>`_ DB. Les assemblages d'entrée sont alignés sur la base de données à l'aide de Minimap2, et Kleborate attribue des pathotypes basés sur une logique adaptée de `EnteroBase <https://enterobase.readthedocs.io/en/latest/pipelines/backend-pipeline-phylotypes.html?highlight=pathovar/>`_.

De plus, Kleborate distingue les espèces *Shigella* en fonction du cluster génétique de biosynthèse du O-antigène spécifique au sérotype. Le module aligne les génomes d'entrée sur une séquence de référence curée dérivée du pipeline de sérotypage *Shigella*, `shigatyper <https://github.com/CFSAN-Biostatistics/shigatyper>`_ en utilisant Minimap2.

Toutes les séquences de référence et les définitions de marqueurs utilisées par ce module sont incluses dans le répertoire **/data** de ce module.


Paramètres pathovar *E. coli*
++++++++++++++++++++++++++++++++++

 
``--escherichia__pathovar_min_identity``

Pourcentage minimal d'identité pour le pathotype (par défaut : 90.0).

``--escherichia_pathovar_min_coverage``

Pourcentage minimal de couverture d'alignement pour le pathotype (par défaut : 80,0).


Sorties pathovar *E. coli*
++++++++++++++++++++++++++++

.. list-table:: 
   :header-rows: 0

   * - ``Pathotype``
     - Pathotype prédit

   * - ``Stx1``, ``Stx2``, ``ST``, ``LT``, ``eae``, ``ipaH``
     - Marqueurs de virulence


.. _escherichia__mlst_lee:

Typage de l'îlot de pathogénicité de *E. coli*
----------------------------------------------

.. code-block:: Python

   -m escherichia__mlst_lee

Le locus de l'effacement entérocytaire (LEE) est un îlot de pathogénicité chromosomique de ~40 kb composé de 41 gènes core organisés en cinq opérons `Elliot et al., 1998 <https://onlinelibrary.wiley.com/doi/10.1046/j.1365-2958.1998.00783.x>`_. Il code (i) une protéine adhésive de membrane externe, connue sous le nom de protéine intimine codée par le gène eae (ii) un système de sécrétion de type III (T3SS), et (iii) le récepteur transloqué (Tir) ainsi que les translocons, chaperones, régulateurs et protéines effecteur sécrétées qui sont liés à la virulence.

Kleborate comprend un module pour le sous-typage de l'îlot de pathogénicité LEE. On trouvera des détails sur les sous-types et les lignées de LEE dans le document `Nature Microbiology <https://www.nature.com/articles/nmicrobiol201510>`_.

La base de données de typage LEE est basée sur l'analyse de plus de 250 génomes *E. coli* contenant le LEE et comprend 7 loci (eae (intimin), tir, espA, espB, espD, espH, espZ). Les données sont fournies sous la forme d'une base de données de type MLST, dans laquelle des combinaisons d'allèles sont attribuées à un sous-type LEE, afin de faciliter une nomenclature commune pour les sous-types LEE. Chaque séquence de la base de données représente un groupe d'allèles étroitement liés qui ont été assignés au même type de locus. Le schéma LEE comprend trois lignées distinctes : la lignée 1 comprend les sous-types LEE 1-2; la lignée 2 comprend les sous-types LEE 3-8; la lignée 3 comprend les sous-types LEE 9-30.

Les séquences de référence et les définitions de profil de type MLST sont incluses dans le répertoire **/data** de ce module.


Paramètres
++++++++++

``--escherichia__mlst_LEE_min_identity``

Pourcentage minimal d'identité d'alignement pour ``escherichia_mlst_LEE``. *Défaut :* ``90.0``

``-escherichia_mlst_LEE_min_coverage``

Pourcentage minimal de couverture d'alignement pour ``escherichia_mlst_LEE``. *Défaut:* ``80.0``

``escherichia__mlst_LEE_mlst_required_exact_matches``

Nombre minimum d'allèles exacts requis pour attribuer un ST. *Défaut :* ``3``


Sorties MLST LEE de *E. coli* 
++++++++++++++++++++++++++++

La sortie du module *E. coli* LEE MLST comprend les colonnes suivantes:

.. list-table::

   * - ``LEE_ST``
     - Type de séquence LEE assignée.

   * - ``LEE_lineage``
     - Lignée associée au ST LEE.

   * - ``LEE_eae``, ``LEE_tir``, ``LEE_espA``, ``LEE_espB``, ``LEE_espD``, ``LEE_espH``, ``LEE_espZ``

     - Numéros d'allèles pour chaque LEE.



Notes complémentaires
----------------

* Kleborate tente de signaler le ST correspondant le plus proche si une correspondance exacte n'est pas trouvée.
* Les correspondances d'allèles imprécises sont indiquées avec un ``*``.
* Les appels de ST imprécis sont indiqués avec ``-nLV``\ , où n indique le nombre de loci qui ne sont pas en accord avec le ST signalé. Par exemple, ``ST10-3LV`` indique une variante à trois locus (SLV) de ST10 (c'est-à-dire 4 sur 7 loci correspondent à ST10).


.. _escherichia__stxtyper:


Stxtyper
-----------

.. code-block:: Python

   -m escherichia__stxtyper

Les toxines de Shiga (Stxs) sont des facteurs de virulence clés des * Escherichia coli* STEC. On les trouve également dans *Shigella dysenteriae 1*. Les stx appartiennent à la famille des toxines de type AB et sont divisés en deux groupes antigéniquement distincts : Stx1 et Stx2. Chaque groupe contient plusieurs variants/sous-types: six pour Stx1 (a, b, c, d, e, f) et sept pour Stx2 (a, b, c, d, e, f et g) [`Yano et al., 2023 <https://www.nature.com/articles/s41598-023-3211-8>`_, `Melton-Celsa, 2014 <https://pmc.ncbi.nlm.nih.gov/articles/PMC4270005/>`_]. Ces toxines sont encodées par des bactériophages lysogènes (Stx phage) et les souches STEC peuvent produire soit un seul sous-type Stx, soit une combinaison de sous-types.

Ce module va exécuter StxTyper pour déterminer le type stx. Voir la documentation `StxTyper <https://github.com/ncbi/stxtyper>`_ pour plus de détails sur son fonctionnement.


Sorties StxTyper
+++++++++++++++++++++

Les résultats StxTyper sont affichés dans les colonnes suivantes :

.. list-table::
   :header-rows: 1

   * - Column Name
     - Description
   * - ``Stx_type``
     - La toxine Shiga. Si l'opéron est complet, le sous-type sera signalé (par exemple, ``stx1a``). Si l'opéron est incomplet ou ambigu, une désignation plus large est utilisée: ``stx1``, ``stx2``, ou simplement ``stx`` si l'algorithme ne peut pas décider.
   * - ``operon``
     - Statut de l'opéron détecté. Valeurs possibles:
       ``COMPLETE`` – Opéron complet.  
       ``PARTIAL`` – Operon incomplet.  
       ``PARTIAL_CONTIG_END`` – Opéron partiel probablement tronqué à la limite du contig.  
       ``EXTENDED`` – La séquence de codage s'étend au-delà du codon stop de la référence pour une ou les deux sous-unités.  
       ``INTERNAL_STOP`` – Une sous-unité contient une mutation nonsense.  
       ``FRAMESHIFT`` – Indel détecté dans la séquence codante.  
       ``AMBIGUOUS`` – Bases ambiguës trouvées dans la séquence.  
       ``COMPLETE_NOVEL`` – Opéron complet qui ne peut pas être typé.
   * - ``identity``
     - Pourcentage d'identité pour les sous-unités A et B.
   * - ``target_start``
     - Position de départ de l'alignement.
   * - ``target_stop``
     - Position finale de l'alignement.
   * - ``target_strand``
     - Orientation de la cible.
   * - ``A_reference``
     - Protéine de référence la plus proche pour la sous-unité A.
   * - ``A_identity``
     - Pourcentage d'identité à la référence pour la sous-unité A.
   * - ``A_reference_subtype``
     - Sous-type attribué à la séquence de référence pour la sous-unité A.
   * - ``A_coverage``
     - Pourcentage de la séquence de référence de la sous-unité A couverte par l'alignement.
   * - ``B_reference``
     - Protéine de référence la plus proche pour la sous-unité B.
   * - ``B_reference_subtype``
     - Sous-type attribué à la séquence de référence pour la sous-unité B.
   * - ``B_identity``
     - Pourcentage d'identité à la référence pour la sous-unité B.
   * - ``B_coverage``
     - Pourcentage de la séquence de référence de la sous-unité B couverte par l'alignement.


.. _escherichia__ectyper:

Sérotypage *E. coli* O:H
----------------------

.. code-block:: Python

   -m escherichia__ectyper

Les sérotypes d'E. coli* sont définis par des combinaisons d'antigènes O (lipopolysaccharide) et H (flagellaire). Actuellement, il y a ~183 groupes O et 53 types H qui ont été définis sérologiquement `Ørskov et Ørskov 1984 <https://www.sciencedirect.com/science/article/abs/pii/S0580951708704471/>`_.


O-antigène
++++++++++

L'antigène O est une composante intégrale du Lipopolysaccharide (LPS) présent dans la membrane externe de la bactérie. Le LPS comprend trois composants : le lipide A, un oligosaccharide central (core) et la chaîne polysaccharidique spécifique O (antigène O). Le domaine O-antigène présente une variabilité significative comprenant 10 à 25 unités d'oligosaccharide répétées, chaque unité contenant deux à sept résidus de sucre `Liu et al., 2020 <https://pmc.ncbi.nlm.nih.gov/articles/PMC7685785/>`_. Les gènes responsables de la synthèse des O-antigènes sont généralement présents en tant que groupe de gènes et se situent entre les deux gènes chromosomiques de ménage galF et gnd/ugd `Iguchi et al 2014 <https://pmc.ncbi.nlm.nih.gov/articles/PMC4379981/>`_. Les principales voies impliquées dans l'assemblage, la synthèse et le transport de l'antigène O comprennent, la voie Wzy, la voie dépendante de Wzx/Wzy, encodée par les gènes wzx (O-antigène flippase) et wzy (O-antigène polymérase), et la voie de transport ABC, encodée par wzm et wzt. Ces gènes sont des biomarqueurs idéaux pour prédire les types d'antigènes O.


Antigènes H
++++++++++

Les antigènes H (flagellaire) sont des protéines de surface composées de molécules répétées de la protéine flagelline, qui facilite la motilité bactérienne. Ces antigènes sont numérotés de H1 à H56 (H13, H22 et H50 ne sont pas utilisés) et sont distincts des antigènes O et K. La flagelline est encodée par le gène fliC sur le locus chromosomique ou ses homologues (gènes de la flagelline non-FliC comme FlkA, FllA et FlmA). Des 53 types d'antigènes H bien connus, 44 sont conférés par l'expression du gène FliC, les 9 autres types H sont encodés par des gènes flagellines non fliC. Plus précisément, les codes H3, H35, H36, H47 et H53 sont codés par flkA, H44 et H55 par fllA, H54 par flmA et H17 par flnA.


Kleborate utilise ECTyper pour le sérotypage in-silico. Voir `ECTyper paper <https://pmc.ncbi.nlm.nih.gov/articles/PMC8767331/>`_. pour plus de détails.

Sorties
+++++++

Les sorties du module ECTyper sont les colonnes suivantes:

.. list-table:: 
   :header-rows: 0

   * - ``O-type``
     - Antigène O prédit.

   * - ``H-type``
     - Antigène H prédit.

   * - ``Serotype``
     - Prédiction combinée O/H.

   * - ``QC``
     - Valeurs de contrôle qualité résumant la confiance globale de la prédiction du sérotype.

   * - ``Evidence``
     - Nombre total d'allèles utilisés pour appeler les antigènes O et H.

   * - ``GeneScores``
     - Scores ECTyper pour les antigènes O et H, allant de 0 à 1.

   * - ``AllelesKeys``
     - Meilleure correspondance des allèles de la base de données ECTyper utilisée pour l'attribution du sérotype.

   * - ``GeneIdentities(%)``
     - Pourcentage des valeurs d'identité des allèles de la requête.

   * - ``GeneCoverages(%)``
     - Pourcentage des valeurs de couverture pour les allèles de la requête.

   * - ``GeneLengths``
     - Longueur des gènes (en paires de bases) des allèles de la requête.

   * - ``Warnings``
     - Messages supplémentaires liés à l'état du contrôle qualité ou à d'autres questions touchant la prédiction du sérotype.


.. _escherichia__ezclermont:


ClermonTyping
----------------------

.. code-block:: Python

   -m escherichia__ezclermont


Le genre *Escherichia* comprend plusieurs clades, dont *Escherichia albertii*, *E. fergusonii*, cinq clades *Escherichia* cryptiques (I–V) et *E. coli* sensu stricto. Au sein de *E. coli*, les souches peuvent être divisées en sept groupes principaux : A, B1, B2, C, D, E et F.

Kleborate attribue les génomes à ces phylogroupes et clades en utilisant l'outil `EzClermont <https://pmc.ncbi.nlm.nih.gov/articles/PMC7656184/>`_, qui est basé sur une logique de test PCR in vitro.


Paramètres
++++++++++

``--escherichia__ezclermont_min_length``

Longueur minimale du contig. *Défaut:* ``500``


Sorties
+++++++

.. list-table:: 
   :header-rows: 0

   * - ``Clermont_type``
     - Phylogroupe ou clade assigné.

   * - ``Clermont_profile``
     - Présence ou absence de produits PCR.


.. _escherichia__amr:


Résistance *Escherichia*
------------------------

.. code-block:: Python

   -m escherichia__amr


Ce module analyse les génomes d'entrée pour des gènes de résistance acquis et des mutations ponctuelles connues associées à la résistance en utilisant l'outil `AMRFinderPlus <https://www.nature.com/articles/s41598-021-91456-0/>`_ . Les déterminants identifiés sont regroupés par classe d'antibiotiques.


Paramètres AMR
++++++++++++++++++

``--organism`` 

Utilisé pour détecter les mutations ponctuelles dans les marqueurs de résistance spécifiques aux espèces.

``-t , --threads`` 

Nombre de threads à utiliser pour l'alignement.


Sorties AMR
++++++++++++++++++

Les résultats du module de résistance *Escherichia* sont regroupés par classe d'antibiotiques :

.. list-table::
   :header-rows: 0

   * - ``Aminoglycoside``
     - Les gènes de résistance aux aminoglycosides.

   * - ``Fluoroquinolone``
     - Les gènes de résistance aux fluoroquinolones.

   * - ``Fosfomycin``
     - Les gènes de résistance à la fosfomycine.

   * - ``Sulfonamide``
     - Les gènes de résistance au sulfamide.

   * - ``Tetracycline``
     - Les gènes de résistance à la tétracycline.

   * - ``Glycopeptide``
     - Les gènes de résistance aux glycopeptides.

   * - ``Colistin``
     - Les gènes de résistance à la colistine.

   * - ``Phénicol``
     - Les gènes de résistance aux phénicols.

   * - ``Macrolide``
     - Les gènes de résistance aux macrolides.

   * - ``Rifamycin``
     - Les gènes de résistance à la rifampine.

   * - ``Trimethoprim``
     - Les gènes de résistance au triméthoprime.

   * - ``BetaLactam``
     - Les gènes de bêta-lactamase.

   * - ``Carbapenem``
     - Les gènes de carbapénemases.

   * - ``Cephalosporin``
     - Les gènes de résistance aux céphalosporines de troisième génération.

   * - ``Methicillin``
     - Les gènes de résistance à la méthicilline.

   * - ``Other Classes``
     - Les gènes de résistance à d'autres catégories d'agents antimicrobiens.
