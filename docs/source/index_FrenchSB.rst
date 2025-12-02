.. Kleborate documentation master file, created by
sphinx-quickstart le jeu avril 25 06:02:56 2024.
Vous pouvez adapter ce fichier complètement à votre goût, mais il devrait au moins
contenir la directive racine "toctree".

.. toctree::
   :maxdepth: 1
   :hidden:

Installation
Utilisation
modules
Création de nouveaux modèles


########################
Présentation de Kleborate v3
########################

Kleborate a été principalement développé pour analyser les assemblages génomiques de *Klebsiella pneumoniae* et du complexe d'espèces *Klebsiella pneumoniae* (KpSC) afin de détecter :

* L'espèce (par exemple : *K. pneumoniae*\ , *K. quasipneumoniae*\ , *K. varicola*\ , etc.)
* Séquençotype MLST pour *K. pneumoniae* 
* Loci de virulence associés à *ICEKp* : yersiniabactine (*ybt*), colibactine (*clb*), salmochéline (*iro*), hypermucosité (*rmpA*)
* Loci associés au plasmide de virulence : salmochéline (\ *iro*\ ), aérobactine (\ *iuc*\ ), hypermucosité (\ *rmpA*\ , *rmpA2*\ )
* Déterminants de résistance aux antimicrobiens : gènes acquis, SNPs, troncations de gènes et β-lactamases intrinsèques
* Prédiction du sérotype K (capsule) et de l'antigène O (LPS), via les allèles *wzi* et `Kaptive <https://github.com/klebgenomics/Kaptive>`_


`Kleborate v3 <https://github.com/klebgenomics/Kleborate>`_ comprend une réécriture du code pour (i) remplacer l'utilisation de BLAST par minimap (plus rapide et plus stable); et (ii) introduire une structure modulaire permettant d'ajouter facilement de nouveaux modules de génotypage, y compris pour d'autres espèces.


Pour le complexe d'espèces *K. pneumoniae*, **Kleborate v3 peut reproduire les sorties de Kleborate v2 en exécutant les modules préréglés pour KpSC via:**
 

.. code-block:: Python

   kleborate -a *.fasta -o kleborate_results -p kpsc --trim_headers

(Notez que la commande a changé depuis Kleborate v2, la commande ci-dessus est équivalente à l'exécution  ``kleborate --all -o results.txt -a *.fasta`` avec Kleborate v2 et comprend toutes les résistances et le typage par Kaptive)

**De nouveaux modules pour d'autres espèces sont en cours d'élaboration**, pour l'instant, il s'agit notamment de schémas MLST pour le complexe des espèces *Klebsiella oxytoca* et pour *Escherichia coli* (voir la page Modules).


Citations
----------

Si vous utilisez Kleborate, veuillez citer la publication : Lam, MMC. et al. A genomic surveillance framework and genotyping tool for *Klebsiella pneumoniae* and its related species complex, *Nature Communications* (2021). `<https://www.nature.com/articles/s41467-021-24448-3>`_

Si vous utilisez Kaptive pour typer des locus K et O, veuillez également citer Wyres, KL. et al. Identification of *Klebsiella* capsule synthesis loci from whole genome data. *Microbial Genomics* (2016). `<http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000102>`_

Les publications ci-dessous fournissent plus d'informations sur les composants de génotypage incorporés dans Kleborate:

..
   
Y.  ersiniabactine et colibactine (*ICEKp*):
   Lam, MMC. et al. Genetic diversity, mobilisation and spread of the yersiniabactin-encoding mobile element *ICEKp* in *Klebsiella pneumoniae* populations. *Microbial Genomics* (2018). `Microbial Genomics <http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000196>`_

   Aérobactine et salmochéline:
   Lam, MMC. et al. Tracking key virulence loci encoding aerobactin and salmochelin siderophore synthesis in *Klebsiella pneumoniae*. *Genome Medicine* (2018). `Genome Medicine <https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-018-0587-5>`_

   Kaptive pour le sérotypage capsulaire (K):
   Wyres, KL. et al. Identification of *Klebsiella* capsule synthesis loci from whole genome data. *Microbial Genomics* (2016). `Microbial Genomics 2 <http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000102>`_

   Kaptive pour le sérotypage de l'antigène O (LPS) :
   Wick, RR et. al. Kaptive Web: user-friendly capsule and lipopolysaccharide serotype prediction for *Klebsiella* genomes. *Journal of Clinical Microbiology* (2018). `Journal of Clinical Microbiology <http://jcm.asm.org/content/56/6/e00197-18>`_


Modifications comparées à la v2
----------------

Lorsque Kleborate v3 est exécuté en utilisant l'option ``-p kpsc`` pour exécuter des modules prédéfinis pour *K. pneumoniae* la même logique est implémentée que Kleborate v2, plus les modifications/mises à jour suivantes:

* Mise à jour des bases de données MLST & virulence (avril 2024)
* La colonne ``Chr_ST`` a été supprimée dans v3, car elle est redondante avec ``ST``
* Base de données AMR mise à jour avec CARD v3.2.19 (juin 2024)
* Ajouté ``$`` pour indiquer quand PmrB ou MgrB ont une mutation dans le codon de départ qui peut perturber la traduction (dans la colonne ``Col_mutations``)
* Vérification supplémentaire de la mutation synonyme dans ompK36 (25 C > T) associée à une résistance accrue aux carbapénèmes (dans la colonne ``Omp_mutations``)
* Mise à jour pour utiliser Kaptive v3, qui a apporté quelques modifications aux noms des variables de sortie:
* ``K_locus_missing_genes`` a été renommé ``K_Missing_expected_genes``
* ``O_locus_missing_genes`` a été renommé ``O_Missing_expected_genes`
* Ajout d'un nouveau rapport de génotypage AMR compatible avec la norme `hAMRonization <https://github.com/pha4ge/hAMRonization/blob/master/schema/PHA4GE%20AMR%20Gene%20%26%20Variant%20Specification.csv>`_. élaborée par Public Health Alliance for Genomic Epidemiology (PHA4GE).
* Ajout d'un module pour la prédiction de la résistance à la Ciprofloxacine
* Les mutations sont signalées en utilisant la nomenclature `HGVS nomenclature <https://github.com/AMRverse/AMRrulesCuration/blob/main/syntax.md>`_.
* Ajout de nouveaux modules pour les espèces *Escherichia*: pathotypage, typage de l'îlot de pathogénicité LEE, ClermonTyping, typage de types stx utilisant StxTyper, sérotypage O:H utilisant ECTyper, typage AMR utilisant AMRFinderPlus



Tutoriel
--------

Un tutoriel étape par étape pour Kleborate v3 est disponible à l'atelier `kleborate-workshop <https://docs.google.com/document/d/1R61bQbBngpiDB2Gl_eXigePBVakYZEjy/edit>`_, couvrant:


* Les caractéristiques de Kleborate et leur justification scientifique
* Comment exécuter Kleborate
* Exemples, illustrant comment exécuter et interpréter les résultats


Rapports publics
----------------

La publication de `Kleborate <https://www.nature.com/articles/s41467-021-24448-3>`_ rapporte les résultats de génotypage de ~10 000 génomes publics filtrés pour éliminer les séquences redondantes (par ex. les clusters épidémiques, identifiés par leur faible distance mash avec la même année d'isolement, le même lieu d'isolement et le même génotype), avec Kleborate v2. Les résultats peuvent être explorés dans `Microreact <https://bit.ly/klebMR>`_ (qui montre l'arbre mash, la sortie de Kleborate et les métadonnées curées) ou `Kleborate-Viz <https://kleborate.erc.monash.edu/>`_ (R shiny app). Kleborate-Viz a également l'ensemble de données EuSCAPE préchargé, et vous pouvez aussi visualiser vos propres résultats Kleborate.

kleborate est également inclus dans `Klebsiella Pathogenwatch <https://pathogen.watch/>`_ qui présente des arbres interactifs, des cartes et des listes d'entrées pour *Klebsiella pneumoniae* et vous permet d'analyser vos propres données dans le contexte des collections publiques. Voir `ce document <https://doi.org/10.1093/cid/ciab784>`_ pour un exemple d'utilisation.

Contactez-nous
----------

Kleborate est développement actif avec de nombreux autres outils et projets d'analyse génomique en cours (voir `github.com/klebgenomics <https://github.com/klebgenomics>`_).

N'hésitez pas à nous contacter via le `tracker de questions de GitHub <https://github.com/klebgenomics/Kleborate/issues>`_ si vous avez des problèmes, des questions ou des idées.

Pour en savoir plus sur notre laboratoire, y compris d'autres logiciels, voir `http://holtlab.net <http://holtlab.net>`_

Licence
-------

`GNU General Public License, version 3 <https://www.gnu.org/licenses/gpl-3.0.html>`_



