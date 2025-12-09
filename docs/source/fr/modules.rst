.. _modules:

**************************************
Modules
**************************************

.. toctree::
   :maxdepth: 1
   :hidden:

   General_modules
   kpsc_modules
   kosc_modules
   Escherichia_modules


Kleborate v3 inclut différents modules pour le typage de genomes bactériens, la plupart desquels sont spécifiques d'une espèce particulière des complexes d'espèces ou espèces (*Klebsiella pneumoniae Species Complex, KpSC*: *Klebsiella oxytoca SC*, *Escherichia coli*). Nous recommandons donc de spécifier ``-p`` (*kpsc*, *kosc*, *escherichia*) ou ``-m`` (liste des modules à exécuter en fonction de l'organisme). Le module de détection des espèces sera exécuté en premier, et si l'espèce correspond à celle spécifiée dans --preset, les modules prédéfinis pour cette espèce seront exécutés (si ce n'est pas le cas, l'espèce sera signalée et les champs restants seront vides).




**Les modules de kleborate sont divisés en :**

1. Modules généraux
2. Modules pour le complexe d'espèces *Klebsiella pneumoniae*
3. Modules pour le complexe d'espèces *Klebsiella oxytoca*
4. Modules pour *Escherichia* *coli*


Résumé des modules disponibles et de leurs colonnes de sortie
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table:: 
   :header-rows: 1
   :widths: 30 70
   :class: colwidths-given

   * - **Module Name**
     - **Columns**
   * - :ref:`enterobacterales__species <species_detection>`
     - species, species_match
   * - :ref:`general__contig_stats <contig_stats>`
     - contig_count, N50, largest_contig, total_size, ambiguous_bases, QC_warnings
   * - :ref:`klebsiella_pneumo_complex__mlst <klebsiella_pneumo_complex_mlst>`
     - ST, gapA, infB, mdh, pgi, phoE, rpoB, tonB
   * - :ref:`klebsiella__ybst <klebsiella__ybst>`
     - YbST, Yersiniabactin, ybtS, ybtX, ybtQ, ybtP, ybtA, irp2, irp1, ybtU, ybtT, ybtE, fyuA
   * - :ref:`klebsiella__cbst <klebsiella__cbst>`
     - CbST, Colibactin, clbA, clbB, clbC, clbD, clbE, clbF, clbG, clbH, clbI, clbL, clbM, clbN, clbO, clbP, clbQ
   * - :ref:`klebsiella__abst <klebsiella__abst>`
     - AbST, Aerobactin, iucA, iucB, iucC, iucD, iutA
   * - :ref:`klebsiella__smst <klebsiella__smst>`
     - Salmochelin, SmST, iroB, iroC, iroD, iroN
   * - :ref:`klebsiella__rmst <klebsiella__rmst>`
     - RmST, RmpADC, rmpA, rmpD, rmpC
   * - :ref:`klebsiella__rmpa2 <klebsiella__rmpa2>`
     - rmpA2
   * - :ref:`klebsiella_pneumo_complex__virulence_score <klebsiella_pneumo_complex__virulence_score>`
     - virulence_score (Score of 0-5)
   * - :ref:`klebsiella_pneumo_complex__amr <klebsiella_pneumo_complex__amr>`
     - AGly_acquired, Col_acquired, Fcyn_acquired, Flq_acquired, Gly_acquired, MLS_acquired, Phe_acquired, Rif_acquired, Sul_acquired, Tet_acquired, Tgc_acquired, Tmt_acquired, Bla_acquired, Bla_ESBL_acquired, Bla_ESBL_inhR_acquired, Bla_Carb_acquired, Bla_chr, SHV_mutations, Omp_mutations, Col_mutations, Flq_mutations, truncated_resistance_hits, spurious_resistance_hits
   * - :ref:`klebsiella_pneumo_complex__resistance_score <Resistance scores and counts>`
     - resistance_score (Score of 0-3)
   * - :ref:`klebsiella_pneumo_complex__resistance_gene_count <Resistance scores and counts>`
     - num_resistance_genes
   * - :ref:`klebsiella_pneumo_complex__resistance_class_count <Resistance scores and counts>`
     - num_resistance_classes
   * - :ref:`klebsiella_pneumo_complex__cipro_prediction <klebsiella_pneumo_complex__cipro_prediction>`
     - Ciprofloxacin_prediction, Ciprofloxacin_profile, Ciprofloxacin_profile_support, Ciprofloxacin_MIC_prediction
   * - :ref:`klebsiella_pneumo_complex__wzi <klebsiella_pneumo_complex__wzi>`
     - wzi allele
   * - :ref:`klebsiella_pneumo_complex__kaptive <klebsiella_pneumo_complex__kaptive>`
     - Best match locus, Best match type, Match confidence, Problems, Identity, Coverage, Length discrepancy, Expected genes in locus, details, Missing expected gene
   * - :ref:`klebsiella_oxytoca_complex__mlst <klebsiella_oxytoca_complex__mlst>`
     - ST, gapA, infB, mdh, pgi, phoE, rpoB, tonB
   * - :ref:`escherichia__mlst_pasteur <escherichia__mlst_achtman>`
     - ST, dinB, icdA, pabB, polB, putP, trpA, trpB, uidA
   * - :ref:`escherichia__mlst_achtman <escherichia__mlst_achtman>`
     - ST, adk, fumC, gyrB, icd, mdh, purA, recA
   * - :ref:`escherichia__pathovar <escherichia__pathovar>`
     - Pathotype, Stx1, Stx2, ST, LT, eae, ipaH
   * - :ref:`escherichia__mlst_lee <escherichia__mlst_lee>`
     - LEE_ST, LEE_lineage, LEE_eae, LEE_tir, LEE_espA, LEE_espB, LEE_espD, LEE_espH, LEE_espZ
   * - :ref:`escherichia__stxtyper <escherichia__stxtyper>`
     - Stx_type, operon, identity, target_start, target_stop, target_strand, A_reference 'A_identity, A_reference_subtype, A_coverage,B_reference, B_reference_subtype, B_identity, B_coverage
   * - :ref:`escherichia__ectyper <escherichia__ectyper>`
     - O-type, H-type, Serotype, QC, Evidence, GeneScores, AllelesKeys, GeneIdentities(%), GeneCoverages(%), GeneLengths, Warnings
   * - :ref:`escherichia__ezclermont <escherichia__ezclermont>`
     - Clermont_type, Clermont_profile
   * - :ref:`escherichia__amr <escherichia__amr>`
     - Aminoglycoside, Fluoroquinolone, Fosfomycin, Sulfonamide, Tetracycline, Glycopeptide, Colistin, Phenicol, Macrolide, Rifamycin, Trimethoprim, BetaLactam, Carbapenem, Cephalosporin,    Methicillin, Other Classes











