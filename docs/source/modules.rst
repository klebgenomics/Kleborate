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


Kleborate v3 includes a range of modules for typing bacterial genomes, most of which are specific to a particular species or complex (*Klebsiella pneumoniae SC*, *Klebsiella oxytoca SC*, *Escherichia coli*). We therefore recommend specifying ``-p`` (*kpsc*, *kosc*, *escherichia*) or ``-m`` (list of modules to run based on the organism). This will run the species detection module first, and if the species matches that specified in --preset, the preset modules for that species will be run (if not, the species will be reported and the remaining fields will be blank). 




**Kleborate modules are divided into:**

1. General Modules
2. Modules for *Klebsiella pneumoniae* species complex
3. Modules for *Klebsiella oxytoca* species complex
4. Modules for *Escherichia* *coli* 


Summary of availabe modules and their output columns
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
   * - :ref:`kpsc__mlst <_kpsc_mlst>`
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
   * - :ref:`kpsc__virulence_score <kpsc__virulence_score>`
     - virulence_score (Score of 0-5)
   * - :ref:`kpsc__amr <kpsc__amr>`
     - AGly_acquired, Col_acquired, Fcyn_acquired, Flq_acquired, Gly_acquired, MLS_acquired, Phe_acquired, Rif_acquired, Sul_acquired, Tet_acquired, Tgc_acquired, Tmt_acquired, Bla_acquired, Bla_ESBL_acquired, Bla_ESBL_inhR_acquired, Bla_Carb_acquired, Bla_chr, SHV_mutations, Omp_mutations, Col_mutations, Flq_mutations, truncated_resistance_hits, spurious_resistance_hits
   * - :ref:`kpsc__resistance_score <Resistance scores and counts>`
     - resistance_score (Score of 0-3)
   * - :ref:`kpsc__resistance_gene_count <Resistance scores and counts>`
     - num_resistance_genes
   * - :ref:`kpsc__resistance_class_count <Resistance scores and counts>`
     - num_resistance_classes
   * - :ref:`kpsc__cipro_prediction <kpsc__cipro_prediction>`
     - Ciprofloxacin_prediction, Ciprofloxacin_profile, Ciprofloxacin_profile_support, Ciprofloxacin_MIC_prediction
   * - :ref:`kpsc__wzi <kpsc__wzi>`
     - wzi allele
   * - :ref:`kpsc__kaptive <kpsc__kaptive>`
     - Best match locus, Best match type, Match confidence, Problems, Identity, Coverage, Length discrepancy, Expected genes in locus, details, Missing expected gene
   * - :ref:`kpsc__cgMLST <kpsc__cgMLST>`
     - cgST,LIN code,Sublineage,Clonal group
   * - :ref:`kpsc__mrk <kpsc__mrk>`
     - ST,mrkA, mrkB, mrkC, mrkD, mrkF, mrkH, mrkI, mrkJ
   * - :ref:`kosc__mlst <kosc__mlst>`
     - ST, gapA, infB, mdh, pgi, phoE, rpoB, tonB
   * - :ref:`kosc__kaptive <kosc__kaptive>`
     - Best match locus, Best match type, Match confidence, Problems, Identity, Coverage, Length discrepancy, Expected genes in locus, details, Missing expected gene
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