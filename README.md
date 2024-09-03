<p align="center"><picture><source srcset="images/logo-dark.png" media="(prefers-color-scheme: dark)"><img src="images/logo.png" alt="Kleborate logo" width="400"></picture></p>

Kleborate v3 includes a rewrite of [Kleborate](https://github.com/klebgenomics/kleborate). It aims to:
* Modularise the code for easier extensibility and maintenance
* Provide functionality for other species, e.g. _Klebsiella oxytoca_ and _Escherichia coli_
* Replace the BLAST dependency with [minimap2](https://lh3.github.io/minimap2/minimap2.html)

**Available modules for Kleborate v3:**

* Basic contigs stats
* Mash-based species detection for *Klebsiella* and closely-related genera
* MLST on the KpSC aerobactin locus (*iuc* genes)
* MLST on the KpSC colibactin locus (*clb* genes)
* MLST on the KpSC Rmp locus (*rmp* genes)
* MLST on the KpSC salmochelin locus (*iro* genes)
* MLST on the KpSC yersiniabactin locus (*ybt* and *irp* genes)
* Chromosomal MLST for *Escherichia coli* using the Achtman scheme
* Chromosomal MLST for *Escherichia coli* using the Pasteur scheme
* Chromosomal MLST for the *Klebsiella oxytoca* species complex
* Chromosomal MLST for the *Klebsiella pneumoniae* species complex
* Virulence score (0-5) for the *Klebsiella pneumoniae* species complex, based on the results of the abst, cbst and ybst modules
* Genotyping acquired genes and mutations for the *Klebsiella pneumoniae* species complex
* Resistance score (0-3) for the *Klebsiella pneumoniae* species complex
* Resistance gene classes count
* Resistance genes count
* Wzi typing for K antigen prediction
* Module to call [Kaptive](https://github.com/klebgenomics/Kaptive) for K (capsule) and O antigen (LPS) serotype prediction


**For information on how to install, and run Kleborate v3 please visit the [Docs](https://kleboratemodular.readthedocs.io/).**
