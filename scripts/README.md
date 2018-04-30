## Table of contents

* [Using Kleborate with Microreact](#using-kleborate-with-microreact)
* [Preparing Kleborate's resistance gene database files](#preparing-kleborates-resistance-gene-database-files)

## Using Kleborate with Microreact

If you have a phylogenetic tree of your samples, you can view Kleborate's results on the tree using [Microreact](https://microreact.org/).

#### Requirements

* Kleborate output table
* Phylogenetic tree (Newick or Nexus format)

Importantly, the sample names in the table and tree must match.

#### Usage

```
kleborate_to_microreact.py --kleborate_in Kleborate_results.txt --tree_in tree.nwk --csv_out microreact.csv --tree_out microreact.nwk
```

If that command completes successfully, you can go to [Microreact's upload page](https://microreact.org/upload) and drag the `microreact.csv` and `microreact.nwk` files onto the page. Then explore the data with the available colours, labels and blocks!



## Preparing Kleborate's resistance gene database files

These instructions describe how to generate Kleborate's resistance gene database files: `ARGannot_r2.fasta` and `ARGannot_clustered80_r2.csv`.

This is necessary for a couple reasons: 

1) The [ARG-ANNOT](http://en.mediterranee-infection.com/article.php?laref=283%26titre=arg-annot) labels beta-lactamase genes all together in the 'Bla' category. However, this is too general for _Klebsiella_, so Kleborate uses multiple categories for beta-lactamases, e.g. 'Bla_ESBL' for extended-spectrum beta-lactamases and 'Bla_Carb' for carbapenemases.
2) Some resistance genes don't apply to _Klebsiella_. These are removed from the database to keep Kleborate's output a bit simpler.



#### Requirements

* Python 3 (with BioPython)
* [SRST2](https://github.com/katholt/srst2) (for the files in its data directory)


#### 1. Download beta-lactamase information

The first file you need is from the [Lahey Clinic beta-lactamase classifications](https://www.lahey.org/Studies/). While that site is no longer active, it still contains useful information that is not available elsewhere. Download the data table [here](ftp://ftp.ncbi.nlm.nih.gov/pathogen/betalactamases/Lahey.tab).

The second file you need is the GenBank from [BioProject PRJNA313047](https://www.ncbi.nlm.nih.gov/bioproject/?term=313047). Go to [this list of nucleotide records in that project](https://www.ncbi.nlm.nih.gov/nuccore?term=313047%5BBioProject%5D), click 'Send to', choose 'File', set format to 'GenBank' and then click 'Create File' to start the download.



#### 2. Make the beta-lactamase information table

Run the `bla_info.py` script to create a table containing a description and class for each beta-lactamase allele. It takes as input the two files you just downloaded:
```
./bla_info.py Lahey.tab sequence.gb > bla_info_table
```


#### 3. Make the ARG-ANNOT csv file

Kleborate uses the `ARGannot_r2.fasta` and `ARGannot_clustered80_r2.csv` files. These are included in [SRST2](https://github.com/katholt/srst2), but Kleborate needs slightly modified versions.

Run the `make_argannot_csv.py` to add the beta-lactamase information to SRST2's csv:
```
./make_argannot_csv.py path/to/srst2/data/ARGannot_clustered80_r2.csv bla_info_table > ARGannot_clustered80_r2.csv
```


#### 4. Make the ARG-ANNOT fasta file

Now you need to make a fasta file corresponding to the csv file you just made. This removes alleles that aren't in the csv:
```
./make_argannot_fasta.py path/to/srst2/data/ARGannot_r2.fasta ARGannot_clustered80_r2.csv > ARGannot_r2.fasta
```

That's it! Put the two new files (`ARGannot_r2.fasta` and `ARGannot_clustered80_r2.csv`) in Kleborate's data directory (`Kleborate/kleborate/data/`).
