# Kleborate with PlasmidFinder Integration
> [!IMPORTANT]
> This document provides a comprehensive guide for installing and using the PlasmidFinder module within **Kleborate 3**. Please follow the instructions carefully to ensure proper setup and functionality.
## Overview

This module provides PlasmidFinder integration for Kleborate 3, enabling plasmid replicon detection in bacterial genomes. 
> [!CAUTION]
> Note that the `general__plasmidfinder` module is not yet officially integrated into Kleborate and requires manual execution through the provided runner script.

## Prerequisites

This module requires:
- Kleborate 3
- PlasmidFinder
- KMA (K-mer Alignment)
- Python 3.9+

## Installation

### Step 1: Create Kleborate environment

Create a new conda environment with all required dependencies:

```bash
conda create -n klebsiella_analysis -c bioconda python=3.9 minimap2 mash \
  ezclermont ectyper stxtyper ncbi-amrfinderplus plasmidfinder kma -y
```

**Note:** If you already have a Kleborate environment, you can add the missing packages:
```bash
conda install -c bioconda plasmidfinder kma
```

### Step 2: Download AMRFinder Database

Initialize the AMRFinder database:

```bash
amrfinder -u
```

### Step 3: Install Kleborate

Install Kleborate from bioconda:

```bash
conda install -c bioconda kleborate
```

### Step 4: Clone Kleborate repository

Clone the Kleborate repository to access the PlasmidFinder module:

```bash
git clone https://github.com/tnmquann/Kleborate.git
```

## Usage

### Running with PlasmidFinder module

Since the `general__plasmidfinder` module is not yet officially integrated, you must use the custom runner script:

```bash
python /path/to/Kleborate/kleborate-runner.py [OPTIONS]
```

#### Setting up environment variables (Optional)

> [!TIP]
> To simplify command invocation, create an alias or environment variable for the runner script:

```bash
# Add to your shell configuration (~/.bashrc, ~/.zshrc, etc.)
export KLEBORATE_RUNNER="/path/to/Kleborate/kleborate-runner.py"

# Usage:
python $KLEBORATE_RUNNER [OPTIONS]
```

## PlasmidFinder module parameters

### View available parameters

Display all module parameters:

```bash
python /path/to/Kleborate/kleborate-runner.py --help_all
```

This will output the PlasmidFinder module options:

```
general__plasmidfinder module:
  --plasmidfinder_dbpath PLASMIDFINDER_DBPATH
                        Path to custom PlasmidFinder database directory
  --plasmidfinder_dbselection PLASMIDFINDER_DBSELECTION
                        Comma-separated list of databases to use (default: all)
  --plasmidfinder_mincov PLASMIDFINDER_MINCOV
                        Minimum coverage for plasmid replicon match (default: 0.60)
  --plasmidfinder_threshold PLASMIDFINDER_THRESHOLD
                        Minimum identity threshold for plasmid replicon match (default: 0.90)
  --plasmidfinder_extendoutput
                        Generate extended output with alignment files and TSV results
                        (output to {outdir}/{strain}/)
```

### Parameter mappings

The following table maps PlasmidFinder module parameters to native PlasmidFinder options:

| Kleborate Parameter | PlasmidFinder Parameter | Description |
|---|---|---|
| `--plasmidfinder_dbpath` | `-p/--databasePath` | Path to custom PlasmidFinder database |
| `--plasmidfinder_dbselection` | `-d/--databases` | Comma-separated database selection |
| `--plasmidfinder_mincov` | `-l/--mincov` | Minimum coverage threshold (0.0-1.0) |
| `--plasmidfinder_threshold` | `-t/--threshold` | Minimum identity threshold (0.0-1.0) |
| `--plasmidfinder_extendoutput` | `-x/--extended_output` | Enable extended output files |

### Parameter details

**`--plasmidfinder_dbpath`**  
Corresponds to PlasmidFinder's `-p/--databasePath` option. Specifies the path to a custom PlasmidFinder database directory. If not provided, the default database from the PlasmidFinder installation will be used.

**`--plasmidfinder_dbselection`**  
Corresponds to PlasmidFinder's `-d/--databases` option. Allows selective analysis using specific databases. Provide a comma-separated list of database names (e.g., `enterobacteriales,gram_positive`). If not specified, all available databases will be used.

**`--plasmidfinder_mincov`**  
Corresponds to PlasmidFinder's `-l/--mincov` option. Sets the minimum coverage threshold for plasmid replicon matches. Values must be between 0 and 1 (default: 0.60). A lower value includes more marginal hits; a higher value increases stringency.

**`--plasmidfinder_threshold`**  
Corresponds to PlasmidFinder's `-t/--threshold` option. Sets the minimum sequence identity threshold for plasmid replicon matches. Values must be between 0 and 1 (default: 0.90). Higher values require closer sequence matches to the reference database.

**`--plasmidfinder_extendoutput`**  
Corresponds to PlasmidFinder's `-x/--extended_output` option. When enabled, generates extended output files including sequence alignments, TSV-formatted results, and detailed reports. Output files are placed in `{outdir}/{strain}/` subdirectories.

For more information, refer to the [PlasmidFinder documentation](https://bitbucket.org/genomicepidemiology/plasmidfinder/src/master/).

## Example workflows

### Basic Analysis with Default parameters

Run PlasmidFinder with default thresholds (60% coverage, 90% identity) and KpSC preset:

```bash
python /path/to/Kleborate/kleborate-runner.py \
  -a /path/to/sample.fasta \
  -o /path/to/output_dir \
  -m general__plasmidfinder \
  -p kpsc
```

### Analysis with extended output

Generate detailed output including alignment files and TSV results:

```bash
python /path/to/Kleborate/kleborate-runner.py \
  -a /path/to/sample.fasta \
  -o /path/to/output_dir \
  -m general__plasmidfinder \
  --plasmidfinder_extendoutput \
  -p kpsc
```

**Output structure:**
```
/path/to/output_dir/
├── general__plasmidfinder_output.txt    (Main results)
└── {strain_name}/                       (Extended output directory)
    ├── plasmidfinder_results_tab.tsv
    ├── plasmidfinder_Hit_in_genome_seq.fsa
    ├── plasmidfinder_Plasmid_seqs.fsa
    ├── plasmidfinder_results.txt
    └── plasmidfinder_data.json
```

### Batch processing

Process multiple FASTA files in a directory:

```bash
python /path/to/Kleborate/kleborate-runner.py \
  -a /path/to/*.fasta \
  -o /path/to/batch_output \
  -m general__plasmidfinder \
  --plasmidfinder_extendoutput \
  -p kpsc
```

### Custom database and stringent thresholds

Use a custom database with stricter matching criteria:

```bash
python /path/to/Kleborate/kleborate-runner.py \
  -a /path/to/sample.fasta \
  -o /path/to/output_dir \
  -m general__plasmidfinder \
  --plasmidfinder_dbpath /path/to/custom/db \
  --plasmidfinder_dbselection enterobacteriales \
  --plasmidfinder_mincov 0.80 \
  --plasmidfinder_threshold 0.95 \
  -p kpsc
```

### Specific database selection

Analyze only specific plasmid databases:

```bash
python /path/to/Kleborate/kleborate-runner.py \
  -a /path/to/sample.fasta \
  -o /path/to/output_dir \
  -m general__plasmidfinder \
  --plasmidfinder_dbselection "enterobacteriales,gram_positive" \
  -p kpsc
```

## Output interpretation

### Main output file

The primary output file `general__plasmidfinder_output.txt` contains:

| Column | Description |
|---|---|
| `plasmidlist` | Identified plasmid replicons (semicolon-separated if multiple) |
| `plasmiddb` | Database category for each plasmid |
| `plasmidident` | Sequence identity percentage for each hit |
| `plasmidcov` | Coverage percentage for each hit |
| `plasmidcontig` | Contig name containing the replicon |
| `plasmidcontigposition` | Position coordinates within the contig |
| `plasmidnote` | Additional classification notes (e.g., VIR for virulence factors) |
| `plasmidaccession` | GenBank accession number of the matched replicon |

### Extended output files

When using `--plasmidfinder_extendoutput`, additional files are generated:

- **results_tab.tsv**: Tabular results with detailed hit information
- **Hit_in_genome_seq.fsa**: FASTA sequences of detected replicon hits
- **Plasmid_seqs.fsa**: Reference plasmid sequences from the database
- **results.txt**: Human-readable summary report
- **data.json**: Complete JSON output from PlasmidFinder

## Troubleshooting

### ModuleNotFoundError

**Error:** `ModuleNotFoundError: general__plasmidfinder`

**Solution:** Ensure you are using the custom `kleborate-runner.py` script from the **cloned repository**, not the system-wide `kleborate` command.

### PlasmidFinder not in PATH

**Error:** `Error: plasmidfinder.py not found in PATH`

**Solution:** Install PlasmidFinder and verify the conda environment is activated:
```bash
conda activate klebsiella_analysis
conda install -c bioconda plasmidfinder
```

### Invalid Parameter Values

**Error:** `Error: --plasmidfinder_mincov must be between 0 and 1`

**Solution:** Ensure threshold parameters are decimal values between 0 and 1:
```bash
# Correct
--plasmidfinder_mincov 0.60 --plasmidfinder_threshold 0.90

# Incorrect
--plasmidfinder_mincov 60 --plasmidfinder_threshold 90
```

### Custom Database Path Error

**Error:** `Error: PlasmidFinder database path does not exist`

**Solution:** Verify the custom database directory path exists and is accessible:
```bash
# Check if directory exists
ls -la /path/to/custom/db
```

## Related Resources

- [Kleborate GitHub Repository](https://github.com/tnmquann/Kleborate)
- [PlasmidFinder Documentation](https://bitbucket.org/genomicepidemiology/plasmidfinder/src/master/)
- [Kleborate Official Documentation](https://kleborate.readthedocs.io/en/latest/index.html)

---

**Last Updated:** 2025  
**Author:** Minh-Quan Ton-Ngoc  
**Repository:** https://github.com/tnmquann/Kleborate
