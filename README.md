# RiboRez

A tool for exploring mRNA transcripts for taxonomic resolution. 

## PREREQUISITES 
```bash
# prerequisites
conda create -n env_name python=3.8
conda install -c bioconda muscle
```

This tool also requires the NCBI Datasets command-line tool to be installed. You can install it using any of these methods (see the [official installation guide](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/) for more details):

## NCBI Automatic Installation
```bash
# Run the automatic installation script
python install_ncbi_datasets.py
```

## NCBI Manual Installation Options

#### Option 1: Using pip
```bash
pip install ncbi-datasets-cli
```

#### Option 2: Using conda
```bash
conda install -c conda-forge ncbi-datasets-cli
```


## RIBOREZ INSTALLATION

### From Source
```bash
git clone https://github.com/bbshockey/RiboRez_V1.git
cd RiboRez_V1
pip install -e .
```

### From GitHub
```bash
pip install git+https://github.com/bbshockey/RiboRez_V1.git
```

**Note:** Make sure to install the NCBI Datasets CLI first (see Prerequisites above).

## USAGE

### Download Taxa Data

Download and unpack NCBI genomes for a specific taxon. The taxon ID must be obtained from the NCBI Taxonomy database (https://www.ncbi.nlm.nih.gov/taxonomy). The taxon name is used for organizing output files and can be any name that helps you identify the dataset:

```bash
# Basic usage
riborez download-taxa --taxon-name Pseudomonas --taxon-id 286

# Force overwrite existing directory (will only trigger this with same taxon-name parameter)
riborez download-taxa --taxon-name Salmonella --taxon-id 590 --force

# Dry run (prints ncbi datasets commands without executing)
riborez download-taxa --taxon-name Bacillus --taxon-id 1386 --dry-run

# Download only 200 genomes (useful for testing or limited resources)
riborez download-taxa --taxon-name Pseudomonas --taxon-id 286 --max-genomes 200
```

### EXTRACT GENES

Extract genes from downloaded datasets:
The user should match their taxon-name from their riborez download-taxa subcommand. If using a custom set of gffs and fastas they should match the data structure of the datasets download. Example structure: taxon_name_NCBI/genomes/ncbi_dataset/data/GCA_043082705.1/genomic.gff and taxon_name_NCBI/genomes/ncbi_dataset/data/GCA_043082705.1/GCA_043082705.1_PDT002164297.1_genomic.fna -This input would be -taxon-name taxon_name.

```bash
# Extract all genes from all genomes 
riborez gene-extract --taxon-name Pseudomonas

# Extract specific rRNA genes
riborez gene-extract --taxon-name Pseudomonas --genes 16S 23S

# Extract genes from a sample of 200 genomes
riborez gene-extract --taxon-name Acinetobacter --sample-size 200

# Extract multiple specific genes from a sample of 100 genomes using a repeatable seed.
riborez gene-extract --taxon-name Bacillus --genes 16S recA tuf ccoN --sample-size 100 --random-seed 123
```

### PRIMER AND AMPLICON ANALYSIS

Design universal primers for each gene using the PMPrimer python tool. Then analyze the resulting amplicons for differentiation. 
```bash
# Basic primer design with default settings
riborez primer-design --input-folder Pseudomonas_AllGenesExtracted_rRNA

# Primer design with custom minimum sequences and thread count
riborez primer-design --input-folder Ecoli_genes --min-sequences 20 --threads 4

# Primer design with automatic amplicon analysis (recommended workflow)
riborez primer-design --input-folder Pseudomonas_genes --run-amplicon-analysis
```

### Command Options

#### download-taxa
- `--taxon-name`: Taxon name (used for output folder name) [required]
- `--taxon-id`: NCBI Taxon ID [required]
- `--output-dir`: Optional custom output directory
- `--rehydrate`: Rehydrate datasets (default)
- `--no-rehydrate`: Skip rehydration
- `--force`: Overwrite output directory if it exists
- `--dry-run`: Print commands without executing
- `--max-genomes`: Maximum number of genomes to download (default: all available)

#### gene-extract
- `--taxon-name`: Taxon name (used to locate downloaded data directory) [required]
- `--genes`: Specific genes to extract (default: all genes). Examples: 16S, 23S, rRNA, gyrA, recA
- `--data-root`: Path to data directory (auto-detected if not provided)
- `--output-dir`: Output directory for extracted genes (auto-generated if not provided)
- `--sample-size`: Number of genomes to sample (default: all available)
- `--random-seed`: Random seed for sampling (default: 42)

#### primer-design
- `--input-folder`: Input folder containing FASTA files (e.g., output from gene-extract) [required]
- `--output-folder`: Output directory for primer design results (auto-generated if not provided)
- `--min-sequences`: Minimum number of sequences required per gene (default: 10)
- `--threads`: Number of threads to use (default: 8)
- `--run-amplicon-analysis`: Automatically run amplicon analysis on the output folder after primer design

## COMPLETE WORKFLOW EXAMPLE

Here's a complete example showing a typical workflow starting with downloading all available Pseudomonas genomes via ncbi datasets, to the extraction of all 16S rRNA from those genomes, to universal primer design and amplicon analysis of the 16S gene. 

```bash
# 1. Download genomes for a taxon
riborez download-taxa --taxon-name Pseudomonas --taxon-id 286

# 2. Extract specific genes (e.g., 16S rRNA)
riborez gene-extract --taxon-name Pseudomonas --genes 16S

# 3. Design primers AND analyze amplicons in one step
riborez primer-design --input-folder Pseudomonas_AllGenesExtracted_rRNA --run-amplicon-analysis
```
