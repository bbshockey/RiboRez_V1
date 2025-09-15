# RiboRez

A tool for exploring mRNA transcripts for taxonomic resolution. 

![RiboRez workflow overview](https://raw.githubusercontent.com/bbshockey/RiboRez_V1/main/assets/Riborez_workflow.jpeg)

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
- `--reference`: Restrict to reference genomes only (default: all RefSeq genomes)

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
- `--faster`: Run only the alignment and primary PMPrimer commands, skipping the full parameter sweep. This significantly speeds up the process, but may reduce the diversity of primer candidates.

## COMPLETE WORKFLOW EXAMPLE

Here's a complete example showing a typical workflow starting with downloading all available Pseudomonas genomes via ncbi datasets, to the extraction of all 16S rRNA from those genomes, to universal primer design and amplicon analysis of the 16S gene. 
Special Case: 'If you already have a gene fasta file and don't want to extract from ncbi.' Put your fasta file in a parent folder and run the command primer-design with your --input-folder being the parent folder location.

```bash
# 1. Download genomes for a taxon
riborez download-taxa --taxon-name Pseudomonas --taxon-id 286

# 2. Extract specific genes (e.g., 16S rRNA)
riborez gene-extract --taxon-name Pseudomonas --genes 16S

# 3. Design primers AND analyze amplicons in one step
riborez primer-design --input-folder Pseudomonas_TranscriptsExtracted --run-amplicon-analysis
```

## Tool Setup and Prerequisites

### Environment Setup
```bash
conda create -n env_name python=3.8
conda install -c bioconda muscle
```

### NCBI Datasets CLI Installation
This tool requires the NCBI Datasets command-line tool. You can install it using any of these methods (see the [official installation guide](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/) for more details):

#### Automatic Installation
```bash
python install_ncbi_datasets.py
```

#### Manual Installation Options
- **Using pip:**
  ```bash
  pip install ncbi-datasets-cli
  ```
- **Using conda:**
  ```bash
  conda install -c conda-forge ncbi-datasets-cli
  ```

### RiboRez Installation

#### From Source
```bash
git clone https://github.com/bbshockey/RiboRez_V1.git
cd RiboRez_V1
pip install -e .
```

#### From GitHub
```bash
pip install git+https://github.com/bbshockey/RiboRez_V1.git
```

