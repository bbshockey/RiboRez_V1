# RiboRez

A tool for ribosomal analysis and data management.

### Prerequisites
```bash
# prerequisites
conda create -n env_name python=3.8
```

This tool requires the NCBI Datasets command-line tool to be installed. You can install it using any of these methods:

### Automatic Installation (Recommended)
```bash
# Run the automatic installation script
python install_ncbi_datasets.py
```

### Manual Installation Options

#### Option 1: Using pip
```bash
pip install ncbi-datasets-cli
```

#### Option 2: Using conda
```bash
conda install -c conda-forge ncbi-datasets-cli
```

#### Option 3: Direct download
**Linux:**
```bash
curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
curl -o dataformat 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat'
chmod +x datasets dataformat
```

**macOS:**
```bash
curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/mac/datasets'
curl -o dataformat 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/mac/dataformat'
chmod +x datasets dataformat
```

**Windows:**
```bash
curl -o datasets.exe "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/win64/datasets.exe"
curl -o dataformat.exe "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/win64/dataformat.exe"
```

## Installation

### From Source
```bash
git clone https://github.com/bbshockey/RiboRez_V1.git
cd RiboRez_V1
pip install -e .
```

### From GitHub (after publishing)
```bash
pip install git+https://github.com/bbshockey/RiboRez_V1.git
```

**Note:** Make sure to install the NCBI Datasets CLI first (see Prerequisites above).

## Usage

### Download Taxa Data

Download and unpack NCBI genomes for a specific taxon:

```bash
# Basic usage
riborez download-taxa --taxon-name Pseudomonas --taxon-id 286

# With custom output directory
riborez download-taxa --taxon-name Ecoli --taxon-id 562 --output-dir my_data

# Force overwrite existing directory
riborez download-taxa --taxon-name Salmonella --taxon-id 590 --force

# Dry run (print commands without executing)
riborez download-taxa --taxon-name Bacillus --taxon-id 1386 --dry-run

# Skip rehydration
riborez download-taxa --taxon-name Staphylococcus --taxon-id 1279 --no-rehydrate
```

### Extract Genes

Extract genes from downloaded datasets:

```bash
# Extract all genes from all genomes
riborez gene-extract --taxon-name Pseudomonas

# Extract specific rRNA genes
riborez gene-extract --taxon-name Pseudomonas --genes 16S 23S

# Extract all rRNA genes
riborez gene-extract --taxon-name Pseudomonas --genes rRNA

# Extract specific protein-coding genes
riborez gene-extract --taxon-name Ecoli --genes gyrA recA

# Extract genes from a sample of 200 genomes
riborez gene-extract --taxon-name Acinetobacter --sample-size 200

# Custom output directory
riborez gene-extract --taxon-name Ecoli --output-dir my_extracted_genes

# Custom data root path
riborez gene-extract --taxon-name Salmonella --data-root /path/to/data

# Custom random seed for sampling
riborez gene-extract --taxon-name Bacillus --sample-size 100 --random-seed 123
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

#### gene-extract
- `--taxon-name`: Taxon name (used to locate downloaded data directory) [required]
- `--genes`: Specific genes to extract (default: all genes). Examples: 16S, 23S, rRNA, gyrA, recA
- `--data-root`: Path to data directory (auto-detected if not provided)
- `--output-dir`: Output directory for extracted genes (auto-generated if not provided)
- `--sample-size`: Number of genomes to sample (default: all available)
- `--random-seed`: Random seed for sampling (default: 42)

## Development

To set up the development environment:

```bash
git clone https://github.com/yourusername/RiboRez_V1.git
cd RiboRez_V1
pip install -e .
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## Issues

If you encounter any issues, please report them on the [GitHub issues page](https://github.com/bbshockey/RiboRez_V1/issues).

## Available Commands

### 1. `download-taxa`
Download and unpack NCBI genomes for a specific taxon.

**Required arguments:**
- `--taxon-name`: Taxon name (used for output folder name)
- `--taxon-id`: NCBI Taxon ID

**Optional arguments:**
- `--output-dir`: Custom output directory
- `--force`: Overwrite output directory if it exists
- `--dry-run`: Print commands without executing
- `--rehydrate`/`--no-rehydrate`: Control dataset rehydration (default: rehydrate)

**Examples:**
```bash
riborez download-taxa --taxon-name Pseudomonas --taxon-id 286
riborez download-taxa --taxon-name Ecoli --taxon-id 562 --output-dir my_data --force
riborez download-taxa --taxon-name Salmonella --taxon-id 590 --dry-run
```

### 2. `gene-extract`
Extract genes from downloaded NCBI datasets.

**Required arguments:**
- `--taxon-name`: Taxon name (used to locate downloaded data directory)

**Optional arguments:**
- `--genes`: Specific genes to extract (default: all genes). Examples: 16S, 23S, rRNA, gyrA, recA
- `--data-root`: Path to data directory (auto-detected if not provided)
- `--output-dir`: Output directory for extracted genes (auto-generated if not provided)
- `--sample-size`: Number of genomes to sample (default: all available)
- `--random-seed`: Random seed for sampling (default: 42)

**Examples:**
```bash
riborez gene-extract --taxon-name Pseudomonas
riborez gene-extract --taxon-name Acinetobacter --sample-size 200
riborez gene-extract --taxon-name Pseudomonas --genes 16S 23S
riborez gene-extract --taxon-name Ecoli --genes rRNA
riborez gene-extract --taxon-name Bacillus --genes gyrA recA
```

### 3. `primer-design`
Design primers for extracted genes using PMPrimer.

**Required arguments:**
- `--input-folder`: Input folder containing FASTA files (e.g., output from gene-extract)

**Optional arguments:**
- `--output-folder`: Output directory for primer design results (auto-generated if not provided)
- `--min-sequences`: Minimum number of sequences required per gene (default: 10)
- `--threads`: Number of threads to use (default: 8)
- `--run-amplicon-analysis`: Automatically run amplicon analysis on the output folder after primer design

**Examples:**
```bash
riborez primer-design --input-folder Pseudomonas_AllGenesExtracted_rRNA
riborez primer-design --input-folder Ecoli_genes --min-sequences 20 --threads 4
riborez primer-design --input-folder Pseudomonas_genes --run-amplicon-analysis
```

### 4. `amplicon-analysis`
Analyze amplicons from primer design results using a comprehensive workflow.

**Required arguments:**
- `--input-folder`: Input folder containing primer design results (output from primer-design)

**Optional arguments:**
- `--output-folder`: Output directory for analysis results (auto-generated if not provided)
- `--threads`: Number of threads to use (default: 8)

**Examples:**
```bash
riborez amplicon-analysis --input-folder Pseudomonas_AllGenesExtracted_rRNA_Primers
riborez amplicon-analysis --input-folder Ecoli_Primers --output-folder my_analysis --threads 8
```

## Complete Workflow Example

Here's a complete example of the typical workflow:

```bash
# 1. Download genomes for a taxon
riborez download-taxa --taxon-name Pseudomonas --taxon-id 286

# 2. Extract specific genes (e.g., 16S rRNA)
riborez gene-extract --taxon-name Pseudomonas --genes 16S

# 3. Design primers for the extracted genes
riborez primer-design --input-folder Pseudomonas_AllGenesExtracted_rRNA

# 4. Analyze the amplicons comprehensively
riborez amplicon-analysis --input-folder Pseudomonas_AllGenesExtracted_rRNA_Primers
```

## Output Directory Structure

### Gene Extraction Output
```
Pseudomonas_AllGenesExtracted_rRNA/
├── 16S.fasta
├── 23S.fasta
├── extraction_log.txt
└── temp/
```

### Primer Design Output
```
Pseudomonas_AllGenesExtracted_rRNA_Primers/
├── 16S/
│   ├── 16S.fasta
│   ├── 16S.filt.mc.fasta
│   ├── 16S_reference_mapping.tsv
│   └── [primer design files]
├── 23S/
│   └── [similar structure]
├── reference_mappings/
└── pmprimer_log.log
```

### Amplicon Analysis Output
```
Pseudomonas_AllGenesExtracted_rRNA_Primers_AmpliconAnalysis/
├── Pseudomonas_AllGenesExtracted_rRNA_Primers_best_amplicon_summary.csv
└── amplicon_analysis_log.log
```