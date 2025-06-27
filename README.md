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

Extract all genes (CDS and rRNA) from downloaded datasets:

```bash
# Extract all genes from all genomes
riborez gene-extract --taxon-name Pseudomonas

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