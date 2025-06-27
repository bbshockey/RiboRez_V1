# RiboRez

A tool for ribosomal analysis and data management.

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

### Command Options

- `--taxon-name`: Taxon name (used for output folder name) [required]
- `--taxon-id`: NCBI Taxon ID [required]
- `--output-dir`: Optional custom output directory
- `--rehydrate`: Rehydrate datasets (default)
- `--no-rehydrate`: Skip rehydration
- `--force`: Overwrite output directory if it exists
- `--dry-run`: Print commands without executing

## Prerequisites

This tool requires the NCBI Datasets command-line tool to be installed:

```bash
# Install NCBI Datasets
pip install ncbi-datasets-cli
```

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