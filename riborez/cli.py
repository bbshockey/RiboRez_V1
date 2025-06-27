#!/usr/bin/env python3
"""
RiboRez Command Line Interface
"""

import argparse
import sys
from .download_taxa import download_taxa
from .gene_extract import extract_genes

def main():
    """Main entry point for the riborez command."""
    parser = argparse.ArgumentParser(
        prog="riborez",
        description="RiboRez - A tool for ribosomal analysis and data management",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  riborez download-taxa --taxon-name Pseudomonas --taxon-id 286
  riborez download-taxa --taxon-name Ecoli --taxon-id 562 --output-dir my_data --force
  riborez download-taxa --taxon-name Salmonella --taxon-id 590 --dry-run
  riborez gene-extract --taxon-name Pseudomonas
  riborez gene-extract --taxon-name Acinetobacter --sample-size 200
        """
    )
    
    subparsers = parser.add_subparsers(
        dest="command", 
        help="Available commands",
        metavar="COMMAND"
    )
    
    # Download-taxa subcommand
    download_parser = subparsers.add_parser(
        "download-taxa",
        help="Download and unpack NCBI genomes for a taxon",
        description="Download genome data from NCBI for a specific taxon ID"
    )
    
    download_parser.add_argument(
        "--taxon-name", 
        required=True, 
        help="Taxon name (used for output folder name)"
    )
    download_parser.add_argument(
        "--taxon-id", 
        required=True, 
        type=int, 
        help="NCBI Taxon ID"
    )
    download_parser.add_argument(
        "--output-dir", 
        help="Optional custom output directory"
    )
    download_parser.add_argument(
        "--rehydrate", 
        dest="rehydrate", 
        action="store_true", 
        help="Rehydrate datasets (default)"
    )
    download_parser.add_argument(
        "--no-rehydrate", 
        dest="rehydrate", 
        action="store_false", 
        help="Skip rehydration"
    )
    download_parser.set_defaults(rehydrate=True)
    download_parser.add_argument(
        "--force", 
        action="store_true", 
        help="Overwrite output directory if it exists"
    )
    download_parser.add_argument(
        "--dry-run", 
        action="store_true", 
        help="Print commands without executing"
    )
    
    # Gene-extract subcommand
    gene_extract_parser = subparsers.add_parser(
        "gene-extract",
        help="Extract genes from downloaded NCBI datasets",
        description="Extract all genes (CDS and rRNA) from downloaded genome datasets"
    )
    
    gene_extract_parser.add_argument(
        "--taxon-name", 
        required=True, 
        help="Taxon name (used to locate downloaded data directory)"
    )
    gene_extract_parser.add_argument(
        "--data-root", 
        help="Path to data directory (auto-detected if not provided)"
    )
    gene_extract_parser.add_argument(
        "--output-dir", 
        help="Output directory for extracted genes (auto-generated if not provided)"
    )
    gene_extract_parser.add_argument(
        "--sample-size", 
        type=int, 
        help="Number of genomes to sample (default: all available)"
    )
    gene_extract_parser.add_argument(
        "--random-seed", 
        type=int, 
        default=42, 
        help="Random seed for sampling (default: 42)"
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    # Handle no command provided
    if not args.command:
        parser.print_help()
        sys.exit(1)
    
    # Execute commands
    try:
        if args.command == "download-taxa":
            download_taxa(
                taxon_name=args.taxon_name,
                taxon_id=args.taxon_id,
                output_dir=args.output_dir,
                rehydrate=args.rehydrate,
                force=args.force,
                dry_run=args.dry_run
            )
        elif args.command == "gene-extract":
            extract_genes(
                taxon_name=args.taxon_name,
                data_root=args.data_root,
                output_dir=args.output_dir,
                sample_size=args.sample_size,
                random_seed=args.random_seed
            )
        else:
            print(f"Unknown command: {args.command}")
            parser.print_help()
            sys.exit(1)
            
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main() 