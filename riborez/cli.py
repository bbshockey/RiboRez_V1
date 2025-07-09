#!/usr/bin/env python3
"""
RiboRez Command Line Interface
"""

import argparse
import sys
from .download_taxa import download_taxa
from .gene_extract import extract_genes
from .primer_design import design_primers
from .amplicon_analysis import analyze_amplicons

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
  riborez download-taxa --taxon-name Pseudomonas --taxon-id 286 --max-genomes 200
  riborez download-taxa --taxon-name Pseudomonas --taxon-id 286 --max-genomes 200 --reference
  riborez gene-extract --taxon-name Pseudomonas
  riborez gene-extract --taxon-name Acinetobacter --sample-size 200
  riborez gene-extract --taxon-name Pseudomonas --genes 16S 23S
  riborez gene-extract --taxon-name Ecoli --genes rRNA
  riborez gene-extract --taxon-name Bacillus --genes gyrA recA
  riborez primer-design --input-folder Pseudomonas_AllGenesExtracted_rRNA
  riborez primer-design --input-folder Ecoli_genes --min-sequences 20 --threads 4
  riborez primer-design --input-folder Pseudomonas_genes --run-amplicon-analysis
  riborez amplicon-analysis --input-folder Pseudomonas_AllGenesExtracted_rRNA_Primers
  riborez amplicon-analysis --input-folder Ecoli_Primers --output-folder my_analysis --threads 8
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
    download_parser.add_argument(
        "--max-genomes", 
        type=int, 
        help="Maximum number of genomes to download (default: all available)"
    )
    download_parser.add_argument(
        "--reference",
        action="store_true",
        help="Restrict to reference genomes only (default: all RefSeq genomes)"
    )
    
    # Gene-extract subcommand
    gene_extract_parser = subparsers.add_parser(
        "gene-extract",
        help="Extract genes from downloaded NCBI datasets",
        description="Extract genes from downloaded genome datasets"
    )
    
    gene_extract_parser.add_argument(
        "--taxon-name", 
        required=True, 
        help="Taxon name (used to locate downloaded data directory)"
    )
    gene_extract_parser.add_argument(
        "--genes", 
        nargs="+", 
        help="Specific genes to extract (default: all genes). Examples: 16S, 23S, rRNA, gyrA, recA"
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
    
    # Primer-design subcommand
    primer_design_parser = subparsers.add_parser(
        "primer-design",
        help="Design primers for extracted genes using PMPrimer",
        description="Design primers for genes using PMPrimer with multiple parameter sets"
    )
    
    primer_design_parser.add_argument(
        "--input-folder", 
        required=True, 
        help="Input folder containing FASTA files (e.g., output from gene-extract)"
    )
    primer_design_parser.add_argument(
        "--output-folder", 
        help="Output directory for primer design results (auto-generated if not provided)"
    )
    primer_design_parser.add_argument(
        "--min-sequences", 
        type=int, 
        default=10, 
        help="Minimum number of sequences required per gene (default: 10)"
    )
    primer_design_parser.add_argument(
        "--threads", 
        type=int, 
        default=8, 
        help="Number of threads to use (default: 8)"
    )
    primer_design_parser.add_argument(
        "--run-amplicon-analysis", 
        action="store_true", 
        help="Automatically run amplicon analysis on the output folder after primer design"
    )
    primer_design_parser.add_argument(
        "--faster",
        action="store_true",
        help="Run only alignment and primary PMPrimer commands (skip parameter sweeps)"
    )
    
    # Amplicon-analysis subcommand
    amplicon_analysis_parser = subparsers.add_parser(
        "amplicon-analysis",
        help="Analyze amplicons from primer design results",
        description="Run comprehensive amplicon analysis workflow including amplification, mapping, and centralization"
    )
    
    amplicon_analysis_parser.add_argument(
        "--input-folder", 
        required=True, 
        help="Input folder containing primer design results (output from primer-design)"
    )
    amplicon_analysis_parser.add_argument(
        "--output-folder", 
        help="Output directory for analysis results (auto-generated if not provided)"
    )
    amplicon_analysis_parser.add_argument(
        "--threads", 
        type=int, 
        default=8, 
        help="Number of threads to use (default: 8)"
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
                dry_run=args.dry_run,
                max_genomes=args.max_genomes,
                reference=args.reference
            )
        elif args.command == "gene-extract":
            extract_genes(
                taxon_name=args.taxon_name,
                data_root=args.data_root,
                output_dir=args.output_dir,
                sample_size=args.sample_size,
                random_seed=args.random_seed,
                genes=args.genes
            )
        elif args.command == "primer-design":
            design_primers(
                input_folder=args.input_folder,
                output_folder=args.output_folder,
                min_sequences=args.min_sequences,
                threads=args.threads,
                run_amplicon_analysis=args.run_amplicon_analysis,
                faster=getattr(args, 'faster', False)
            )
        elif args.command == "amplicon-analysis":
            analyze_amplicons(
                input_folder=args.input_folder,
                output_folder=args.output_folder,
                threads=args.threads
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