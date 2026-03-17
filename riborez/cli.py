#!/usr/bin/env python3
"""
RiboRez Command Line Interface
"""

import argparse
import sys
from .download_taxa import download_taxa, download_taxa_multi
from .gene_extract import extract_genes
from .primer_design import design_primers
from .amplicon_analysis import analyze_amplicons
from .ribozyme_design import design_ribozymes

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
  riborez primer-design --input-folder Pseudomonas_RNAextracted
  riborez primer-design --input-folder Ecoli_RNAextracted --min-sequences 20 --threads 4
  riborez primer-design --input-folder Pseudomonas_RNAextracted --run-amplicon-analysis
  riborez amplicon-analysis --input-folder Pseudomonas_Primers
  riborez amplicon-analysis --input-folder Ecoli_Primers --output-folder my_analysis --threads 8
  riborez run --taxon-name Pseudomonas --taxon-id 286 --genes 16S 23S --threads 22
  riborez run --taxon-name Ecoli --taxon-id 562 --genes rRNA --max-genomes 200 --assembly-level complete --threads 22 --faster
  riborez run --taxon-name Pseudomonas --taxon-id 286 --skip-download --genes 16S --threads 22
  riborez ribozyme-design --input-folder Pseudomonas_Primers/16S
  riborez ribozyme-design --input-folder Pseudomonas_Primers/acnB --max-igs-mismatches 10
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
        help="NCBI Taxon ID or comma-separated list of IDs (e.g., 286,590)"
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
        help="Restrict to NCBI-designated reference genomes only (very few, highly curated). For broader high-quality filtering, use --assembly-level instead."
    )
    download_parser.add_argument(
        "--assembly-level",
        choices=["complete", "chromosome", "scaffold", "contig"],
        help=(
            "Filter by assembly quality level. "
            "'complete' = fully assembled, no gaps (recommended for well-curated sets). "
            "'chromosome' = assembled to chromosome level. "
            "'scaffold' = scaffold-level assembly. "
            "'contig' = raw contigs (lowest quality). "
            "Default: all levels included."
        )
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
        help="Must match the name used in download-taxa. Used to locate the downloaded data folder (e.g., 'Pseudomonas' looks for Pseudomonas_NCBI/). Use --data-root if you renamed the folder."
    )
    gene_extract_parser.add_argument(
        "--genes",
        nargs="+",
        help="Specific genes to extract (default: all genes). Examples: 16S, 23S, rRNA, gyrA, recA"
    )
    gene_extract_parser.add_argument(
        "--data-root",
        help=(
            "Path to your genome data directory. "
            "Supports NCBI dataset format (subdirectories per genome, each with FASTA+GFF) "
            "or flat format (FASTA and GFF files placed directly in the directory, paired by matching filename stem). "
            "Use this when you renamed the download folder or are working with your own genome files."
        )
    )
    gene_extract_parser.add_argument(
        "--output-dir", 
        help="Output directory for extracted genes (auto-generated if not provided)"
    )
    gene_extract_parser.add_argument(
        "--min-per-gene",
        type=int,
        default=5,
        help="Minimum sequences required to write a gene FASTA (default: 5)"
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

    # Ribozyme-design subcommand
    ribozyme_parser = subparsers.add_parser(
        "ribozyme-design",
        help="Design group-I-intron ribozymes paired with primer pairs",
        description=(
            "Find group-I-intron T-sites near each reverse primer in a gene directory, "
            "score IGS/EGS conservation, and produce a combined ranking of primer pairs "
            "that maximize both taxonomic resolution and ribozyme cleavage reliability.\n\n"
            "Requires 'riborez primer-design' output. For combined ranking, also requires "
            "'riborez amplicon-analysis' output (amplicon.summary.csv or *_processed_amplicon.csv)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    ribozyme_parser.add_argument(
        "--input-folder",
        required=True,
        help="Single gene directory from primer-design output (e.g., Pseudomonas_Primers/16S/)"
    )
    ribozyme_parser.add_argument(
        "--output-folder",
        help="Output directory (auto-named {gene}_RibozymeDesign/ next to input folder if not provided)"
    )
    ribozyme_parser.add_argument(
        "--window",
        type=int,
        default=10,
        help="bp to extend T-site search beyond reverse primer region (default: 10)"
    )
    ribozyme_parser.add_argument(
        "--egs-start",
        type=int,
        default=11,
        help="EGS region start offset from T-site (default: 11)"
    )
    ribozyme_parser.add_argument(
        "--egs-end",
        type=int,
        default=60,
        help="EGS region end offset from T-site (default: 60)"
    )
    ribozyme_parser.add_argument(
        "--p1-loop",
        default="TAACCACA",
        help="Fixed P1 loop sequence (default: TAACCACA)"
    )
    ribozyme_parser.add_argument(
        "--max-amplicon-length",
        type=int,
        default=500,
        help="Discard combined ranking candidates with amplicon >= this bp (default: 500)"
    )
    ribozyme_parser.add_argument(
        "--max-igs-mismatches",
        type=int,
        default=None,
        help="Discard candidates where total IGS mismatches across all sequences exceeds this (default: no filter)"
    )

    # Run (full pipeline) subcommand
    run_parser = subparsers.add_parser(
        "run",
        help="Run the full pipeline: download → extract → primer-design → amplicon-analysis",
        description=(
            "Execute all four RiboRez steps in sequence for a single taxon.\n"
            "Outputs follow the standard naming convention:\n"
            "  {taxon-name}_NCBI/  →  {taxon-name}_RNAextracted/  →  {taxon-name}_Primers/  →  {taxon-name}_AmpliconAnalysis/"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # --- Required ---
    run_parser.add_argument(
        "--taxon-name",
        required=True,
        help="Taxon name used for all output folder names (e.g., 'Pseudomonas')"
    )
    run_parser.add_argument(
        "--taxon-id",
        required=True,
        help="NCBI Taxon ID or comma-separated list of IDs (e.g., 286 or 286,590)"
    )

    # --- Download options ---
    run_parser.add_argument(
        "--max-genomes",
        type=int,
        help="Maximum number of genomes to download (default: all available)"
    )
    run_parser.add_argument(
        "--assembly-level",
        choices=["complete", "chromosome", "scaffold", "contig"],
        help="Filter downloads by assembly quality level (recommended: 'complete')"
    )
    run_parser.add_argument(
        "--reference",
        action="store_true",
        help="Restrict download to NCBI-designated reference genomes only (very few per species)"
    )
    run_parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing download directory if it exists"
    )

    # --- Extract options ---
    run_parser.add_argument(
        "--genes",
        nargs="+",
        help="Genes to extract (default: all). Examples: 16S 23S rRNA gyrA"
    )
    run_parser.add_argument(
        "--min-per-gene",
        type=int,
        default=5,
        help="Minimum sequences required to write a gene FASTA (default: 5)"
    )
    run_parser.add_argument(
        "--sample-size",
        type=int,
        help="Number of genomes to sample for extraction (default: all)"
    )

    # --- Primer design options ---
    run_parser.add_argument(
        "--min-sequences",
        type=int,
        default=10,
        help="Minimum sequences per gene required for primer design (default: 10)"
    )
    run_parser.add_argument(
        "--faster",
        action="store_true",
        help="Skip full PMPrimer parameter sweep (faster but fewer primer candidates)"
    )

    # --- Shared ---
    run_parser.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Number of threads for primer design and amplicon analysis (default: 8)"
    )

    # --- Skip flags ---
    run_parser.add_argument(
        "--skip-download",
        action="store_true",
        help="Skip the download-taxa step (use an existing {taxon-name}_NCBI/ folder)"
    )
    run_parser.add_argument(
        "--skip-extract",
        action="store_true",
        help="Skip the gene-extract step (use an existing {taxon-name}_RNAextracted/ folder)"
    )
    run_parser.add_argument(
        "--skip-primer-design",
        action="store_true",
        help="Skip the primer-design step (use an existing {taxon-name}_Primers/ folder)"
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
            # Support single ID or comma-separated list
            try:
                taxon_ids = [int(x.strip()) for x in str(args.taxon_id).split(',') if x.strip()]
            except ValueError:
                raise ValueError("--taxon-id must be an integer or a comma-separated list of integers")

            if len(taxon_ids) == 1:
                download_taxa(
                    taxon_name=args.taxon_name,
                    taxon_id=taxon_ids[0],
                    output_dir=args.output_dir,
                    rehydrate=args.rehydrate,
                    force=args.force,
                    dry_run=args.dry_run,
                    max_genomes=args.max_genomes,
                    reference=args.reference,
                    assembly_level=args.assembly_level
                )
            else:
                download_taxa_multi(
                    taxon_name=args.taxon_name,
                    taxon_ids=taxon_ids,
                    output_dir=args.output_dir,
                    rehydrate=args.rehydrate,
                    force=args.force,
                    dry_run=args.dry_run,
                    max_genomes=args.max_genomes,
                    reference=args.reference,
                    assembly_level=args.assembly_level
                )
        elif args.command == "gene-extract":
            extract_genes(
                taxon_name=args.taxon_name,
                data_root=args.data_root,
                output_dir=args.output_dir,
                min_per_gene=args.min_per_gene,
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
        elif args.command == "ribozyme-design":
            design_ribozymes(
                input_folder=args.input_folder,
                output_folder=args.output_folder,
                window=args.window,
                egs_start=args.egs_start,
                egs_end=args.egs_end,
                p1_loop=args.p1_loop,
                max_amplicon_length=args.max_amplicon_length,
                max_igs_mismatches=args.max_igs_mismatches,
            )
        elif args.command == "run":
            try:
                taxon_ids = [int(x.strip()) for x in str(args.taxon_id).split(',') if x.strip()]
            except ValueError:
                raise ValueError("--taxon-id must be an integer or a comma-separated list of integers")

            # Step 1: Download
            if args.skip_download:
                print(f"\n[SKIP] Step 1/4: download-taxa (--skip-download)")
            else:
                print(f"\n{'='*60}")
                print(f"  Step 1/4: download-taxa")
                print(f"{'='*60}")
                if len(taxon_ids) == 1:
                    download_taxa(
                        taxon_name=args.taxon_name,
                        taxon_id=taxon_ids[0],
                        output_dir=None,
                        rehydrate=True,
                        force=args.force,
                        dry_run=False,
                        max_genomes=args.max_genomes,
                        reference=args.reference,
                        assembly_level=args.assembly_level,
                    )
                else:
                    download_taxa_multi(
                        taxon_name=args.taxon_name,
                        taxon_ids=taxon_ids,
                        output_dir=None,
                        rehydrate=True,
                        force=args.force,
                        dry_run=False,
                        max_genomes=args.max_genomes,
                        reference=args.reference,
                        assembly_level=args.assembly_level,
                    )

            # Step 2: Gene extraction
            if args.skip_extract:
                print(f"\n[SKIP] Step 2/4: gene-extract (--skip-extract)")
            else:
                print(f"\n{'='*60}")
                print(f"  Step 2/4: gene-extract")
                print(f"{'='*60}")
                extract_output = extract_genes(
                    taxon_name=args.taxon_name,
                    data_root=None,
                    output_dir=None,
                    min_per_gene=args.min_per_gene,
                    sample_size=args.sample_size,
                    random_seed=42,
                    genes=args.genes,
                )

            # Step 3: Primer design
            if args.skip_primer_design:
                print(f"\n[SKIP] Step 3/4: primer-design (--skip-primer-design)")
            else:
                print(f"\n{'='*60}")
                print(f"  Step 3/4: primer-design")
                print(f"{'='*60}")
                extract_folder = f"{args.taxon_name}_RNAextracted"
                primers_output = design_primers(
                    input_folder=extract_folder,
                    output_folder=None,
                    min_sequences=args.min_sequences,
                    threads=args.threads,
                    run_amplicon_analysis=False,
                    faster=args.faster,
                )

            # Step 4: Amplicon analysis
            print(f"\n{'='*60}")
            print(f"  Step 4/4: amplicon-analysis")
            print(f"{'='*60}")
            primers_folder = f"{args.taxon_name}_Primers"
            analyze_amplicons(
                input_folder=primers_folder,
                output_folder=None,
                threads=args.threads,
            )

            print(f"\n{'='*60}")
            print(f"  Full pipeline complete for: {args.taxon_name}")
            print(f"{'='*60}")
        else:
            print(f"Unknown command: {args.command}")
            parser.print_help()
            sys.exit(1)
            
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main() 