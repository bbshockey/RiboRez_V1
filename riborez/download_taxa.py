import os
import shutil
import subprocess
import sys
from pathlib import Path
from .install_dependencies import ensure_dependencies

def check_datasets_available():
    """Check if the datasets command is available."""
    try:
        result = subprocess.run(['datasets', '--version'], 
                              capture_output=True, text=True)
        return result.returncode == 0
    except FileNotFoundError:
        return False

def run_command(cmd_list, dry_run):
    print("[CMD]", ' '.join(cmd_list))
    if not dry_run:
        subprocess.run(cmd_list, check=True)

def download_taxa(taxon_name, taxon_id, output_dir, rehydrate, force, dry_run, max_genomes=None, reference=False):
    """
    Download and unpack NCBI genomes for a taxon.
    
    Args:
        taxon_name (str): Taxon name (used for output folder name)
        taxon_id (int): NCBI Taxon ID
        output_dir (str): Optional custom output directory
        rehydrate (bool): Whether to rehydrate datasets
        force (bool): Overwrite output directory if it exists
        dry_run (bool): Print commands without executing
        max_genomes (int): Maximum number of genomes to download (None for all)
        reference (bool): Restrict to reference genomes only (default: False)
    """
    # Check if datasets command is available and install if needed
    if not check_datasets_available():
        print("[INFO] NCBI Datasets CLI not found. Attempting automatic installation...")
        if not ensure_dependencies():
            print("[ERROR] Failed to install NCBI Datasets CLI automatically")
            print("[INFO] Please install it manually using one of these methods:")
            print("  1. Run: python install_ncbi_datasets.py")
            print("  2. Install via conda: conda install -c conda-forge ncbi-datasets-cli")
            print("  3. Download manually from: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/")
            sys.exit(1)
        print("[SUCCESS] NCBI Datasets CLI installed successfully!")
    
    if not output_dir:
        output_dir = f"{taxon_name}_NCBI"

    output_path = Path(output_dir)
    zip_path = output_path / "genomes.zip"
    unzip_dir = output_path / "genomes"

    if output_path.exists():
        if force:
            print(f"[INFO] Removing existing directory: {output_path}")
            if not dry_run:
                shutil.rmtree(output_path)
        else:
            raise FileExistsError(f"Output directory {output_path} exists. Use --force to overwrite.")

    print(f"[INFO] Creating output directory: {output_path}")
    if not dry_run:
        output_path.mkdir(parents=True, exist_ok=True)

    if max_genomes:
        print(f"[INFO] Limiting download to {max_genomes} genomes")
        
        # Step 1: Get summary with limit to get accession IDs
        accessions_file = output_path / "accessions.txt"
        accessions_clean_file = output_path / "accessions_clean.txt"
        
        summary_cmd = [
            "datasets", "summary", "genome", "taxon", str(taxon_id),
            "--limit", str(max_genomes),
            "--assembly-source", "RefSeq",
            "--report", "ids_only",
            "--as-json-lines"
        ]
        if reference:
            summary_cmd.append("--reference")
        
        # Run summary and pipe to dataformat
        if not dry_run:
            with open(accessions_file, 'w') as f:
                subprocess.run(summary_cmd, stdout=f, check=True)
            
            # Convert to TSV format
            dataformat_cmd = ["dataformat", "tsv", "genome", "--fields", "accession"]
            with open(accessions_file, 'r') as infile, open(accessions_clean_file, 'w') as outfile:
                subprocess.run(dataformat_cmd, stdin=infile, stdout=outfile, check=True)
            
            # Remove header line (tail -n +2 equivalent)
            with open(accessions_clean_file, 'r') as f:
                lines = f.readlines()
            with open(accessions_clean_file, 'w') as f:
                f.writelines(lines[1:])  # Skip first line (header)
        else:
            print(f"[CMD] {' '.join(summary_cmd)} | dataformat tsv genome --fields accession > {accessions_file}")
            print(f"[CMD] tail -n +2 {accessions_file} > {accessions_clean_file}")
        
        # Step 2: Download genomes using accession file
        download_cmd = [
            "datasets", "download", "genome", "accession",
            "--inputfile", str(accessions_clean_file),
            "--dehydrated",
            "--include", "genome,gff3",
            "--filename", str(zip_path)
        ]
        run_command(download_cmd, dry_run)
        
        # Clean up temporary files
        if not dry_run:
            accessions_file.unlink(missing_ok=True)
            accessions_clean_file.unlink(missing_ok=True)
    else:
        # Original workflow: download all genomes
        print("[INFO] Downloading all available genomes")
        download_cmd = [
            "datasets", "download", "genome", "taxon", str(taxon_id),
            "--dehydrated",
            "--include", "genome,gff3",
            "--filename", str(zip_path)
        ]
        run_command(download_cmd, dry_run)

    # Step 3: Unzip
    unzip_cmd = ["unzip", "-o", str(zip_path), "-d", str(unzip_dir)]
    run_command(unzip_cmd, dry_run)

    # Step 4: Rehydrate (optional)
    if rehydrate:
        rehydrate_cmd = ["datasets", "rehydrate", "--directory", str(unzip_dir)]
        run_command(rehydrate_cmd, dry_run) 


def download_taxa_multi(taxon_name, taxon_ids, output_dir, rehydrate, force, dry_run, max_genomes=None, reference=False):
    """Download multiple taxon IDs into the same output directory.

    Each taxon ID is downloaded in sequence into the same zip and unzip locations.
    The first ID will create the directory (respecting --force); subsequent IDs
    will reuse the existing directory and append contents.
    """
    # Ensure dependencies once
    if not check_datasets_available():
        print("[INFO] NCBI Datasets CLI not found. Attempting automatic installation...")
        if not ensure_dependencies():
            print("[ERROR] Failed to install NCBI Datasets CLI automatically")
            sys.exit(1)
        print("[SUCCESS] NCBI Datasets CLI installed successfully!")

    if not output_dir:
        output_dir = f"{taxon_name}_NCBI"

    output_path = Path(output_dir)
    zip_path = output_path / "genomes.zip"
    unzip_dir = output_path / "genomes"

    # For multi, only clear directory once if requested
    if output_path.exists():
        if force:
            print(f"[INFO] Removing existing directory: {output_path}")
            if not dry_run:
                shutil.rmtree(output_path)
        else:
            print(f"[INFO] Output directory exists, appending new downloads: {output_path}")

    print(f"[INFO] Creating output directory: {output_path}")
    if not dry_run:
        output_path.mkdir(parents=True, exist_ok=True)

    for idx, taxon_id in enumerate(taxon_ids, start=1):
        print(f"\n[INFO] ({idx}/{len(taxon_ids)}) Downloading taxon {taxon_id} into {output_path}")

        # Download either limited or full
        if max_genomes:
            print(f"[INFO] Limiting download to {max_genomes} genomes for taxon {taxon_id}")
            accessions_file = output_path / f"accessions_{taxon_id}.txt"
            accessions_clean_file = output_path / f"accessions_{taxon_id}_clean.txt"

            summary_cmd = [
                "datasets", "summary", "genome", "taxon", str(taxon_id),
                "--limit", str(max_genomes),
                "--assembly-source", "RefSeq",
                "--report", "ids_only",
                "--as-json-lines"
            ]
            if reference:
                summary_cmd.append("--reference")

            if not dry_run:
                with open(accessions_file, 'w') as f:
                    subprocess.run(summary_cmd, stdout=f, check=True)
                dataformat_cmd = ["dataformat", "tsv", "genome", "--fields", "accession"]
                with open(accessions_file, 'r') as infile, open(accessions_clean_file, 'w') as outfile:
                    subprocess.run(dataformat_cmd, stdin=infile, stdout=outfile, check=True)
                with open(accessions_clean_file, 'r') as f:
                    lines = f.readlines()
                with open(accessions_clean_file, 'w') as f:
                    f.writelines(lines[1:])
            else:
                print(f"[CMD] {' '.join(summary_cmd)} | dataformat tsv genome --fields accession > {accessions_file}")
                print(f"[CMD] tail -n +2 {accessions_file} > {accessions_clean_file}")

            download_cmd = [
                "datasets", "download", "genome", "accession",
                "--inputfile", str(accessions_clean_file),
                "--dehydrated",
                "--include", "genome,gff3",
                "--filename", str(zip_path)
            ]
            run_command(download_cmd, dry_run)

            if not dry_run:
                accessions_file.unlink(missing_ok=True)
                accessions_clean_file.unlink(missing_ok=True)
        else:
            print("[INFO] Downloading all available genomes")
            download_cmd = [
                "datasets", "download", "genome", "taxon", str(taxon_id),
                "--dehydrated",
                "--include", "genome,gff3",
                "--filename", str(zip_path)
            ]
            run_command(download_cmd, dry_run)

        # Unzip and rehydrate per batch
        run_command(["unzip", "-o", str(zip_path), "-d", str(unzip_dir)], dry_run)
        if rehydrate:
            run_command(["datasets", "rehydrate", "--directory", str(unzip_dir)], dry_run)