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

def download_taxa(taxon_name, taxon_id, output_dir, rehydrate, force, dry_run, max_genomes=None):
    """
    Download and unpack NCBI genomes for a taxon.
    
    Args:
        taxon_name (str): Taxon name (used for output folder name)
        taxon_id (int): NCBI Taxon ID
        output_dir (str): Optional custom output directory
        rehydrate (bool): Whether to rehydrate datasets
        force (bool): Overwrite output directory if it exists
        dry_run (bool): Print commands without executing
        max_genomes (int): Maximum number of genomes to download (optional)
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

    # Step 1: Download dehydrated genome data
    download_cmd = [
        "datasets", "download", "genome", "taxon", str(taxon_id),
        "--dehydrated",
        "--include", "genome,gff3",
        "--filename", str(zip_path)
    ]
    
    # Add max-genomes limit if specified
    if max_genomes is not None:
        download_cmd.extend(["--max-genomes", str(max_genomes)])
        print(f"[INFO] Limiting download to {max_genomes} genomes")
    
    run_command(download_cmd, dry_run)

    # Step 2: Unzip
    unzip_cmd = ["unzip", "-o", str(zip_path), "-d", str(unzip_dir)]
    run_command(unzip_cmd, dry_run)

    # Step 3: Rehydrate (optional)
    if rehydrate:
        rehydrate_cmd = ["datasets", "rehydrate", "--directory", str(unzip_dir)]
        run_command(rehydrate_cmd, dry_run) 