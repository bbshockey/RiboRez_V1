#!/usr/bin/env python3
"""
Primer design functionality for RiboRez using PMPrimer
"""

import os
import shutil
import subprocess
import datetime
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO

def log_message(message, log_file):
    """Log a message with timestamp."""
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(log_file, "a") as log:
        log.write(f"[{timestamp}] {message}\n")
    print(f"[{timestamp}] {message}")

def count_sequences(fasta_path):
    """Count sequences in a FASTA file."""
    try:
        return sum(1 for line in open(fasta_path) if line.startswith(">"))
    except Exception as e:
        print(f"Error counting sequences in {fasta_path}: {str(e)}")
        return 0

def create_reference_mapping(gene_dir, gene_name, reference_mapping_dir):
    """Create reference mapping for a gene."""
    try:
        original_fasta = os.path.join(gene_dir, f"{gene_name}.fasta")
        filtered_fasta = os.path.join(gene_dir, f"{gene_name}.filt.fasta")
        
        if not os.path.exists(filtered_fasta):
            print(f"Filtered FASTA not found for {gene_name}, skipping reference mapping")
            return None
            
        rep_dict = {str(rec.seq).upper(): rec.id for rec in SeqIO.parse(filtered_fasta, "fasta")}
        mapping = defaultdict(list)
        for rec in SeqIO.parse(original_fasta, "fasta"):
            seq = str(rec.seq).upper()
            if seq in rep_dict:
                mapping[rep_dict[seq]].append(rec.id)

        mapping_file = os.path.join(gene_dir, f"{gene_name}_reference_mapping.tsv")
        with open(mapping_file, "w") as out_f:
            out_f.write("Representative\tTotal_Mapped\tRedundant_Count\tMapped_Headers\n")
            for rep, headers in mapping.items():
                out_f.write(f"{rep}\t{len(headers)}\t{len(headers) - 1}\t{';'.join(headers)}\n")

        shutil.copy(mapping_file, os.path.join(reference_mapping_dir, os.path.basename(mapping_file)))
        return mapping_file
    except Exception as e:
        print(f"Error creating reference mapping for {gene_name}: {str(e)}")
        return None

def run_pmprimer_in_subdir(fasta_path, output_folder, reference_mapping_dir, min_sequences, log_file):
    """Run PMPrimer on a single FASTA file."""
    filename = os.path.basename(fasta_path)
    gene_name = os.path.splitext(filename)[0]
    seq_count = count_sequences(fasta_path)
    
    if seq_count < min_sequences:
        log_message(f"Skipping {filename}: Only {seq_count} sequences (minimum: {min_sequences})", log_file)
        return f"{filename}: SKIPPED"

    gene_dir = os.path.join(output_folder, gene_name)
    os.makedirs(gene_dir, exist_ok=True)
    copied_fasta_path = os.path.join(gene_dir, filename)
    
    if not os.path.exists(copied_fasta_path):
        shutil.copy(fasta_path, copied_fasta_path)

    mc_filt_path = os.path.join(gene_dir, f"{gene_name}.filt.mc.fasta")
    
    # Run alignment
    alignment_cmd = ["pmprimer", "-f", copied_fasta_path, "-p", "default", "-a", "muscle", "-e", "save"]
    log_message(f"Running alignment for {filename}: {' '.join(alignment_cmd)}", log_file)
    
    try:
        subprocess.run(alignment_cmd, check=True, cwd=gene_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        log_message(f"Alignment failed for {filename}: {e.stderr}", log_file)

    if not os.path.exists(mc_filt_path):
        log_message(f"Alignment failed or not produced for {filename}; falling back to raw fasta", log_file)
        mc_filt_path = copied_fasta_path

    # Primary primer design
    primary_cmd = ["pmprimer", "-f", mc_filt_path, "-p", "default", "-a", "threshold:0.5", "gaps:30", "tm:35", "primer2", "-e", "hpcnt:1000", "maxlen:1000", "save", "-d", "2"]
    fallback_cmd = ["pmprimer", "-f", mc_filt_path, "-a", "primer2", "-e", "save"]
    
    try:
        subprocess.run(primary_cmd, check=True, cwd=gene_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        log_message(f"Primary primer design successful for {filename}", log_file)
    except subprocess.CalledProcessError as e:
        log_message(f"Primary primer design failed for {filename}: {e.stderr}", log_file)

    try:
        subprocess.run(fallback_cmd, check=True, cwd=gene_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        log_message(f"Fallback primer design successful for {filename}", log_file)
    except subprocess.CalledProcessError as e:
        log_message(f"Fallback primer design failed for {filename}: {e.stderr}", log_file)

    # Additional primer design commands
    additional_cmds = [
        ["pmprimer", "-f", mc_filt_path, "-p", "default", "-a", "threshold:0.3", "gaps:90", "tm:35", "minlen:15", "maxlen:800", "merge", "primer2", "-e", "hpcnt:1000", "minlen:50", "save", "-d", "2"],
        ["pmprimer", "-f", mc_filt_path, "-p", "default", "-a", "threshold:0.3", "gaps:50", "tm:30", "primer2", "-e", "hpcnt:2000", "maxlen:2000", "save", "-d", "2"],
        ["pmprimer", "-f", mc_filt_path, "-a", "threshold:0.2", "merge", "primer2", "-e", "hpcnt:3000", "save"],
        ["pmprimer", "-f", mc_filt_path, "-a", "threshold:0.3", "minlen:5", "merge", "primer2", "-e", "save"]
    ]
    
    for cmd in additional_cmds:
        try:
            subprocess.run(cmd, check=True, cwd=gene_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            log_message(f"Additional command successful: {' '.join(cmd)}", log_file)
        except subprocess.CalledProcessError as e:
            log_message(f"Additional command failed: {' '.join(cmd)}; {e.stderr}", log_file)

    create_reference_mapping(gene_dir, gene_name, reference_mapping_dir)
    return f"{filename}: COMPLETED"

def check_pmprimer_installed():
    """Check if PMPrimer is installed and available."""
    try:
        result = subprocess.run(["pmprimer", "--help"], 
                              capture_output=True, text=True)
        return result.returncode == 0
    except FileNotFoundError:
        return False

def design_primers(input_folder, output_folder=None, min_sequences=10, threads=8):
    """
    Design primers for genes in the input folder using PMPrimer.
    
    Args:
        input_folder (str): Path to folder containing FASTA files
        output_folder (str): Output directory (auto-generated if None)
        min_sequences (int): Minimum number of sequences required
        threads (int): Number of threads to use
    """
    # Check if PMPrimer is installed
    if not check_pmprimer_installed():
        raise RuntimeError("PMPrimer is not installed. Please install PMPrimer first.")
    
    # Setup paths
    input_folder = Path(input_folder)
    if not input_folder.exists():
        raise FileNotFoundError(f"Input folder not found: {input_folder}")
    
    if output_folder is None:
        output_folder = input_folder.parent / f"{input_folder.name}_Primers"
    else:
        output_folder = Path(output_folder)
    
    output_folder.mkdir(parents=True, exist_ok=True)
    reference_mapping_dir = output_folder / "reference_mappings"
    reference_mapping_dir.mkdir(exist_ok=True)
    
    log_file = output_folder / "pmprimer_log.log"
    
    # Initialize log
    with open(log_file, "w") as log:
        log.write(f"=== PMPrimer Processing Log Started at {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} ===\n\n")
    
    print(f"[INFO] Starting primer design for {input_folder}")
    print(f"[INFO] Output directory: {output_folder}")
    print(f"[INFO] Minimum sequences: {min_sequences}")
    print(f"[INFO] Threads: {threads}")
    
    # Find FASTA files
    fasta_files = [f for f in input_folder.iterdir() if f.suffix == ".fasta"]
    log_message(f"Found {len(fasta_files)} FASTA files to process", log_file)
    
    if not fasta_files:
        raise ValueError(f"No FASTA files found in {input_folder}")
    
    # Filter files by sequence count
    filtered_files = []
    skipped_count = 0
    for fasta_file in fasta_files:
        if count_sequences(fasta_file) >= min_sequences:
            filtered_files.append(fasta_file)
        else:
            skipped_count += 1
            if skipped_count % 100 == 0:
                log_message(f"Skipped {skipped_count} files so far, latest: {fasta_file.name}", log_file)
    
    log_message(f"Pre-filtering complete: {skipped_count} skipped, {len(filtered_files)} to process", log_file)
    
    if not filtered_files:
        log_message("No files meet the minimum sequence requirement", log_file)
        return output_folder
    
    # Process files in parallel
    successful = 0
    failed = 0
    
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {
            executor.submit(run_pmprimer_in_subdir, f, output_folder, reference_mapping_dir, min_sequences, log_file): f 
            for f in filtered_files
        }
        
        for future in as_completed(futures):
            try:
                result = future.result()
                if "COMPLETED" in result:
                    successful += 1
                else:
                    failed += 1
            except Exception as e:
                failed += 1
                log_message(f"Unhandled error: {str(e)}", log_file)
    
    # Post-processing reference mappings
    log_message("=== Post-processing reference mappings ===", log_file)
    count = 0
    for d in output_folder.iterdir():
        if d.is_dir() and d.name != "reference_mappings":
            gene_dir = d
            gene_name = d.name
            mapping_file = gene_dir / f"{gene_name}_reference_mapping.tsv"
            if not mapping_file.exists():
                if create_reference_mapping(gene_dir, gene_name, reference_mapping_dir):
                    count += 1
    
    log_message(f"Created {count} additional reference mappings", log_file)
    log_message(f"Processing complete: {successful} successful, {failed} failed", log_file)
    
    print(f"[SUCCESS] Primer design completed!")
    print(f"[INFO] Successful: {successful}, Failed: {failed}")
    print(f"[INFO] Output directory: {output_folder}")
    print(f"[INFO] Log file: {log_file}")
    
    return output_folder 