#!/usr/bin/env python3
"""
Amplicon analysis functionality for RiboRez
Chains multiple analysis scripts together for comprehensive primer evaluation
"""

import os
import subprocess
import pandas as pd
import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

def log_message(message, log_file):
    """Log a message with timestamp."""
    from datetime import datetime
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(log_file, "a") as log:
        log.write(f"[{timestamp}] {message}\n")
    print(f"[{timestamp}] {message}")

def run_amplification_analysis(input_dir, log_file):
    """
    Run the amplification and analysis workflow for a single gene directory.
    This corresponds to Post_alignjson_many_working_Server.py functionality.
    """
    try:
        # Find the best JSON file (most amplicons ≤ 500 bp)
        all_files = os.listdir(input_dir)
        candidate_jsons = []
        
        for f in all_files:
            if f.endswith(".json"):
                stem = f.replace(".json", "")
                csv_file = stem + ".csv"
                csv_path = os.path.join(input_dir, csv_file)
                if os.path.exists(csv_path):
                    candidate_jsons.append((f, csv_file))

        if not candidate_jsons:
            log_message(f"No valid JSON/CSV pairs found in {input_dir}", log_file)
            return False

        # Find the best JSON based on amplicon count
        best_json = None
        max_amplicons = -1

        for json_file, csv_file in candidate_jsons:
            try:
                csv_path = os.path.join(input_dir, csv_file)
                df = pd.read_csv(csv_path)
                df.columns = [col.strip() for col in df.columns]

                length_col = [
                    col for col in df.columns
                    if "length" in col.lower() and ("amplicon" in col.lower() or "effective" in col.lower())
                ]
                if not length_col:
                    continue
                col = length_col[0]
                count = df[col].apply(pd.to_numeric, errors='coerce').dropna().le(500).sum()
                if count > max_amplicons:
                    max_amplicons = count
                    best_json = os.path.join(input_dir, json_file)
            except Exception as e:
                log_message(f"Failed to read {csv_file}: {e}", log_file)

        # Find FASTA file
        fasta_file = next((os.path.join(input_dir, f) for f in all_files if f.endswith(".filt.mc.fasta")), None)
        
        if not best_json or not fasta_file:
            log_message(f"Skipping {os.path.basename(input_dir)}: no valid json or fasta found", log_file)
            return False

        log_message(f"Selected {os.path.basename(best_json)} ({max_amplicons} good amplicons) for {os.path.basename(input_dir)}", log_file)

        # Create output directory
        output_folder = os.path.join(input_dir, "output")
        os.makedirs(output_folder, exist_ok=True)

        # Run in silico amplification
        from .InSilico_amplification_hammingfixed_Server import main as amplify_main
        amplify_main(best_json, fasta_file, output_folder)
        
        # Run amplicon analysis
        summary_csv = os.path.join(input_dir, "amplicon.summary.csv")
        from .analyze_ampliconV2_Server import main as analyze_main
        analyze_main(output_folder, summary_csv)
        
        return True
        
    except Exception as e:
        log_message(f"Error processing {input_dir}: {e}", log_file)
        return False

def run_sequence_mapping(input_dir, log_file):
    """
    Run sequence mapping for a single gene directory.
    This corresponds to synonomousSequenceMAPPINGV5.py functionality.
    """
    try:
        from .synonomousSequenceMAPPINGV5 import process_subfolder
        process_subfolder(input_dir)
        return True
    except Exception as e:
        log_message(f"Error in sequence mapping for {input_dir}: {e}", log_file)
        return False

def run_data_centralization(parent_folder, log_file):
    """
    Run data centralization for the entire project.
    This corresponds to Data_centralization_server.py functionality.
    """
    try:
        from .Data_centralization_server import extract_best_row
        extract_best_row(parent_folder)
        return True
    except Exception as e:
        log_message(f"Error in data centralization: {e}", log_file)
        return False

def analyze_amplicons(input_folder, output_folder=None, threads=8):
    """
    Run comprehensive amplicon analysis workflow.
    
    Args:
        input_folder (str): Path to primer design output folder
        output_folder (str): Output directory (auto-generated if None)
        threads (int): Number of threads to use
    """
    # Setup paths
    input_folder = Path(input_folder)
    if not input_folder.exists():
        raise FileNotFoundError(f"Input folder not found: {input_folder}")
    
    if output_folder is None:
        output_folder = input_folder.parent / f"{input_folder.name}_AmpliconAnalysis"
    else:
        output_folder = Path(output_folder)
    
    output_folder.mkdir(parents=True, exist_ok=True)
    log_file = output_folder / "amplicon_analysis_log.log"
    
    # Initialize log
    with open(log_file, "w") as log:
        log.write(f"=== Amplicon Analysis Started at {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')} ===\n\n")
    
    print(f"[INFO] Starting amplicon analysis for {input_folder}")
    print(f"[INFO] Output directory: {output_folder}")
    print(f"[INFO] Threads: {threads}")
    
    # Find gene directories (skip reference_mappings)
    gene_dirs = [d for d in input_folder.iterdir() if d.is_dir() and d.name != "reference_mappings"]
    
    if not gene_dirs:
        raise ValueError(f"No gene directories found in {input_folder}")
    
    log_message(f"Found {len(gene_dirs)} gene directories to process", log_file)
    
    # Step 1: Run amplification and analysis for each gene
    print("[INFO] Step 1: Running amplification and analysis...")
    successful_amplification = 0
    failed_amplification = 0
    
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {
            executor.submit(run_amplification_analysis, str(gene_dir), log_file): gene_dir 
            for gene_dir in gene_dirs
        }
        
        for future in as_completed(futures):
            gene_dir = futures[future]
            try:
                if future.result():
                    successful_amplification += 1
                else:
                    failed_amplification += 1
            except Exception as e:
                failed_amplification += 1
                log_message(f"Unhandled error in amplification for {gene_dir.name}: {str(e)}", log_file)
    
    log_message(f"Amplification analysis complete: {successful_amplification} successful, {failed_amplification} failed", log_file)
    
    # Step 2: Run sequence mapping for each gene
    print("[INFO] Step 2: Running sequence mapping...")
    successful_mapping = 0
    failed_mapping = 0
    
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {
            executor.submit(run_sequence_mapping, str(gene_dir), log_file): gene_dir 
            for gene_dir in gene_dirs
        }
        
        for future in as_completed(futures):
            gene_dir = futures[future]
            try:
                if future.result():
                    successful_mapping += 1
                else:
                    failed_mapping += 1
            except Exception as e:
                failed_mapping += 1
                log_message(f"Unhandled error in mapping for {gene_dir.name}: {str(e)}", log_file)
    
    log_message(f"Sequence mapping complete: {successful_mapping} successful, {failed_mapping} failed", log_file)
    
    # Step 3: Run data centralization
    print("[INFO] Step 3: Running data centralization...")
    if run_data_centralization(str(input_folder), log_file):
        log_message("Data centralization completed successfully", log_file)
    else:
        log_message("Data centralization failed", log_file)
    
    # Copy results to output folder
    try:
        # Copy the centralized summary
        centralized_file = input_folder / f"{input_folder.name}_best_amplicon_summary.csv"
        if centralized_file.exists():
            import shutil
            shutil.copy2(centralized_file, output_folder / centralized_file.name)
            log_message(f"Copied centralized summary to {output_folder}", log_file)
    except Exception as e:
        log_message(f"Error copying results: {e}", log_file)
    
    # Print summary
    print(f"[SUCCESS] Amplicon analysis completed!")
    print(f"[INFO] Amplification analysis: {successful_amplification} successful, {failed_amplification} failed")
    print(f"[INFO] Sequence mapping: {successful_mapping} successful, {failed_mapping} failed")
    print(f"[INFO] Output directory: {output_folder}")
    print(f"[INFO] Log file: {log_file}")
    
    return output_folder 