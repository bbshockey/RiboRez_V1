#!/usr/bin/env python3
"""
Gene extraction functionality for RiboRez
"""

from pathlib import Path
from Bio import SeqIO
import os
import statistics
import random
from collections import defaultdict

def parse_attributes(attr_str):
    """Parse GFF/GTF attributes string into a dictionary."""
    attrs = {}
    for item in attr_str.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        
        # Handle GFF format: key=value
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key.strip()] = value.strip()
        # Handle GTF format: key "value"
        elif ' "' in item:
            parts = item.split(' "', 1)
            if len(parts) == 2:
                key = parts[0].strip()
                value = parts[1].rstrip('"').strip()
                attrs[key] = value
    
    return attrs

def extract_genes(taxon_name, data_root=None, output_dir=None, sample_size=None, random_seed=42, genes=None):
    """
    Extract genes from downloaded NCBI datasets.
    
    Args:
        taxon_name (str): Taxon name (used for directory structure)
        data_root (str): Path to the data directory (auto-detected if None)
        output_dir (str): Output directory for extracted genes (auto-generated if None)
        sample_size (int): Number of genomes to sample (None for all)
        random_seed (int): Random seed for sampling
        genes (list): List of specific genes to extract (None for all genes)
    """
    # Auto-detect data root if not provided
    if data_root is None:
        data_root = Path.cwd() / f"{taxon_name}_NCBI" / "genomes" / "ncbi_dataset" / "data"
    else:
        data_root = Path(data_root)
    
    if not data_root.exists():
        raise FileNotFoundError(f"Data directory not found: {data_root}")
    
    # Auto-generate output directory if not provided
    if output_dir is None:
        output_dir = Path.cwd() / f"{taxon_name}_AllGenesExtracted_rRNA"
    else:
        output_dir = Path(output_dir)
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Setup logging
    log_path = output_dir / "extraction_log.txt"
    log = open(log_path, "w")
    
    # Log gene extraction mode
    if genes:
        log.write(f"Extracting specific genes: {', '.join(genes)}\n")
        print(f"[INFO] Extracting specific genes: {', '.join(genes)}")
    else:
        log.write("Extracting all genes (CDS and rRNA)\n")
        print("[INFO] Extracting all genes (CDS and rRNA)")
    
    # Temporary directory for processing
    temp_dir = output_dir / "temp"
    temp_dir.mkdir(parents=True, exist_ok=True)
    
    # Get all genome directories
    all_dirs = [d for d in data_root.iterdir() if d.is_dir()]
    
    if not all_dirs:
        raise ValueError(f"No genome directories found in {data_root}")
    
    # Sample genomes if requested
    if sample_size and sample_size < len(all_dirs):
        random.seed(random_seed)
        sampled_dirs = random.sample(all_dirs, sample_size)
        print(f"[INFO] Sampling {sample_size} genomes from {len(all_dirs)} available")
    else:
        sampled_dirs = all_dirs
        print(f"[INFO] Processing all {len(all_dirs)} genomes")
    
    # Dictionary to store sequence lengths by gene (global across all genomes)
    gene_lengths = defaultdict(list)
    gene_headers = defaultdict(list)
    gene_sequences = defaultdict(list)
    
    # Process genomes and extract genes
    processed_count = 0
    for subdir in sampled_dirs:
        print(f"[INFO] Processing: {subdir.name}")
        log.write(f"\n[{subdir.name}]\n")
        
        # Look for FASTA files with various extensions
        fna_files = list(subdir.glob("*.fna"))
        fasta_files = list(subdir.glob("*.fasta"))
        fa_files = list(subdir.glob("*.fa"))
        
        # Combine all FASTA files found
        all_fasta_files = fna_files + fasta_files + fa_files
        
        gff_files = list(subdir.glob("*.gff"))
        gtf_files = list(subdir.glob("*.gtf"))
        
        # Debug output
        print(f"[DEBUG] {subdir.name}: Found {len(all_fasta_files)} FASTA files, {len(gff_files)} GFF files, {len(gtf_files)} GTF files")
        if all_fasta_files:
            print(f"[DEBUG] FASTA files: {[f.name for f in all_fasta_files]}")
        if gff_files:
            print(f"[DEBUG] GFF files: {[f.name for f in gff_files]}")
        if gtf_files:
            print(f"[DEBUG] GTF files: {[f.name for f in gtf_files]}")
        
        if not all_fasta_files:
            print(f"[WARNING] Missing FASTA files in {subdir.name}")
            continue
        
        # Use GFF if available, otherwise use GTF
        annotation_files = gff_files + gtf_files
        if not annotation_files:
            print(f"[WARNING] Missing GFF or GTF annotation files in {subdir.name}")
            continue
        
        fasta_path = all_fasta_files[0]
        gff_path = annotation_files[0]
        seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
        
        # Create chromosome name mapping for GFF files
        # GFF uses assembly IDs like "assembly_750c07d9cdfd47c5_1"
        # FASTA uses contig names like "ATCC_12022_contig_1" with assembly_id in description
        chrom_mapping = {}
        for seq_id, seq_record in seq_dict.items():
            # Extract assembly_id from description
            description = seq_record.description
            if 'assembly_id=' in description:
                assembly_id = description.split('assembly_id="')[1].split('"')[0]
                gff_chrom_name = f"assembly_{assembly_id}_1"
                # Store all contigs that map to this GFF chromosome
                if gff_chrom_name not in chrom_mapping:
                    chrom_mapping[gff_chrom_name] = []
                chrom_mapping[gff_chrom_name].append(seq_id)
                print(f"[DEBUG] {subdir.name}: Mapped GFF '{gff_chrom_name}' -> FASTA '{seq_id}'")
        
        # Debug: Show available chromosome names in FASTA
        print(f"[DEBUG] {subdir.name}: FASTA chromosomes: {list(seq_dict.keys())[:3]}...")  # Show first 3
        
        genes_found = 0
        with gff_path.open() as gff:
            for line in gff:
                if line.startswith("#"):
                    continue
                cols = line.strip().split("\t")
                if len(cols) < 9:
                    continue
                
                attr_dict = parse_attributes(cols[8])
                feature_type = cols[2]
                
                # Skip non-CDS and non-rRNA features
                if feature_type not in ["CDS", "rRNA"]:
                    continue
                
                genes_found += 1
                
                # Debug: Show first few attributes for troubleshooting
                if genes_found <= 3:
                    print(f"[DEBUG] {subdir.name} feature {genes_found}: type={feature_type}, attrs={dict(list(attr_dict.items())[:3])}")
                
                # Extract and normalize annotation fields
                product = attr_dict.get("product", "").lower()
                gene_field = (attr_dict.get("gene") or "").lower()
                locus_tag = (attr_dict.get("locus_tag") or "").lower()
                id_field = (attr_dict.get("ID") or "").lower()
                
                # Assign unified names to rRNAs
                if "16s" in product or "16s" in gene_field or "16s" in locus_tag or "16s" in id_field or "16S ribosomal RNA" in attr_dict.get("product", ""):
                    gene_name = "16S"
                elif "23s" in product or "23s" in gene_field or "23s" in locus_tag or "23s" in id_field or "23S ribosomal RNA" in attr_dict.get("product", ""):
                    gene_name = "23S"
                elif "5s" in product or "5s" in gene_field or "5s" in locus_tag or "5s" in id_field or "5S ribosomal RNA" in attr_dict.get("product", ""):
                    gene_name = "5S"
                else:
                    gene_name = (
                        attr_dict.get("gene") or
                        attr_dict.get("locus_tag") or
                        attr_dict.get("ID") or
                        f"{feature_type}_{cols[0]}_{cols[3]}_{cols[4]}"
                    )
                
                # Debug: Show gene naming for first few features
                if genes_found <= 3:
                    print(f"[DEBUG] {subdir.name} gene {genes_found}: name='{gene_name}', product='{product}', gene='{gene_field}', locus='{locus_tag}'")
                
                # Filter by specific genes if requested
                if genes:
                    # Check if this gene matches any of the requested genes
                    gene_matches = False
                    for requested_gene in genes:
                        # Handle rRNA requests
                        if requested_gene.lower() in ["rrna", "rna"]:
                            if gene_name in ["16S", "23S", "5S"]:
                                gene_matches = True
                                break
                        # Handle specific rRNA requests
                        elif requested_gene.upper() in ["16S", "23S", "5S"]:
                            if gene_name == requested_gene.upper():
                                gene_matches = True
                                break
                        # Handle gene name matching (case-insensitive)
                        elif requested_gene.lower() in gene_name.lower():
                            gene_matches = True
                            break
                    
                    if not gene_matches:
                        continue
                
                chrom = cols[0]
                start, end = int(cols[3]), int(cols[4])
                strand = cols[6]
                
                # Try to find the sequence using chromosome mapping
                seq_record = seq_dict.get(chrom)
                if not seq_record and chrom in chrom_mapping:
                    # Try each mapped contig to find the one containing this feature
                    mapped_contigs = chrom_mapping[chrom]
                    for mapped_chrom in mapped_contigs:
                        potential_seq = seq_dict.get(mapped_chrom)
                        if potential_seq and len(potential_seq.seq) >= end:
                            seq_record = potential_seq
                            if genes_found <= 3:
                                print(f"[DEBUG] {subdir.name} gene {genes_found}: Mapped GFF '{chrom}' -> FASTA '{mapped_chrom}' (length: {len(potential_seq.seq)})")
                            break
                
                if not seq_record:
                    if genes_found <= 3:
                        print(f"[DEBUG] {subdir.name} gene {genes_found}: Chromosome '{chrom}' not found in FASTA")
                    continue
                
                sub_seq = seq_record.seq[start - 1:end]
                if strand == "-":
                    sub_seq = sub_seq.reverse_complement()
                
                if genes_found <= 3:
                    print(f"[DEBUG] {subdir.name} gene {genes_found}: Extracted sequence length {len(sub_seq)} from {chrom}:{start}-{end}")
                
                seq_length = len(sub_seq)
                gene_lengths[gene_name].append(seq_length)
                
                strain_info = seq_record.description
                full_header = f"{subdir.name}|{seq_record.id}|{start}-{end}|{strain_info}"
                gene_headers[gene_name].append(full_header)
                gene_sequences[gene_name].append(str(sub_seq))
                
                log.write(f"  {gene_name}: {full_header} (length: {seq_length})\n")
        
        print(f"[DEBUG] {subdir.name}: Found {genes_found} CDS/rRNA features, extracted {len(gene_headers)} unique genes")
        print(f"[DEBUG] {subdir.name}: Gene names found: {list(gene_headers.keys())[:5]}...")  # Show first 5 gene names
        print(f"[DEBUG] Total genes across all genomes so far: {sum(len(headers) for headers in gene_headers.values())}")
        processed_count += 1
    
    # Calculate median lengths and apply length filter
    print("[INFO] Applying length filter...")
    log.write("\n\n--- Length Filtering ---\n")
    
    filtered_count = 0
    total_count = 0
    small_gene_groups = 0
    saved_gene_groups = 0
    
    for gene_name, lengths in gene_lengths.items():
        if len(lengths) < 1:
            continue
        
        total_count += len(lengths)
        median_length = statistics.median(lengths)
        threshold = median_length / 2
        
        # Special length filter for 16S sequences
        if gene_name == "16S":
            min_16s_length = 1000  # Lowered from 1400 to 1000 bp for 16S
            log.write(f"{gene_name}: median={median_length}, threshold={threshold}, 16S_min_length={min_16s_length}\n")
        else:
            log.write(f"{gene_name}: median={median_length}, threshold={threshold}\n")
        
        filtered_sequences = []
        for i in range(len(lengths)):
            # Apply length filtering
            passes_length_filter = True
            
            # Special 16S length filter
            if gene_name == "16S" and lengths[i] < min_16s_length:
                passes_length_filter = False
                log.write(f"  Filtered 16S: {gene_headers[gene_name][i]} (length: {lengths[i]} < {min_16s_length} bp)\n")
            # Standard median-based filter
            elif lengths[i] < threshold:
                passes_length_filter = False
                log.write(f"  Filtered: {gene_headers[gene_name][i]} (length: {lengths[i]} < {threshold})\n")
            
            if passes_length_filter:
                filtered_sequences.append((gene_headers[gene_name][i], gene_sequences[gene_name][i]))
            else:
                filtered_count += 1
        
        if len(filtered_sequences) >= 5:
            saved_gene_groups += 1
            with open(output_dir / f"{gene_name}.fasta", "w") as f:
                for header, sequence in filtered_sequences:
                    f.write(f">{header} | {gene_name}\n{sequence}\n")
        else:
            small_gene_groups += 1
            log.write(f"  Skipped creating file for {gene_name}: only {len(filtered_sequences)} sequences (< 5 required)\n")
    
    # Close log
    log.close()
    
    # Print summary
    print(f"[SUCCESS] Gene extraction completed!")
    print(f"[INFO] Processed {processed_count} genomes")
    print(f"[INFO] Total sequences: {total_count}")
    print(f"[INFO] Filtered out {filtered_count} sequences (< 50% of median length)")
    print(f"[INFO] Skipped {small_gene_groups} gene groups with fewer than 5 sequences")
    print(f"[INFO] Created {saved_gene_groups} gene FASTA files with 5+ sequences")
    print(f"[INFO] Output directory: {output_dir}")
    print(f"[INFO] Log file: {log_path}")
    
    return output_dir 