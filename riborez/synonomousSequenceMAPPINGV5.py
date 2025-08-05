import os
import sys
import pandas as pd

def process_subfolder(folder_path):
    # Find mapping file
    mapping_file = [f for f in os.listdir(folder_path) if f.endswith("_reference_mapping.tsv")]
    if not mapping_file:
        print(f"Skipping {folder_path}: No mapping file found.")
        return
    mapping_file = os.path.join(folder_path, mapping_file[0])
    gene_name = os.path.basename(mapping_file).split('_')[0]
    mapping_df = pd.read_csv(mapping_file, sep='\t')
    rep_to_mapped = dict(zip(mapping_df['Representative'], mapping_df['Mapped_Headers']))

    # Find amplicon file
    amplicon_file = [f for f in os.listdir(folder_path) if f.endswith("amplicon.summary.csv")]
    if not amplicon_file:
        print(f"Skipping {folder_path}: No amplicon summary file found.")
        return
    amplicon_file = os.path.join(folder_path, amplicon_file[0])
    df = pd.read_csv(amplicon_file)

    # Replace representatives with mapped headers + gene name
    def replace_rep(cell):
        if isinstance(cell, str):
            for rep, mapped in rep_to_mapped.items():
                # Skip if rep is identical to mapped (no real mapping)
                if rep == mapped:
                    continue
                if rep in cell:
                    mapped_entries = mapped.split(';')
                    formatted = ', '.join([entry + f" | {gene_name}" for entry in mapped_entries])
                    cell = cell.replace(rep, formatted)
        return cell


    df = df.applymap(replace_rep)

    # Update sequence counts and restructure columns
    def update_counts_and_restructure(df, gene_name, folder_path):
        # Get original sequence counts
        original_fasta = os.path.join(folder_path, f"{gene_name}.fasta")
        filtered_fasta = os.path.join(folder_path, f"{gene_name}.filt.fasta")
        
        # Count sequences in original FASTA
        original_count = 0
        if os.path.exists(original_fasta):
            with open(original_fasta, 'r') as f:
                original_count = sum(1 for line in f if line.startswith('>'))
        
        # Count sequences in filtered FASTA
        non_redundant_count = 0
        if os.path.exists(filtered_fasta):
            with open(filtered_fasta, 'r') as f:
                non_redundant_count = sum(1 for line in f if line.startswith('>'))
        
        # Count successful amplifications and unique bacteria from amplicon CSV files
        successful_amplifications = []
        unique_bacteria_counts = []
        for _, row in df.iterrows():
            primer_csv = row.get('PrimerPairCSV', '')
            if primer_csv:
                # Look for the amplicon CSV file in the output directory
                amplicon_csv = os.path.join(folder_path, "output", primer_csv)
                if os.path.exists(amplicon_csv):
                    try:
                        # Count rows in the CSV file (excluding header)
                        with open(amplicon_csv, 'r') as f:
                            # Count lines and subtract 1 for header
                            line_count = sum(1 for line in f)
                            amp_count = max(0, line_count - 1)  # Subtract header row
                        successful_amplifications.append(amp_count)
                        
                        # Count unique bacteria by parsing headers
                        unique_bacteria = 0
                        try:
                            amplicon_df = pd.read_csv(amplicon_csv)
                            if 'Header' in amplicon_df.columns:
                                # Extract GCF identifiers from headers
                                gcf_identifiers = set()
                                for header in amplicon_df['Header']:
                                    if isinstance(header, str):
                                        # Extract GCF identifier (e.g., GCF_000005845.2_EColi -> GCF_000005845.2)
                                        if 'GCF_' in header:
                                            gcf_part = header.split('|')[0]  # Get first part before |
                                            if '_' in gcf_part:
                                                gcf_id = gcf_part.split('_')[:-1]  # Remove strain name
                                                gcf_id = '_'.join(gcf_id)
                                                gcf_identifiers.add(gcf_id)
                                unique_bacteria = len(gcf_identifiers)
                        except Exception as e:
                            print(f"Warning: Could not parse unique bacteria for {primer_csv}: {e}")
                            unique_bacteria = 0
                        
                        unique_bacteria_counts.append(unique_bacteria)
                    except:
                        successful_amplifications.append(0)
                        unique_bacteria_counts.append(0)
                else:
                    successful_amplifications.append(0)
                    unique_bacteria_counts.append(0)
            else:
                successful_amplifications.append(0)
                unique_bacteria_counts.append(0)
        
        # Remove problematic columns
        columns_to_remove = ['NumInputSequences', 'NumberOfUniqueBacteria']
        for col in columns_to_remove:
            if col in df.columns:
                df = df.drop(columns=[col])
        
        # Add new columns at the beginning (after PrimerPairCSV)
        df.insert(1, 'Original#ofSequences', original_count)
        df.insert(2, 'nonRedundantOriginal#ofSequences', non_redundant_count)
        df.insert(3, 'SuccessfulAmplifications', successful_amplifications)
        df.insert(4, 'UniqueBacteria', unique_bacteria_counts)
        
        return df



    df = update_counts_and_restructure(df, gene_name, folder_path)

    # Save updated file
    output_file = os.path.join(folder_path, f"{gene_name}_processed_amplicon.csv")
    df.to_csv(output_file, index=False)
    print(f"Processed {folder_path} -> {output_file}")

def main(parent_folder):
    for subdir in os.listdir(parent_folder):
        full_path = os.path.join(parent_folder, subdir)
        if os.path.isdir(full_path):
            process_subfolder(full_path)

# Main runner
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python synonomousSequenceMAPPINGV4.py <parent_folder>")
        sys.exit(1)

    parent_folder = sys.argv[1]
    main(parent_folder)
