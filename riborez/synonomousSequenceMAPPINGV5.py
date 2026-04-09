import os
import sys
import re
import pandas as pd

def extract_gcf_id(header):
    """Extract GCF/GCA accession from a FASTA header string."""
    if not isinstance(header, str):
        return None
    m = re.match(r'(GC[FA]_\d+\.\d+)', header)
    if m:
        return m.group(1)
    return None

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

    df = df.map(replace_rep)

    # Update sequence counts and restructure columns
    def update_counts_and_restructure(df, gene_name, folder_path, rep_to_mapped):
        # Get original sequence counts
        original_fasta = os.path.join(folder_path, f"{gene_name}.fasta")
        filtered_fasta = os.path.join(folder_path, "pmprimer_outputs", f"{gene_name}.filt.fasta")

        # Count sequences and unique input genomes from original FASTA
        original_count = 0
        input_genomes_count = 0
        if os.path.exists(original_fasta):
            input_gcf_ids = set()
            with open(original_fasta, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        original_count += 1
                        gcf = extract_gcf_id(line[1:].strip())
                        if gcf:
                            input_gcf_ids.add(gcf)
            input_genomes_count = len(input_gcf_ids)

        # Count sequences in filtered FASTA
        non_redundant_count = 0
        if os.path.exists(filtered_fasta):
            with open(filtered_fasta, 'r') as f:
                non_redundant_count = sum(1 for line in f if line.startswith('>'))

        # Count successful amplifications and bacteria metrics from amplicon CSV files
        successful_amplifications = []
        bacteria_amplified = []
        unique_bacteria = []

        for _, row in df.iterrows():
            primer_csv = row.get('PrimerPairCSV', '')
            if primer_csv:
                amplicon_csv = os.path.join(folder_path, "output", primer_csv)
                if os.path.exists(amplicon_csv):
                    try:
                        amplicon_df_pre = pd.read_csv(amplicon_csv)
                        if 'ErrorStatus' in amplicon_df_pre.columns:
                            amp_count = int((amplicon_df_pre['ErrorStatus'] == 'OK').sum())
                            amplicon_df = amplicon_df_pre[amplicon_df_pre['ErrorStatus'] == 'OK']
                        else:
                            amp_count = max(0, len(amplicon_df_pre))
                            amplicon_df = amplicon_df_pre
                        successful_amplifications.append(amp_count)

                        bact = 0
                        uniq = 0

                        try:
                            if 'Header' in amplicon_df.columns and 'AmpliconSequence' in amplicon_df.columns:
                                # Count unique genomes (GCF IDs) across all successfully amplified sequences
                                gcf_set = set()
                                for header in amplicon_df['Header']:
                                    gcf = extract_gcf_id(str(header) if isinstance(header, str) else '')
                                    if gcf:
                                        gcf_set.add(gcf)
                                bact = len(gcf_set)

                                # UniqueBacteria: genomes whose amplicon sequence is unique to that genome
                                seq_to_gcf = {}
                                for _, row_data in amplicon_df.iterrows():
                                    seq = row_data['AmpliconSequence']
                                    hdr = row_data['Header']
                                    if isinstance(seq, str) and isinstance(hdr, str):
                                        gcf = extract_gcf_id(hdr)
                                        if gcf:
                                            seq_to_gcf.setdefault(seq, set()).add(gcf)
                                uniq = sum(1 for gcf_s in seq_to_gcf.values() if len(gcf_s) == 1)

                        except Exception as e:
                            print(f"Warning: Could not parse bacteria metrics for {primer_csv}: {e}")

                        bacteria_amplified.append(bact)
                        unique_bacteria.append(uniq)

                    except Exception:
                        successful_amplifications.append(0)
                        bacteria_amplified.append(0)
                        unique_bacteria.append(0)
                else:
                    successful_amplifications.append(0)
                    bacteria_amplified.append(0)
                    unique_bacteria.append(0)
            else:
                successful_amplifications.append(0)
                bacteria_amplified.append(0)
                unique_bacteria.append(0)

        # Remove old columns if present
        columns_to_remove = ['NumInputSequences', 'NumberOfUniqueBacteria']
        for col in columns_to_remove:
            if col in df.columns:
                df = df.drop(columns=[col])

        # Insert new columns after PrimerPairCSV (index 0)
        df.insert(1, 'InputGenomes', input_genomes_count)
        df.insert(2, 'Original#ofSequences', original_count)
        df.insert(3, 'nonRedundantOriginal#ofSequences', non_redundant_count)
        df.insert(4, 'SequencesSuccessfullyAmplified', successful_amplifications)
        df.insert(5, 'BacteriaAmplified', bacteria_amplified)
        df.insert(6, 'UniqueBacteria', unique_bacteria)

        return df

    df = update_counts_and_restructure(df, gene_name, folder_path, rep_to_mapped)

    # --- Split into stats file and ASV mapping file ---

    # Identify ASV columns (ASV_i, ASV_i_Sequence, ASV_i_Headers)
    asv_col_pattern = re.compile(r'^ASV_\d+')
    asv_cols = [c for c in df.columns if asv_col_pattern.match(c)]
    stats_cols = [c for c in df.columns if not asv_col_pattern.match(c)]

    # 1. Stats file: one row per primer pair, no ASV data
    stats_df = df[stats_cols]
    stats_file = os.path.join(folder_path, f"{gene_name}_amplicon_stats.csv")
    stats_df.to_csv(stats_file, index=False)

    # 2. ASV mapping file: long format — one row per ASV per primer pair
    # Columns: PrimerPairCSV, ASV_Number, Sequence, GenomeHeaders
    asv_rows = []
    # Find max ASV index present
    asv_indices = sorted(set(
        int(re.match(r'^ASV_(\d+)', c).group(1))
        for c in asv_cols
    ))
    for _, row in df.iterrows():
        primer_csv = row['PrimerPairCSV']
        for i in asv_indices:
            seq_col = f"ASV_{i}_Sequence"
            hdr_col = f"ASV_{i}_Headers"
            if seq_col not in df.columns:
                continue
            seq = row.get(seq_col, '')
            headers = row.get(hdr_col, '')
            # Skip empty ASV slots (sparse rows with fewer ASVs than max)
            if not isinstance(seq, str) or seq.strip() == '':
                continue
            asv_rows.append({
                "PrimerPairCSV": primer_csv,
                "ASV_Number": i,
                "Sequence": seq,
                "GenomeHeaders": headers if isinstance(headers, str) else ''
            })

    asv_df = pd.DataFrame(asv_rows, columns=["PrimerPairCSV", "ASV_Number", "Sequence", "GenomeHeaders"])
    asv_file = os.path.join(folder_path, f"{gene_name}_asv_mapping.tsv")
    asv_df.to_csv(asv_file, sep='\t', index=False)

    print(f"Processed {folder_path}")
    print(f"  -> {stats_file} ({len(stats_df)} primer pairs)")
    print(f"  -> {asv_file} ({len(asv_df)} ASV entries)")

def main(parent_folder):
    for subdir in os.listdir(parent_folder):
        full_path = os.path.join(parent_folder, subdir)
        if os.path.isdir(full_path):
            process_subfolder(full_path)

# Main runner
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python synonomousSequenceMAPPINGV5.py <parent_folder>")
        sys.exit(1)

    parent_folder = sys.argv[1]
    main(parent_folder)
