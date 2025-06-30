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

    # Update NumInputSequences and NumUniqueASVs
    def update_counts(df, gene_name):
        num_input_sequences = []
        unique_bacteria_counts = []
        asv_cols = [col for col in df.columns if 'ASV_' in col and '_Headers' in col]

        for _, row in df.iterrows():
            total_hits = 0
            unique_hits = 0
            for col in asv_cols:
                val = row[col]
                if isinstance(val, str):
                    count = val.count(f'| {gene_name}')
                    total_hits += count
                    if count == 1:
                        unique_hits += 1
            num_input_sequences.append(total_hits)
            unique_bacteria_counts.append(unique_hits)

        df['NumInputSequences'] = num_input_sequences
        df['NumberOfUniqueBacteria'] = unique_bacteria_counts

        # Reorder columns to make NumberOfUniqueBacteria the 3rd column
        cols = df.columns.tolist()
        if 'NumberOfUniqueBacteria' in cols:
            cols.insert(2, cols.pop(cols.index('NumberOfUniqueBacteria')))
            df = df[cols]

        return df



    df = update_counts(df, gene_name)

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
