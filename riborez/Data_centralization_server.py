import os
import sys
import pandas as pd

def extract_best_row(parent_folder, output_dir=None):
    """
    Args:
        parent_folder: root folder containing per-gene subdirectories
        output_dir:    where to write outputs (summary CSV + best_asvs/).
                       Defaults to parent_folder when None (standalone use).
    """
    if output_dir is None:
        output_dir = parent_folder

    output_rows = []
    folder_basename = os.path.basename(os.path.normpath(parent_folder))
    output_csv_path = os.path.join(output_dir, f"{folder_basename}_best_amplicon_summary.csv")

    # Dedicated output folder for per-gene best-primer ASV tables
    best_asvs_dir = os.path.join(output_dir, "best_asvs")
    os.makedirs(best_asvs_dir, exist_ok=True)

    for subdir, _, files in os.walk(parent_folder):
        for file in files:
            if file.endswith("_amplicon_stats.csv"):
                file_path = os.path.join(subdir, file)
                try:
                    df = pd.read_csv(file_path, low_memory=False)
                    df_filtered = df[df["AmpliconLength"] < 500]
                    if not df_filtered.empty:
                        # Find the maximum number of unique ASVs
                        max_asvs = df_filtered["NumUniqueASVs"].max()
                        # Get all rows with the maximum ASV count
                        max_asv_rows = df_filtered[df_filtered["NumUniqueASVs"] == max_asvs]
                        
                        # If there's a tie, choose the one with highest MedianHammingDistance
                        if len(max_asv_rows) > 1:
                            best_row = max_asv_rows.loc[max_asv_rows["MedianHammingDistance"].idxmax()]
                        else:
                            best_row = max_asv_rows.iloc[0]
                        
                        # Check if new columns exist, otherwise use fallback values
                        has_new_columns = all(col in df.columns for col in ["Original#ofSequences", "nonRedundantOriginal#ofSequences", "SequencesSuccessfullyAmplified"])
                        has_input_genomes = "InputGenomes" in df.columns
                        has_bacteria_amplified = "BacteriaAmplified" in df.columns
                        has_unique_bacteria = "UniqueBacteria" in df.columns
                        has_old_columns = all(col in df.columns for col in ["NumInputSequences", "NumberOfUniqueBacteria"])
                        if has_new_columns:
                            output_row = {
                                "Filename": file,
                                "PrimerPairCSV": best_row.get("PrimerPairCSV", ""),
                                "InputGenomes": best_row.get("InputGenomes", 0) if has_input_genomes else 0,
                                "Original#ofSequences": best_row.get("Original#ofSequences", 0),
                                "nonRedundantOriginal#ofSequences": best_row.get("nonRedundantOriginal#ofSequences", 0),
                                "SequencesSuccessfullyAmplified": best_row.get("SequencesSuccessfullyAmplified", 0),
                                "BacteriaAmplified": best_row.get("BacteriaAmplified", 0) if has_bacteria_amplified else 0,
                                "UniqueBacteria": best_row.get("UniqueBacteria", 0) if has_unique_bacteria else 0,
                                "NumUniqueASVs": best_row.get("NumUniqueASVs", 0),
                                "MedianHammingDistance": best_row.get("MedianHammingDistance", 0),
                                "AmpliconLength": best_row.get("AmpliconLength", 0)
                            }
                        elif has_old_columns:
                            # Fallback for old format
                            output_row = {
                                "Filename": file,
                                "PrimerPairCSV": best_row.get("PrimerPairCSV", ""),
                                "InputGenomes": 0,
                                "Original#ofSequences": best_row.get("NumInputSequences", 0),
                                "nonRedundantOriginal#ofSequences": 0,
                                "SequencesSuccessfullyAmplified": 0,
                                "BacteriaAmplified": 0,
                                "UniqueBacteria": 0,
                                "NumUniqueASVs": best_row.get("NumUniqueASVs", 0),
                                "MedianHammingDistance": best_row.get("MedianHammingDistance", 0),
                                "AmpliconLength": best_row.get("AmpliconLength", 0)
                            }
                        else:
                            # Minimal fallback if neither format is available
                            output_row = {
                                "Filename": file,
                                "PrimerPairCSV": best_row.get("PrimerPairCSV", ""),
                                "InputGenomes": 0,
                                "Original#ofSequences": 0,
                                "nonRedundantOriginal#ofSequences": 0,
                                "SequencesSuccessfullyAmplified": 0,
                                "BacteriaAmplified": 0,
                                "UniqueBacteria": 0,
                                "NumUniqueASVs": best_row.get("NumUniqueASVs", 0),
                                "MedianHammingDistance": best_row.get("MedianHammingDistance", 0),
                                "AmpliconLength": best_row.get("AmpliconLength", 0)
                            }
                        output_rows.append(output_row)

                        # ── Export best-primer ASV mapping ────────────────────
                        gene_name = file.replace("_amplicon_stats.csv", "")
                        asv_map_path = os.path.join(subdir, f"{gene_name}_asv_mapping.tsv")
                        best_primer_csv = best_row.get("PrimerPairCSV", "")
                        if best_primer_csv and os.path.exists(asv_map_path):
                            try:
                                asv_df = pd.read_csv(asv_map_path, sep="\t")
                                best_asvs = asv_df[asv_df["PrimerPairCSV"] == best_primer_csv]
                                out_tsv = os.path.join(best_asvs_dir,
                                                        f"{gene_name}_best_asvs.tsv")
                                best_asvs.to_csv(out_tsv, sep="\t", index=False)
                                print(f"  → {gene_name}_best_asvs.tsv "
                                      f"({len(best_asvs)} ASVs, primer: {best_primer_csv})")
                            except Exception as e2:
                                print(f"  [WARN] Could not write {gene_name}_best_asvs.tsv: {e2}")

                except Exception as e:
                    print(f"Failed to process {file_path}: {e}")

    if output_rows:
        out_df = pd.DataFrame(output_rows)
        out_df.to_csv(output_csv_path, index=False)
        print(f"Output written to {output_csv_path}")
    else:
        print("No valid rows found in any subdirectories.")

def main(parent_folder, output_dir=None):
    extract_best_row(parent_folder, output_dir=output_dir)

# Main runner
if __name__ == "__main__":
    if len(sys.argv) not in (2, 3):
        print("Usage: python Data_centralization.py <parent_folder> [output_dir]")
        sys.exit(1)

    parent_folder = sys.argv[1]
    output_dir    = sys.argv[2] if len(sys.argv) == 3 else None
    main(parent_folder, output_dir=output_dir)
