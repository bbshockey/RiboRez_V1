import os
import sys
import pandas as pd

def extract_best_row(parent_folder):
    output_rows = []
    folder_basename = os.path.basename(os.path.normpath(parent_folder))
    output_csv_path = os.path.join(parent_folder, f"{folder_basename}_best_amplicon_summary.csv")

    for subdir, _, files in os.walk(parent_folder):
        for file in files:
            if file.endswith("processed_amplicon.csv"):
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
                        has_new_columns = all(col in df.columns for col in ["Original#ofSequences", "nonRedundantOriginal#ofSequences", "SuccessfulAmplifications"])
                        has_old_columns = all(col in df.columns for col in ["NumInputSequences", "NumberOfUniqueBacteria"])
                        if has_new_columns:
                            output_row = {
                                "Filename": file,
                                "PrimerPairCSV": best_row.get("PrimerPairCSV", ""),
                                "Original#ofSequences": best_row.get("Original#ofSequences", 0),
                                "nonRedundantOriginal#ofSequences": best_row.get("nonRedundantOriginal#ofSequences", 0),
                                "SuccessfulAmplifications": best_row.get("SuccessfulAmplifications", 0),
                                "NumUniqueASVs": best_row.get("NumUniqueASVs", 0),
                                "MedianHammingDistance": best_row.get("MedianHammingDistance", 0),
                                "AmpliconLength": best_row.get("AmpliconLength", 0)
                            }
                        elif has_old_columns:
                            # Fallback for old format
                            output_row = {
                                "Filename": file,
                                "PrimerPairCSV": best_row.get("PrimerPairCSV", ""),
                                "Original#ofSequences": best_row.get("NumInputSequences", 0),  # Use old column as fallback
                                "nonRedundantOriginal#ofSequences": 0,  # Not available in old format
                                "SuccessfulAmplifications": 0,  # Not available in old format
                                "NumUniqueASVs": best_row.get("NumUniqueASVs", 0),
                                "MedianHammingDistance": best_row.get("MedianHammingDistance", 0),
                                "AmpliconLength": best_row.get("AmpliconLength", 0)
                            }
                        else:
                            # Minimal fallback if neither format is available
                            output_row = {
                                "Filename": file,
                                "PrimerPairCSV": best_row.get("PrimerPairCSV", ""),
                                "Original#ofSequences": 0,
                                "nonRedundantOriginal#ofSequences": 0,
                                "SuccessfulAmplifications": 0,
                                "NumUniqueASVs": best_row.get("NumUniqueASVs", 0),
                                "MedianHammingDistance": best_row.get("MedianHammingDistance", 0),
                                "AmpliconLength": best_row.get("AmpliconLength", 0)
                            }
                        output_rows.append(output_row)
                except Exception as e:
                    print(f"Failed to process {file_path}: {e}")

    if output_rows:
        out_df = pd.DataFrame(output_rows)
        out_df.to_csv(output_csv_path, index=False)
        print(f"Output written to {output_csv_path}")
    else:
        print("No valid rows found in any subdirectories.")

def main(parent_folder):
    extract_best_row(parent_folder)

# Main runner
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python Data_centralization.py <parent_folder>")
        sys.exit(1)

    parent_folder = sys.argv[1]
    main(parent_folder)
