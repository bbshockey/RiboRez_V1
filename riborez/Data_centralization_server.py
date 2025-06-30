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
                        best_row = df_filtered.loc[df_filtered["NumUniqueASVs"].idxmax()]
                        output_row = {
                            "Filename": file,
                            "PrimerPairCSV": best_row.get("PrimerPairCSV", ""),
                            "NumInputSequences": best_row.get("NumInputSequences", 0),
                            "NumberOfUniqueBacteria": best_row.get("NumberOfUniqueBacteria", 0),
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
