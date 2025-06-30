import os
import subprocess
import pandas as pd
import sys

# Change working directory to the script's directory
os.chdir(os.path.dirname(os.path.abspath(__file__)))

top_level_input_dir = "/home/bs128/ESKAPE_subsets/Staphylococcus_aureus_200Random_AllGenesExtracted_Primers"
#amplify_script = "/Users/bjornshockey/Desktop/rResolve_p00/scripts/InSilico_amplification.py"
amplify_script = "InSilico_amplification_hammingfixed_Server.py"
analyze_script = "analyze_ampliconV2_Server.py"

for subdir in os.listdir(top_level_input_dir):
    sub_path = os.path.join(top_level_input_dir, subdir)
    if not os.path.isdir(sub_path):
        continue

    print(f"\n[Processing]: {sub_path}")
    all_files = os.listdir(sub_path)
    print(f"  Files in dir: {all_files}")

    candidate_jsons = []
    for f in all_files:
        f = f.strip()
        if f.endswith(".json"):
            stem = f.replace(".json", "")
            csv_file = stem + ".csv"
            csv_path = os.path.join(sub_path, csv_file)
            print(f"  Found JSON: {f} — Checking for CSV: {csv_file} → {'FOUND' if os.path.exists(csv_path) else 'MISSING'}")
            if os.path.exists(csv_path):
                candidate_jsons.append((f, csv_file))

    best_json = None
    max_amplicons = -1

    for json_file, csv_file in candidate_jsons:
        try:
            csv_path = os.path.join(sub_path, csv_file)
            df = pd.read_csv(csv_path)
            df.columns = [col.strip() for col in df.columns]  # normalize column names
            print(f"  Reading: {csv_file} — Columns: {df.columns.tolist()}")

            length_col = [
                col for col in df.columns
                if "length" in col.lower() and ("amplicon" in col.lower() or "effective" in col.lower())
            ]
            print(f"    Detected length columns: {length_col}")
            if not length_col:
                continue
            col = length_col[0]
            count = df[col].apply(pd.to_numeric, errors='coerce').dropna().le(500).sum()
            print(f"    {csv_file} → {count} amplicons ≤ 500 bp")
            if count > max_amplicons:
                max_amplicons = count
                best_json = os.path.join(sub_path, json_file)
        except Exception as e:
            print(f"  Failed to read {csv_file}: {e}")

    fasta_file = next((os.path.join(sub_path, f) for f in all_files if f.endswith(".filt.mc.fasta")), None)
    print(f"  Fasta found: {fasta_file if fasta_file else 'None'}")

    if not best_json or not fasta_file:
        print(f"  Skipping {subdir}: no valid json or fasta found.")
        continue

    print(f"  → Selected: {os.path.basename(best_json)} ({max_amplicons} good amplicons)")

    output_folder = os.path.join(sub_path, "output")
    os.makedirs(output_folder, exist_ok=True)

    subprocess.run(["python", amplify_script, best_json, fasta_file, output_folder], check=True)
    summary_csv = os.path.join(sub_path, "amplicon.summary.csv")
    subprocess.run(["python", analyze_script, output_folder, summary_csv], check=True)

