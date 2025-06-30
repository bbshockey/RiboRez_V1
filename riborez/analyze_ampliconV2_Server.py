#!/usr/bin/env python3
import os
import csv
import glob
import statistics
import subprocess
import tempfile
import sys
from Bio import SeqIO

def check_muscle(muscle_path):
    """Check if MUSCLE is installed and executable."""
    if not os.path.exists(muscle_path) or not os.access(muscle_path, os.X_OK):
        raise FileNotFoundError(f"MUSCLE not found or not executable: {muscle_path}")

def align_sequences_muscle(input_fasta, output_fasta, muscle_path):
    """Align sequences using MUSCLE and save as FASTA using proven command-line options."""
    check_muscle(muscle_path)
    command = [muscle_path, "-align", input_fasta, "-output", output_fasta]
    try:
        subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(f"Alignment saved to {output_fasta}")
    except subprocess.CalledProcessError as e:
        print(f"Error running MUSCLE: {e.stderr.decode()}")
        raise

def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Sequences must be equal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def analyze_csv(file_path, muscle_exe=None):
    rows = []
    with open(file_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append((row["Header"], row["AmpliconSequence"]))
    num_input = len(rows)
    if num_input == 0:
        return 0, 0, 0, {}, {}
    if num_input == 1:
        header, seq = rows[0]
        asv_to_headers = {seq: [header]}
        aligned_dict = {"seq0": seq}
        return 1, 1, 0, asv_to_headers, aligned_dict

    asv_to_headers = {}
    aligned_dict = {}
    for header, seq in rows:
        asv_to_headers.setdefault(seq, []).append(header)
        aligned_dict[header] = seq

    num_unique_asvs = len(asv_to_headers)
    distances = []
    aligned_list = list(asv_to_headers.keys())
    for i in range(len(aligned_list)):
        for j in range(i + 1, len(aligned_list)):
            try:
                d = hamming_distance(aligned_list[i], aligned_list[j])
                distances.append(d)
            except ValueError:
                continue
    median_hd = statistics.median(distances) if distances else 0
    return num_input, num_unique_asvs, median_hd, asv_to_headers, aligned_dict


def main(input_folder, output_summary_csv, muscle_exe=None):
    csv_files = glob.glob(os.path.join(input_folder, "*.csv"))
    summary_data = []
    max_asv_count = 0
    # Create folder for aligned FASTA files.
    aligned_folder = os.path.join(os.path.dirname(output_summary_csv), "aligned_amplicons")
    os.makedirs(aligned_folder, exist_ok=True)
    for csv_file in csv_files:
        basename = os.path.basename(csv_file)
        num_input, num_unique_asvs, median_hd, asv_to_headers, aligned_dict = analyze_csv(csv_file, muscle_exe)
        if len(asv_to_headers) > max_asv_count:
            max_asv_count = len(asv_to_headers)
        # Compute amplicon length from the first aligned sequence (assuming all are equal length)
        if aligned_dict:
            amplicon_length = len(next(iter(aligned_dict.values())))
        else:
            amplicon_length = 0

        row_dict = {
            "PrimerPairCSV": basename,
            "NumInputSequences": num_input,
            "NumUniqueASVs": num_unique_asvs,
            "MedianHammingDistance": median_hd,
            "AmpliconLength": amplicon_length
        }
        # Create separate columns for each ASV.
        for idx, (asv_seq, headers) in enumerate(asv_to_headers.items(), start=1):
            row_dict[f"ASV_{idx}"] = f"ASV_{idx}"
            row_dict[f"ASV_{idx}_Sequence"] = asv_seq
            row_dict[f"ASV_{idx}_Headers"] = ", ".join(headers)
        summary_data.append(row_dict)
        # Write the aligned FASTA for this primer group.
        aligned_fasta_path = os.path.join(aligned_folder, f"aligned_{basename}.fasta")
        with open(aligned_fasta_path, "w") as fout:
            for seq_id, seq in aligned_dict.items():
                fout.write(f">{seq_id}\n{seq}\n")
    # Build header columns for the summary CSV.
    fixed_columns = ["PrimerPairCSV", "NumInputSequences", "NumUniqueASVs", "MedianHammingDistance", "AmpliconLength"]
    asv_columns = []
    for i in range(1, max_asv_count + 1):
        asv_columns.extend([f"ASV_{i}", f"ASV_{i}_Sequence", f"ASV_{i}_Headers"])
    all_columns = fixed_columns + asv_columns
    with open(output_summary_csv, "w", newline="") as out_f:
        writer = csv.DictWriter(out_f, fieldnames=all_columns)
        writer.writeheader()
        for row in summary_data:
            writer.writerow(row)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: analyze_amplicons.py <input_folder> <output_summary_csv>")
        sys.exit(1)
    input_folder = sys.argv[1]
    output_summary_csv = sys.argv[2]
    main(input_folder, output_summary_csv)


