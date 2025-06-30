#!/usr/bin/env python3
import sys
import os
import csv
import json


def read_fasta(filename):
    seqs = {}
    header = None
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                header = line[1:].strip()
                seqs[header] = []
            else:
                seqs[header].append(line)
    return {h: "".join(seq) for h, seq in seqs.items()}


def reverse_complement(seq):
    base_map = str.maketrans("ACGT", "TGCA")
    return seq.translate(base_map)[::-1]


def main(json_file, aligned_fasta, out_folder):
    with open(json_file) as jf:
        data = json.load(jf)

    seq_dict = read_fasta(aligned_fasta)
    os.makedirs(out_folder, exist_ok=True)

    summary = []

    for i, (amplicon_key, amplicon_val) in enumerate(data.items(), start=1):
        forward_deg = amplicon_val[0]
        forward_variants = list(amplicon_val[1].keys())
        reverse_deg = amplicon_val[2]
        reverse_variants = list(amplicon_val[3].keys())

        region_str = amplicon_key.strip("()[] ")
        fregion, rregion = region_str.split("], [")
        fstart, fend = map(int, fregion.split(','))
        rstart, rend = map(int, rregion.split(','))

        print(f"\nAmplicon {i}")
        print(f"Forward region: {fstart}-{fend}, Reverse region: {rstart}-{rend}")

        matches = 0

        preserved_filename = f"amplicon_{i}.csv"
        dashless_filename = f"amplicon_{i}_dashless.csv"

        preserved_path = os.path.join(out_folder, preserved_filename)
        dashless_path = os.path.join(os.path.dirname(out_folder), dashless_filename)
  # MAIN directory

        with open(preserved_path, "w", newline="") as preserved_file, open(dashless_path, "w", newline="") as dashless_file:
            preserved_writer = csv.writer(preserved_file)
            dashless_writer = csv.writer(dashless_file)

            preserved_writer.writerow([
                "Header", "ForwardPrimer", "ReversePrimer",
                "AmpliconStart", "AmpliconEnd", "AmpliconSequence",
                "ForwardVariant", "ReverseVariant"
            ])
            dashless_writer.writerow([
                "Header", "ForwardPrimer", "ReversePrimer",
                "AmpliconStart", "AmpliconEnd", "AmpliconSequence",
                "ForwardVariant", "ReverseVariant"
            ])

            for header, seq in seq_dict.items():
                f_slice = seq[max(0, fstart - 1):fend + 2].replace("-", "")
                r_slice = seq[max(0, rstart - 1):rend + 2].replace("-", "")

                print(f"\nSequence: {header}")
                print(f"Forward slice: {f_slice}")
                print(f"Reverse slice: {r_slice}")

                fmatch = next((v for v in forward_variants if v in f_slice), None)
                if not fmatch:
                    print("No forward match")
                    continue

                rmatch = next((v for v in reverse_variants if reverse_complement(v) in r_slice), None)
                if not rmatch:
                    print("No reverse match")
                    continue

                amp_start = fend
                amp_end = rstart
                amplicon_seq_raw = seq[amp_start:amp_end] if amp_start < amp_end else ""
                amplicon_seq = amplicon_seq_raw.replace("-", "")

                preserved_writer.writerow([
                    header, forward_deg, reverse_deg,
                    amp_start, amp_end, amplicon_seq_raw,
                    fmatch, rmatch
                ])
                dashless_writer.writerow([
                    header, forward_deg, reverse_deg,
                    amp_start, amp_end, amplicon_seq,
                    fmatch, rmatch
                ])

                matches += 1

        summary.append((preserved_filename, matches))

    print("\nMatch Summary:")
    for fname, count in summary:
        print(f"{fname}: {count} matching sequences")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: script.py <pmprimer_json> <aligned_fasta> <out_folder>")
        sys.exit(1)
    json_file = sys.argv[1]
    aligned_fasta = sys.argv[2]
    out_folder = sys.argv[3]
    main(json_file, aligned_fasta, out_folder)

