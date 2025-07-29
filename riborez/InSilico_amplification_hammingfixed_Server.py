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


def find_primer_positions(primer_seq, seq, start_range, end_range, is_reverse=False):
    """
    Find the exact start and end positions of a primer within its designated range.
    Handles gaps in aligned sequences properly.
    
    Args:
        primer_seq: The primer sequence to search for
        seq: The sequence to search in (may contain gaps)
        start_range: Start of the search range
        end_range: End of the search range
        is_reverse: Whether this is a reverse primer (search for reverse complement)
    
    Returns:
        tuple: (primer_start, primer_end, found) where found is a boolean
    """
    # Try multiple approaches for finding the primer
    approaches = []
    
    # Approach 1: Search for reverse complement (original approach)
    if is_reverse:
        approaches.append(reverse_complement(primer_seq))
    
    # Approach 2: Search for the primer as-is (in case PMPrimer already provides correct orientation)
    approaches.append(primer_seq)
    
    # Approach 3: If reverse primer, also try the original sequence
    if is_reverse:
        approaches.append(primer_seq)
    
    # Try each approach
    for search_seq in approaches:
        # Extract the range to search in (keep gaps for accurate positioning)
        range_seq = seq[start_range:end_range + 1]
        
        # Find the primer in the range (gaps are ignored in the search)
        found_pos = -1
        actual_start = None
        actual_end = None
        
        # Search through the range, skipping gaps when matching
        for i in range(len(range_seq) - len(search_seq) + 1):
            match_found = True
            
            # Check if primer matches starting at position i
            for j, primer_char in enumerate(search_seq):
                seq_char = range_seq[i + j]
                
                # Skip gaps and continue matching
                if seq_char == '-':
                    continue
                    
                # Check if characters match
                if primer_char != seq_char:
                    match_found = False
                    break
            
            if match_found:
                found_pos = i
                # Calculate actual positions accounting for gaps
                actual_start = start_range + i
                actual_end = start_range + i + len(search_seq) - 1
                break
        
        if found_pos != -1:
            return actual_start, actual_end, True
    
    # If no match found, try with a more expanded range (±20 positions)
    # This accounts for potential gap-induced shifts in the alignment
    expanded_start = max(0, start_range - 20)
    expanded_end = min(len(seq) - 1, end_range + 20)
    
    for search_seq in approaches:
        range_seq = seq[expanded_start:expanded_end + 1]
        
        for i in range(len(range_seq) - len(search_seq) + 1):
            match_found = True
            
            for j, primer_char in enumerate(search_seq):
                seq_char = range_seq[i + j]
                
                if seq_char == '-':
                    continue
                    
                if primer_char != seq_char:
                    match_found = False
                    break
            
            if match_found:
                actual_start = expanded_start + i
                actual_end = expanded_start + i + len(search_seq) - 1
                return actual_start, actual_end, True
    
    # If still no match found, try searching the entire sequence as a last resort
    # This handles cases where the PMPrimer range is significantly off
    for search_seq in approaches:
        # Search the entire sequence
        for i in range(len(seq) - len(search_seq) + 1):
            match_found = True
            
            for j, primer_char in enumerate(search_seq):
                seq_char = seq[i + j]
                
                if seq_char == '-':
                    continue
                    
                if primer_char != seq_char:
                    match_found = False
                    break
            
            if match_found:
                actual_start = i
                actual_end = i + len(search_seq) - 1
                return actual_start, actual_end, True
    
    return None, None, False


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
                "ForwardVariant", "ReverseVariant", "ErrorStatus"
            ])
            dashless_writer.writerow([
                "Header", "ForwardPrimer", "ReversePrimer",
                "AmpliconStart", "AmpliconEnd", "AmpliconSequence",
                "ForwardVariant", "ReverseVariant", "ErrorStatus"
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

                # Find exact primer positions within their ranges
                f_actual_start, f_actual_end, f_found = find_primer_positions(
                    forward_deg, seq, fstart, fend, is_reverse=False
                )
                r_actual_start, r_actual_end, r_found = find_primer_positions(
                    reverse_deg, seq, rstart, rend, is_reverse=True
                )

                # Debug information
                print(f"  Forward primer: {forward_deg}")
                print(f"  Forward range: {fstart}-{fend}")
                print(f"  Forward found: {f_found} at {f_actual_start}-{f_actual_end}")
                print(f"  Reverse primer: {reverse_deg}")
                print(f"  Reverse range: {rstart}-{rend}")
                print(f"  Reverse found: {r_found} at {r_actual_start}-{r_actual_end}")
                if not r_found:
                    print(f"  WARNING: Reverse primer not found in designated range {rstart}-{rend}")
                    print(f"  Searched expanded range: {max(0, rstart-20)}-{min(len(seq)-1, rend+20)}")
                    print(f"  Also searched entire sequence as fallback")

                # Check for edge cases
                error_status = "OK"
                if not f_found:
                    error_status = "FORWARD_PRIMER_NOT_FOUND"
                elif not r_found:
                    error_status = "REVERSE_PRIMER_NOT_FOUND"
                elif f_actual_end >= r_actual_start:
                    error_status = "PRIMERS_OVERLAP_OR_WRONG_ORDER"

                # Calculate amplicon coordinates using exact primer positions
                if f_found and r_found and f_actual_end < r_actual_start:
                    # Amplicon should be the region between primers (not including primers)
                    # Start after the forward primer ends, end at the start of the reverse primer
                    amp_start = f_actual_end + 1
                    amp_end = r_actual_start
                else:
                    # Fallback to range boundaries if there's an error
                    amp_start = fend + 1
                    amp_end = rstart
                    if error_status == "OK":
                        error_status = "USING_RANGE_BOUNDARIES"

                amplicon_seq_raw = seq[amp_start:amp_end] if amp_start < amp_end else ""
                amplicon_seq = amplicon_seq_raw.replace("-", "")

                preserved_writer.writerow([
                    header, forward_deg, reverse_deg,
                    amp_start, amp_end, amplicon_seq_raw,
                    fmatch, rmatch, error_status
                ])
                dashless_writer.writerow([
                    header, forward_deg, reverse_deg,
                    amp_start, amp_end, amplicon_seq,
                    fmatch, rmatch, error_status
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

