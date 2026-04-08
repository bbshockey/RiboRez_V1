#!/usr/bin/env python3
"""
Ribozyme-Primer Integration for RiboRez
========================================
Finds group-I-intron T-sites near each reverse primer in a RiboRez gene
directory, scores IGS/EGS conservation, and produces a combined ranking of
primer pairs that maximize both taxonomic resolution (amplicon diversity) and
ribozyme cleavage reliability.

Inputs (auto-discovered from a single gene directory):
  - Aligned FASTA  : *.filt.mc.fasta  (from primer-design)
  - PMPrimer JSON  : *_recommand_region_primer.json  (from primer-design)
  - Amplicon CSV   : *_processed_amplicon.csv or amplicon.summary.csv
                     (from amplicon-analysis — required for combined ranking)

Outputs:
  - {gene}_ribozyme_designs.tsv   : one row per (primer_pair, T-site)
  - {gene}_combined_ranked.tsv    : one row per primer pair, rank-1 T-site only
  - ribozyme_design_log.log
"""

import csv
import glob
import json
import os
import re
import statistics
from collections import Counter
from datetime import datetime
from pathlib import Path


# ── Logging ──────────────────────────────────────────────────────────────────

def _log(message, log_file):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(log_file, "a") as fh:
        fh.write(f"[{timestamp}] {message}\n")
    print(f"[{timestamp}] {message}")


# ── FASTA helpers ─────────────────────────────────────────────────────────────

def _read_fasta(path):
    """Return (ordered_names, {name: sequence}) from a FASTA file."""
    seqs, order = {}, []
    cur, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if cur:
                    seqs[cur] = "".join(buf).upper()
                    order.append(cur)
                cur, buf = line[1:].split()[0], []
            else:
                buf.append(line)
    if cur:
        seqs[cur] = "".join(buf).upper()
        order.append(cur)
    return order, seqs


# ── Alignment column helpers ──────────────────────────────────────────────────

def _plurality(chars):
    non_gap = [c for c in chars if c != "-"]
    return Counter(non_gap).most_common(1)[0][0] if non_gap else None


def _consensus_str(seq_list, start, length, aln_len):
    out = []
    for col in range(start, start + length):
        b = _plurality([s[col] for s in seq_list]) if col < aln_len else "?"
        out.append(b or "-")
    return "".join(out)


def _score_p1ext(seq_list, T_pos, aln_len):
    """
    Score the 3bp P1 extension (T_pos+1 to T_pos+3).
    Returns (mm_cols, per_seq, total_mismatches) using the same plurality logic as IGS/EGS.
    Sequences with gaps at any P1 column count as mismatches.
    """
    per_seq = [0] * len(seq_list)
    mm_cols = 0
    for col in range(T_pos + 1, T_pos + 4):
        if col >= aln_len:
            break
        chars = [s[col] for s in seq_list]
        plur = _plurality(chars)
        if plur is None:
            continue
        hit = False
        for i, c in enumerate(chars):
            if c == "-" or c != plur:
                per_seq[i] += 1
                hit = True
        if hit:
            mm_cols += 1
    return mm_cols, per_seq, sum(per_seq)


# ── IGS / EGS scoring ────────────────────────────────────────────────────────

def _score_igs(seq_list, names, T_pos):
    """Score the 5bp IGS region (T_pos-5 to T_pos-1)."""
    per_seq = [0] * len(seq_list)
    mm_cols = 0
    bad = set()
    for col in range(T_pos - 5, T_pos):
        if col < 0:
            continue
        chars = [s[col] for s in seq_list]
        plur = _plurality(chars)
        if plur is None:
            continue
        hit = False
        for i, c in enumerate(chars):
            if c != "-" and c != plur:
                per_seq[i] += 1
                bad.add(i)
                hit = True
        if hit:
            mm_cols += 1
    total = sum(per_seq)
    bad_names = [names[i] for i in sorted(bad)]
    return mm_cols, per_seq, total, bad_names


def _score_egs(seq_list, T_pos, aln_len, egs_start, egs_end):
    """Score the EGS region (T_pos+egs_start to T_pos+egs_end)."""
    per_seq = [0] * len(seq_list)
    mm_cols = 0
    for offset in range(egs_start, egs_end + 1):
        col = T_pos + offset
        if col >= aln_len:
            break
        chars = [s[col] for s in seq_list]
        plur = _plurality(chars)
        if plur is None:
            continue
        hit = False
        for i, c in enumerate(chars):
            if c != "-" and c != plur:
                per_seq[i] += 1
                hit = True
        if hit:
            mm_cols += 1
    total = sum(per_seq)
    mean = statistics.mean(per_seq) if per_seq else 0
    cv = round(statistics.stdev(per_seq) / mean, 4) if mean > 0 else 0.0
    return mm_cols, cv, per_seq, total


# ── T-site search ─────────────────────────────────────────────────────────────

def _find_t_sites(seq_list, names, aln_len, window_start, window_end,
                  egs_start, egs_end, min_t_conservation=1.0):
    """
    Scan alignment columns in [window_start, window_end] for valid T-sites.

    A position is a valid T-site if:
      1. The fraction of non-gap characters equal to 'T' >= min_t_conservation.
         Default 1.0 = 100% (all sequences must have T). Lower values (e.g. 0.9)
         allow a small proportion of non-T sequences, useful for large datasets.
      2. There are at least 5 upstream columns available for the IGS.
      3. The EGS region fits within the alignment.

    P1 extension conservation (T+1 to T+3) is scored rather than hard-filtered,
    so imperfect P1 sites are included but ranked lower.

    Ranking within each primer pair (ascending = better):
      1. igs_per_seq_sum  (IGS conservation — primary)
      2. p1_per_seq_sum   (P1 extension conservation — secondary)
      3. egs_per_seq_sum  (EGS conservation — tertiary)
      4. egs_mm_cols      (tiebreaker)

    Returns a list of candidate dicts, sorted best-first.
    """
    candidates = []
    lo = max(5, window_start)
    hi = min(window_end, aln_len - 1)

    for T_pos in range(lo, hi + 1):
        if T_pos + egs_start >= aln_len:
            continue
        t_chars = [s[T_pos] for s in seq_list if s[T_pos] != "-"]
        if not t_chars:
            continue
        t_fraction = sum(1 for c in t_chars if c == "T") / len(t_chars)
        if t_fraction < min_t_conservation:
            continue

        igs_mm_cols, igs_per_seq, igs_per_seq_sum, mm_seqs = _score_igs(
            seq_list, names, T_pos
        )
        p1_mm_cols, p1_per_seq, p1_per_seq_sum = _score_p1ext(seq_list, T_pos, aln_len)
        egs_mm_cols, _cv, egs_per_seq, egs_per_seq_sum = _score_egs(
            seq_list, T_pos, aln_len, egs_start, egs_end
        )
        candidates.append({
            "aln_pos": T_pos + 1,
            "t_conservation": round(t_fraction, 4),
            "igs": _consensus_str(seq_list, T_pos - 5, 5, aln_len),
            "p1_ext": _consensus_str(seq_list, T_pos + 1, 3, aln_len),
            "igs_mm_cols": igs_mm_cols,
            "igs_per_seq_sum": igs_per_seq_sum,
            "n_igs_mismatch_seqs": len(mm_seqs),
            "per_seq_igs": ",".join(str(x) for x in igs_per_seq),
            "igs_mismatch_seqs": ";".join(mm_seqs) if mm_seqs else "none",
            "p1_mm_cols": p1_mm_cols,
            "p1_per_seq_sum": p1_per_seq_sum,
            "per_seq_p1": ",".join(str(x) for x in p1_per_seq),
            "egs_mm_cols": egs_mm_cols,
            "egs_per_seq_sum": egs_per_seq_sum,
            "per_seq_egs": ",".join(str(x) for x in egs_per_seq),
            "egs_consensus": _consensus_str(
                seq_list, T_pos + egs_start, egs_end - egs_start + 1, aln_len
            ),
        })

    # Ranking: T (hard filter) > IGS > P1 extension > EGS
    candidates.sort(key=lambda r: (
        r["igs_per_seq_sum"],
        r["p1_per_seq_sum"],
        r["egs_per_seq_sum"],
        r["egs_mm_cols"],
    ))
    return candidates


# ── PMPrimer JSON parsing ─────────────────────────────────────────────────────

def _parse_pmprimer_json(json_path):
    with open(json_path) as fh:
        data = json.load(fh)
    pairs = []
    for i, (key, val) in enumerate(data.items(), start=1):
        region_str = key.strip("()[] ")
        fregion, rregion = region_str.split("], [")
        fstart, fend = map(int, fregion.split(","))
        rstart, rend = map(int, rregion.split(","))
        pairs.append({
            "amplicon_id": i,
            "fwd_degenerate": val[0],
            "rev_degenerate": val[2],
            "fstart": fstart, "fend": fend,
            "rstart": rstart, "rend": rend,
        })
    return pairs


# ── Amplicon CSV parsing ──────────────────────────────────────────────────────

def _parse_amplicon_csv(path):
    amplicons = {}
    with open(path) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            primer_csv = row.get("PrimerPairCSV", "")
            m = re.search(r"amplicon_(\d+)", primer_csv)
            if not m:
                continue
            amp_id = int(m.group(1))
            record = {"amplicon_id": amp_id, "PrimerPairCSV": primer_csv}
            for col in ["NumUniqueASVs", "MedianHammingDistance", "AmpliconLength"]:
                try:
                    record[col] = float(row.get(col, "0"))
                except (ValueError, TypeError):
                    record[col] = 0.0
            for col in ["Original#ofSequences", "nonRedundantOriginal#ofSequences",
                        "SequencesSuccessfullyAmplified", "BacteriaAmplified", "UniqueBacteria"]:
                try:
                    record[col] = int(float(row.get(col, "")))
                except (ValueError, TypeError):
                    record[col] = None
            if record.get("SequencesSuccessfullyAmplified") is None:
                try:
                    record["SequencesSuccessfullyAmplified"] = int(float(row.get("NumInputSequences", "0")))
                except (ValueError, TypeError):
                    record["SequencesSuccessfullyAmplified"] = 0
            amplicons[amp_id] = record
    return amplicons


# ── File discovery ────────────────────────────────────────────────────────────

def _find_gene_inputs(gene_dir):
    """
    Auto-discover aligned FASTA, PMPrimer JSON/CSV pairs, and amplicon CSV
    in a gene directory. Searches both the gene dir and pmprimer_outputs/.
    """
    gene_dir = str(gene_dir)
    pmprimer_subdir = os.path.join(gene_dir, "pmprimer_outputs")
    search_dirs = [gene_dir]
    if os.path.isdir(pmprimer_subdir):
        search_dirs.append(pmprimer_subdir)

    # Aligned FASTA — prefer .filt.mc.fasta
    fasta = None
    for sd in search_dirs:
        for f in os.listdir(sd):
            if f.endswith(".filt.mc.fasta"):
                fasta = os.path.join(sd, f)
                break
        if fasta:
            break

    # PMPrimer JSON+CSV pairs
    json_csv_pairs = []
    for sd in search_dirs:
        for f in sorted(os.listdir(sd)):
            if f.endswith(".json"):
                csv_path = os.path.join(sd, f[:-5] + ".csv")
                if os.path.exists(csv_path):
                    json_csv_pairs.append((os.path.join(sd, f), csv_path))

    # Amplicon summary CSV: prefer *_processed_amplicon.csv, fall back to amplicon.summary.csv
    amplicon_csv = None
    processed = glob.glob(os.path.join(gene_dir, "*_processed_amplicon.csv"))
    if processed:
        amplicon_csv = processed[0]
    else:
        raw = os.path.join(gene_dir, "amplicon.summary.csv")
        if os.path.exists(raw):
            amplicon_csv = raw

    return fasta, json_csv_pairs, amplicon_csv


def _match_json_to_amplicon_csv(json_csv_pairs, amplicon_csv_path):
    """
    Return the JSON whose primer pair count matches the number of rows in the
    amplicon CSV (i.e. the JSON that was used for amplicon analysis).
    Falls back to the JSON with the most primer pairs if no exact match found.
    """
    if not json_csv_pairs:
        return None
    if not amplicon_csv_path:
        # No amplicon CSV — return the largest JSON
        best, best_count = None, -1
        for json_path, _ in json_csv_pairs:
            try:
                n = len(_parse_pmprimer_json(json_path))
                if n > best_count:
                    best_count, best = n, json_path
            except Exception:
                continue
        return best

    amplicons = _parse_amplicon_csv(amplicon_csv_path)
    target_count = len(amplicons)

    best, best_count = None, -1
    for json_path, _ in json_csv_pairs:
        try:
            n = len(_parse_pmprimer_json(json_path))
            if n == target_count:
                return json_path   # exact match
            if n > best_count:
                best_count, best = n, json_path
        except Exception:
            continue
    return best   # fallback: largest JSON


# ── Output field definitions ──────────────────────────────────────────────────

_DESIGNS_FIELDS = [
    "amplicon_id", "fwd_primer", "rev_primer",
    "rev_region", "t_site_rank", "aln_pos", "t_conservation",
    "igs", "T", "p1_ext", "p1_loop",
    "igs_mm_cols", "igs_per_seq_sum", "n_igs_mismatch_seqs",
    "per_seq_igs", "igs_mismatch_seqs",
    "p1_mm_cols", "p1_per_seq_sum", "per_seq_p1",
    "egs_mm_cols", "egs_per_seq_sum", "per_seq_egs", "egs_consensus",
]

_AMPLICON_FIELDS = [
    "amplicon_id", "PrimerPairCSV",
    "NumUniqueASVs", "MedianHammingDistance", "AmpliconLength",
]
_EXTENDED_FIELDS = ["SequencesSuccessfullyAmplified", "BacteriaAmplified", "UniqueBacteria"]
_RIBO_OUT_FIELDS = [
    "fwd_primer", "rev_primer", "rev_region", "aln_pos",
    "igs", "T", "p1_ext", "p1_loop",
    "igs_mm_cols", "igs_per_seq_sum", "n_igs_mismatch_seqs",
    "per_seq_igs", "igs_mismatch_seqs",
    "p1_mm_cols", "p1_per_seq_sum", "per_seq_p1",
    "egs_mm_cols", "egs_per_seq_sum", "per_seq_egs", "egs_consensus",
]


# ── Core steps ────────────────────────────────────────────────────────────────

def _run_t_site_search(fasta_path, json_path, output_tsv,
                       window, egs_start, egs_end, p1_loop, log_file,
                       min_t_conservation=1.0):
    """Write ribozyme_designs.tsv — one row per (primer_pair, T-site)."""
    names, seqs = _read_fasta(fasta_path)
    seq_list = [seqs[n] for n in names]
    aln_len = len(seq_list[0])
    primer_pairs = _parse_pmprimer_json(json_path)

    _log(f"  Sequences: {len(names)}, alignment length: {aln_len}", log_file)
    _log(f"  Primer pairs: {len(primer_pairs)}", log_file)

    pairs_with_sites = 0
    total_sites = 0

    with open(output_tsv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=_DESIGNS_FIELDS, delimiter="\t")
        writer.writeheader()
        for pp in primer_pairs:
            amp_id = pp["amplicon_id"]
            rstart, rend = pp["rstart"], pp["rend"]
            win_lo = max(0, rstart - window)
            win_hi = min(aln_len - 1, rend + window)
            candidates = _find_t_sites(
                seq_list, names, aln_len, win_lo, win_hi, egs_start, egs_end,
                min_t_conservation=min_t_conservation,
            )
            if not candidates:
                continue
            pairs_with_sites += 1
            for rank, site in enumerate(candidates, 1):
                total_sites += 1
                writer.writerow({
                    "amplicon_id": amp_id,
                    "fwd_primer": pp["fwd_degenerate"],
                    "rev_primer": pp["rev_degenerate"],
                    "rev_region": f"{rstart}-{rend}",
                    "t_site_rank": rank,
                    "aln_pos": site["aln_pos"],
                    "igs": site["igs"],
                    "T": "T",
                    "p1_ext": site["p1_ext"],
                    "p1_loop": p1_loop,
                    "igs_mm_cols": site["igs_mm_cols"],
                    "igs_per_seq_sum": site["igs_per_seq_sum"],
                    "n_igs_mismatch_seqs": site["n_igs_mismatch_seqs"],
                    "per_seq_igs": site["per_seq_igs"],
                    "igs_mismatch_seqs": site["igs_mismatch_seqs"],
                    "p1_mm_cols": site["p1_mm_cols"],
                    "p1_per_seq_sum": site["p1_per_seq_sum"],
                    "per_seq_p1": site["per_seq_p1"],
                    "egs_mm_cols": site["egs_mm_cols"],
                    "egs_per_seq_sum": site["egs_per_seq_sum"],
                    "per_seq_egs": site["per_seq_egs"],
                    "egs_consensus": site["egs_consensus"],
                })

    _log(
        f"  T-site search: {pairs_with_sites}/{len(primer_pairs)} primer pairs "
        f"had T-sites, {total_sites} total designs",
        log_file,
    )
    return pairs_with_sites, total_sites


def _run_combined_ranking(amplicon_csv, ribozyme_tsv, output_tsv,
                          max_amplicon_length, max_igs_mismatches, log_file):
    """
    Join amplicon resolution data with ribozyme T-site designs, filter, and rank.
    Returns the number of candidates written.
    """
    amplicons = _parse_amplicon_csv(amplicon_csv)

    # Parse ribozyme TSV — keep only rank-1 T-site per primer pair
    sites_by_amp = {}
    with open(ribozyme_tsv) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            amp_id = int(row["amplicon_id"])
            sites_by_amp.setdefault(amp_id, []).append(row)

    ribozymes = {}
    for amp_id, sites in sites_by_amp.items():
        rank1 = [s for s in sites if s["t_site_rank"] == "1"]
        ribozymes[amp_id] = rank1[0] if rank1 else sites[0]

    joined = [
        (amplicons[amp_id], ribozymes[amp_id])
        for amp_id in sorted(set(amplicons) & set(ribozymes))
    ]

    filtered = []
    for amp, ribo in joined:
        if amp["AmpliconLength"] >= max_amplicon_length:
            continue
        if max_igs_mismatches is not None and int(ribo["igs_per_seq_sum"]) > max_igs_mismatches:
            continue
        filtered.append((amp, ribo))

    filtered.sort(key=lambda x: (
        -x[0]["NumUniqueASVs"],
        -x[0]["MedianHammingDistance"],
        int(x[1]["igs_per_seq_sum"]),
        int(x[1]["egs_per_seq_sum"]),
    ))

    _log(
        f"  Combined ranking: {len(amplicons)} amplicons + {len(ribozymes)} designs "
        f"→ {len(joined)} joined → {len(filtered)} after filters",
        log_file,
    )

    if not filtered:
        _log(
            "  No candidates passed filters — try relaxing --max-amplicon-length "
            "or --max-igs-mismatches",
            log_file,
        )
        return 0

    has_extended = any(filtered[0][0].get(f) is not None for f in _EXTENDED_FIELDS)
    out_fields = ["combined_rank"] + _AMPLICON_FIELDS
    if has_extended:
        out_fields += _EXTENDED_FIELDS
    out_fields += _RIBO_OUT_FIELDS

    with open(output_tsv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=out_fields, delimiter="\t")
        writer.writeheader()
        for rank, (amp, ribo) in enumerate(filtered, 1):
            row = {"combined_rank": rank}
            for f in _AMPLICON_FIELDS:
                val = amp.get(f, "")
                if isinstance(val, float) and val == int(val):
                    val = int(val)
                row[f] = val
            if has_extended:
                for f in _EXTENDED_FIELDS:
                    val = amp.get(f)
                    row[f] = val if val is not None else ""
            for f in _RIBO_OUT_FIELDS:
                row[f] = ribo.get(f, "")
            writer.writerow(row)

    # Print top-10 summary to console
    print()
    print("=" * 80)
    print("TOP 10 COMBINED RANKINGS")
    print(f"{'Rank':>4} {'AmpID':>5} {'ASVs':>5} {'HammDist':>8} {'AmpLen':>6} "
          f"| {'Pos':>5} {'IGS':>6} {'IGS_sum':>7} {'EGS_sum':>7}")
    print("-" * 80)
    for rank, (amp, ribo) in enumerate(filtered[:10], 1):
        print(f"{rank:>4} {amp['amplicon_id']:>5} "
              f"{int(amp['NumUniqueASVs']):>5} "
              f"{amp['MedianHammingDistance']:>8.1f} "
              f"{int(amp['AmpliconLength']):>6} "
              f"| {ribo['aln_pos']:>5} "
              f"{ribo['igs']:>6} "
              f"{ribo['igs_per_seq_sum']:>7} "
              f"{ribo['egs_per_seq_sum']:>7}")
    print("=" * 80)

    return len(filtered)


# ── Public entry point ────────────────────────────────────────────────────────

def design_ribozymes(input_folder, output_folder=None, window=10,
                     egs_start=11, egs_end=60, p1_loop="TAACCACA",
                     max_amplicon_length=500, max_igs_mismatches=None,
                     min_t_conservation=1.0):
    """
    Design group-I-intron ribozymes paired with RiboRez primer pairs.

    Args:
        input_folder      : Single gene directory (e.g. Pseudomonas_Primers/16S/)
        output_folder     : Output dir — auto-named {gene}_RibozymeDesign/ if None
        window            : bp to extend search beyond reverse primer (default 10)
        egs_start         : EGS region start offset from T-site (default 11)
        egs_end           : EGS region end offset from T-site (default 60)
        p1_loop           : P1 loop sequence (default TAACCACA)
        max_amplicon_length : Discard candidates with amplicon >= this bp (default 500)
        max_igs_mismatches  : Discard where igs_per_seq_sum > this (default: no filter)
        min_t_conservation  : Minimum fraction of sequences with T at cleavage site
                              (0.0–1.0, default 1.0 = 100% required)

    Returns:
        Path to the output folder.
    """
    input_folder = Path(input_folder)
    if not input_folder.exists():
        raise FileNotFoundError(f"Input folder not found: {input_folder}")

    gene_name = input_folder.name

    if output_folder is None:
        output_folder = input_folder.parent / f"{gene_name}_RibozymeDesign"
    else:
        output_folder = Path(output_folder)

    output_folder.mkdir(parents=True, exist_ok=True)
    log_file = str(output_folder / "ribozyme_design_log.log")

    with open(log_file, "w") as fh:
        fh.write(
            f"=== Ribozyme Design Started at "
            f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} ===\n\n"
        )

    print(f"[INFO] Ribozyme design for gene: {gene_name}")
    print(f"[INFO] Input:  {input_folder}")
    print(f"[INFO] Output: {output_folder}")
    print(f"[INFO] T-site window: reverse primer ±{window} bp")
    print(f"[INFO] EGS window: T+{egs_start} to T+{egs_end}")
    print(f"[INFO] Max amplicon length: {max_amplicon_length} bp")
    print(f"[INFO] Min T-site conservation: {min_t_conservation:.0%}")
    if max_igs_mismatches is not None:
        print(f"[INFO] Max IGS mismatches filter: {max_igs_mismatches}")

    # ── Discover inputs ───────────────────────────────────────────────────────
    fasta, json_csv_pairs, amplicon_csv = _find_gene_inputs(input_folder)

    if not fasta:
        raise FileNotFoundError(
            f"No .filt.mc.fasta found in {input_folder} or {input_folder}/pmprimer_outputs/.\n"
            "Make sure 'riborez primer-design' has been run for this gene."
        )
    if not json_csv_pairs:
        raise FileNotFoundError(
            f"No PMPrimer JSON/CSV pairs found in {input_folder} or "
            f"{input_folder}/pmprimer_outputs/.\n"
            "Make sure 'riborez primer-design' has been run for this gene."
        )

    _log(f"Aligned FASTA: {os.path.basename(fasta)}", log_file)
    _log(f"PMPrimer JSON files found: {len(json_csv_pairs)}", log_file)

    if amplicon_csv:
        _log(f"Amplicon CSV: {os.path.basename(amplicon_csv)}", log_file)
    else:
        _log(
            "WARNING: No amplicon summary CSV found. Run 'riborez amplicon-analysis' "
            "first to enable combined ranking. Ribozyme designs TSV will still be produced.",
            log_file,
        )
        print("[WARN] No amplicon summary CSV found — combined_ranked.tsv will not be produced.")
        print("[WARN] Run 'riborez amplicon-analysis' on this gene folder first.")

    json_path = _match_json_to_amplicon_csv(json_csv_pairs, amplicon_csv)
    if not json_path:
        raise FileNotFoundError("Could not identify a valid PMPrimer JSON to use.")
    _log(f"Using JSON: {os.path.basename(json_path)}", log_file)

    # ── Step 1: T-site search ─────────────────────────────────────────────────
    print(f"\n[INFO] Step 1/2: T-site search...")
    designs_tsv = output_folder / f"{gene_name}_ribozyme_designs.tsv"
    pairs_with_sites, total_sites = _run_t_site_search(
        fasta_path=fasta,
        json_path=json_path,
        output_tsv=str(designs_tsv),
        window=window,
        egs_start=egs_start,
        egs_end=egs_end,
        p1_loop=p1_loop.upper(),
        log_file=log_file,
        min_t_conservation=min_t_conservation,
    )

    if total_sites == 0:
        _log(
            "No T-sites found. Consider increasing --window or checking that the "
            "gene has conserved T positions near reverse primer binding sites.",
            log_file,
        )
        print("[WARN] No T-sites found — ribozyme designs TSV is empty.")
        print(f"[INFO] Output: {output_folder}")
        print(f"[INFO] Log:    {log_file}")
        return output_folder

    print(f"[INFO] {total_sites} ribozyme designs across {pairs_with_sites} primer pairs")
    print(f"[INFO] Ribozyme designs TSV: {designs_tsv}")

    # ── Step 2: Combined ranking ──────────────────────────────────────────────
    if not amplicon_csv:
        print(f"\n[SUCCESS] Ribozyme design complete (no combined ranking — run amplicon-analysis first).")
        print(f"[INFO] Output: {output_folder}")
        print(f"[INFO] Log:    {log_file}")
        return output_folder

    print(f"\n[INFO] Step 2/2: Combined ranking...")
    combined_tsv = output_folder / f"{gene_name}_combined_ranked.tsv"
    n_candidates = _run_combined_ranking(
        amplicon_csv=amplicon_csv,
        ribozyme_tsv=str(designs_tsv),
        output_tsv=str(combined_tsv),
        max_amplicon_length=max_amplicon_length,
        max_igs_mismatches=max_igs_mismatches,
        log_file=log_file,
    )

    print(f"\n[SUCCESS] Ribozyme design complete!")
    print(f"[INFO] {n_candidates} candidates in combined ranking")
    print(f"[INFO] Ribozyme designs TSV:  {designs_tsv}")
    if n_candidates > 0:
        print(f"[INFO] Combined ranked TSV:   {combined_tsv}")
    print(f"[INFO] Log:                   {log_file}")

    return output_folder
