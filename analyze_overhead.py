#!/usr/bin/env python3
"""
analyze_overhead.py
===================
Determines the optimal --overhead value to pass to pattern_detector.py.

The overhead parameter controls which patterns survive the final filter:
    A pattern is KEPT only if:
        count × (seq_len - TOKEN_BYTES)  >  seq_len + overhead

    TOKEN_BYTES = 3  because every substitution uses the format <X>
    (1 token char + 2 delimiter chars '<' and '>').

    NOTE: pattern_detector.py currently uses (seq_len - 1) in its formula,
    which underestimates the token cost by 2 bytes per occurrence.
    This script uses the corrected (seq_len - 3) formula.

    where `overhead` = the fixed JSON metadata cost per entry
    (everything that is NOT the sequence string itself).

This script:
  1. Loads the .seq.txt file and computes k-mer frequencies.
  2. Expands tandem repeat variants deterministically (no LLM).
  3. Measures the REAL JSON entry size for entries at varying sequence lengths.
  4. Sweeps overhead 1 → MAX_OVERHEAD and computes, for every value:
       · patterns retained
       · gross potential savings (bytes, using corrected seq_len - 3)
       · real dictionary cost (bytes, measured from actual JSON structure)
       · NET savings = gross - dict_cost
  5. Reports the overhead that maximises net savings.
  6. Also reports the EMPIRICAL overhead = measured avg fixed JSON metadata cost.

Usage:
    python analyze_overhead.py --file ./data/ERR15993673_5000.seq.txt
    python analyze_overhead.py --file ./data/my.seq.txt --max-overhead 300 --top 30
"""

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

# ──────────────────────────────────────────────────────────────────────────────
# Config
# ──────────────────────────────────────────────────────────────────────────────

DEFAULT_MAX_OVERHEAD = 200
DEFAULT_TOP_KMERS    = 20

# Each pattern substitution produces <TOKEN> in the compressed output.
# '<' + single_char_token + '>' = 3 bytes.
# pattern_detector.py currently uses (seq_len - 1), ignoring the 2 delimiter bytes.
TOKEN_BYTES = 3   # cost of <X> in the compressed file

# These mirror the JSON fields written by optimize_patterns() in pattern_detector.py
# We measure empirically how many bytes the non-sequence fields cost.
SAMPLE_JSON_TEMPLATE = {
    "sequence":          "",        # will be filled
    "priority":          1,
    "count":             9999,
    "potential_savings": 999999,
}


# ──────────────────────────────────────────────────────────────────────────────
# Helpers – copied / adapted from pattern_detector.py (no LLM dependency)
# ──────────────────────────────────────────────────────────────────────────────

def load_sequences(path: str) -> list[str]:
    seqs = []
    try:
        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                s = line.strip()
                if s:
                    seqs.append(s)
    except FileNotFoundError:
        print(f"ERROR: file not found: {path}", file=sys.stderr)
        sys.exit(1)
    return seqs


def compute_top_kmers(sequences: list[str],
                      k_sizes=(4, 6, 8, 10, 12),
                      top_n: int = DEFAULT_TOP_KMERS) -> dict:
    kmer_counts = {k: Counter() for k in k_sizes}
    for seq in sequences:
        for k in k_sizes:
            for i in range(len(seq) - k + 1):
                kmer = seq[i:i + k]
                if kmer.count("N") <= k // 4:
                    kmer_counts[k][kmer] += 1
    return {k: kmer_counts[k].most_common(top_n) for k in k_sizes}


def validate_and_count(pattern: str, sequences: list[str]) -> int:
    total = 0
    for seq in sequences:
        start = 0
        while True:
            pos = seq.find(pattern, start)
            if pos == -1:
                break
            total += 1
            start = pos + 1
    return total


def find_tandem_unit(sequence: str, max_unit: int = 12):
    n = len(sequence)
    for unit_len in range(1, min(max_unit + 1, n // 2 + 1)):
        unit = sequence[:unit_len]
        repeats = n // unit_len
        if unit * repeats == sequence[:unit_len * repeats] and n % unit_len == 0 and repeats >= 2:
            return unit, repeats
    return None, 0


def canonical_unit(unit: str) -> str:
    rotations = [unit[i:] + unit[:i] for i in range(len(unit))]
    return min(rotations)


def build_candidate_patterns(sequences: list[str],
                              top_kmers: dict,
                              min_repeats: int = 2,
                              max_repeats: int = 12) -> dict[str, int]:
    """
    Build candidate patterns by expanding tandem repeats from the top k-mers.
    Returns {sequence: count}.
    """
    candidates: dict[str, int] = {}
    checked_units: set[str] = set()

    for k_size, items in top_kmers.items():
        for kmer, _ in items:
            unit, _ = find_tandem_unit(kmer)
            base = unit if unit else kmer

            if base in checked_units:
                continue
            checked_units.add(base)

            # Try all repeat counts
            for r in range(max_repeats, min_repeats - 1, -1):
                candidate = base * r
                count = validate_and_count(candidate, sequences)
                if count > 0:
                    if candidate not in candidates or candidates[candidate] < count:
                        candidates[candidate] = count

    return candidates


def deduplicate_phase_variants(patterns: dict[str, int]) -> dict[str, int]:
    groups: dict[str, list] = {}
    for seq, count in patterns.items():
        unit, _ = find_tandem_unit(seq)
        canon = canonical_unit(unit) if unit else seq
        groups.setdefault(canon, []).append((seq, count))

    result: dict[str, int] = {}
    for canon, members in groups.items():
        by_length: dict[int, tuple] = {}
        for seq, count in members:
            ln = len(seq)
            if ln not in by_length or count > by_length[ln][1]:
                by_length[ln] = (seq, count)
        for seq, count in by_length.values():
            result[seq] = count
    return result


# ──────────────────────────────────────────────────────────────────────────────
# JSON entry cost measurement
# ──────────────────────────────────────────────────────────────────────────────

def measure_json_entry_bytes(seq_len: int) -> int:
    """
    Measure the actual number of bytes a single JSON entry occupies when written
    with json.dumps(..., indent=4), using the same fields as optimize_patterns().
    """
    entry = dict(SAMPLE_JSON_TEMPLATE)
    entry["sequence"] = "A" * seq_len
    # Simulate a top-level dict with one entry (key = "!" = 1 char token)
    wrapper = {"!": entry}
    raw = json.dumps(wrapper, indent=4)
    # Subtract the outer braces `{\n` and `\n}` (2 + 1 + 1 = 4 chars, plus newlines)
    # We want only the cost of the inner entry including key
    inner = raw[2:-2]   # strips leading "{\n" and "\n}"
    return len(inner.encode("utf-8"))


def compute_empirical_fixed_overhead() -> float:
    """
    Measure JSON entry bytes at multiple sequence lengths, then derive the
    fixed overhead = total_bytes - seq_len (i.e. the non-sequence portion).
    Returns the average fixed overhead across lengths.
    """
    lengths = [8, 12, 16, 20, 24, 32, 40, 48, 64, 80, 96, 128]
    diffs = []
    for ln in lengths:
        total = measure_json_entry_bytes(ln)
        fixed = total - ln          # subtract the sequence characters
        diffs.append(fixed)
    return sum(diffs) / len(diffs)


# ──────────────────────────────────────────────────────────────────────────────
# Simulation
# ──────────────────────────────────────────────────────────────────────────────

def simulate_overhead(patterns: dict[str, int],
                      overhead: int) -> dict:
    """
    Mirror of optimize_patterns() from pattern_detector.py.
    Returns stats for a given overhead value.
    """
    gross_savings  = 0
    dict_cost      = 0
    retained       = 0

    for seq, count in patterns.items():
        seq_len          = len(seq)
        # Corrected: each substitution costs TOKEN_BYTES ("<X>"), not 1 byte
        savings_per_hit  = seq_len - TOKEN_BYTES          # bytes saved per occurrence
        potential        = count * savings_per_hit
        threshold        = seq_len + overhead
        # Patterns of length <= TOKEN_BYTES can never save bytes — skip immediately
        if savings_per_hit <= 0:
            continue
        if potential > threshold:
            retained      += 1
            gross_savings += potential
            dict_cost     += measure_json_entry_bytes(seq_len)

    net_savings = gross_savings - dict_cost
    return {
        "overhead":      overhead,
        "retained":      retained,
        "gross_savings": gross_savings,
        "dict_cost":     dict_cost,
        "net_savings":   net_savings,
    }


# ──────────────────────────────────────────────────────────────────────────────
# Reporting
# ──────────────────────────────────────────────────────────────────────────────

def bar(value: int, max_value: int, width: int = 40) -> str:
    if max_value == 0:
        return ""
    filled = round(width * value / max_value)
    return "█" * filled + "░" * (width - filled)


def report(results: list[dict],
           empirical_overhead: float,
           optimal_overhead: int,
           total_raw_bytes: int) -> None:

    sep = "=" * 72
    print(f"\n{sep}")
    print("  OVERHEAD ANALYSIS REPORT")
    print(sep)
    print(f"  Total raw bytes in sequence file : {total_raw_bytes:,}")
    print(f"  Token format                     : <X>  ({TOKEN_BYTES} bytes per substitution)")
    print(f"  Savings formula used             : count × (seq_len - {TOKEN_BYTES})  [corrected]")
    print(f"  NOTE: pattern_detector.py uses count × (seq_len - 1)  ← underestimates by 2×count bytes")
    print(f"  Empirical fixed JSON overhead    : {empirical_overhead:.1f} bytes/entry")
    print(f"  ► Recommended --overhead value   : {optimal_overhead}")
    print(f"    (default in pattern_detector   : 5)")
    print(sep)

    # Table header
    print(f"\n{'OVH':>5}  {'Kept':>5}  {'Gross Savings':>14}  {'Dict Cost':>10}  {'Net Savings':>12}  Chart")
    print("-" * 72)

    max_net = max(r["net_savings"] for r in results) if results else 1

    for r in results:
        marker = " ◄ OPTIMAL" if r["overhead"] == optimal_overhead else ""
        b = bar(max(r["net_savings"], 0), max_net, width=20)
        print(
            f"{r['overhead']:>5}  "
            f"{r['retained']:>5}  "
            f"{r['gross_savings']:>14,}  "
            f"{r['dict_cost']:>10,}  "
            f"{r['net_savings']:>12,}  "
            f"{b}{marker}"
        )

    print(sep)
    print("\nINTERPRETATION")
    print("-" * 72)
    print(
        f"  Each pattern substitution writes <TOKEN> ({TOKEN_BYTES} bytes) in the compressed\n"
        f"  file, so the real saving per occurrence is (seq_len - {TOKEN_BYTES}), not (seq_len - 1).\n"
        f"\n"
        f"  overhead={optimal_overhead} keeps the patterns whose corrected savings fully\n"
        f"  cover their real JSON dictionary cost (~{empirical_overhead:.0f} bytes fixed metadata\n"
        f"  + sequence length), maximising net compression gain.\n"
        f"\n"
        f"  Use it like this:\n"
        f"    python pattern_detector.py ... --overhead {optimal_overhead}\n"
        f"\n"
        f"  ALSO RECOMMENDED: update optimize_patterns() in pattern_detector.py:\n"
        f"    Change:  p['potential_savings'] = p['count'] * (seq_len - 1)\n"
        f"    To:      p['potential_savings'] = p['count'] * (seq_len - 3)\n"
    )
    print(sep)


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Find the optimal --overhead value for pattern_detector.py",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--file", "-f", required=True,
        help="Input .seq.txt file (one sequence per line)",
    )
    p.add_argument(
        "--max-overhead", type=int, default=DEFAULT_MAX_OVERHEAD,
        help="Upper bound of the overhead sweep range",
    )
    p.add_argument(
        "--step", type=int, default=5,
        help="Step size for the overhead sweep (smaller = finer resolution)",
    )
    p.add_argument(
        "--top", type=int, default=DEFAULT_TOP_KMERS,
        help="Top-N k-mers per k-size to use as seed patterns",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    print(f"\nLoading sequences from: {args.file}")
    sequences = load_sequences(args.file)
    total_raw_bytes = sum(len(s) for s in sequences)
    print(f"  {len(sequences):,} sequences loaded  ({total_raw_bytes:,} raw bytes)")

    print("\nComputing k-mer frequencies (k = 4, 6, 8, 10, 12)...")
    top_kmers = compute_top_kmers(sequences, top_n=args.top)
    for k in (4, 8, 12):
        top3 = ", ".join(f"{km}({c})" for km, c in top_kmers[k][:3])
        print(f"  k={k}: {top3}")

    print("\nExpanding tandem repeat candidates (no LLM)...")
    candidates = build_candidate_patterns(sequences, top_kmers)
    print(f"  {len(candidates):,} raw candidates found")

    print("Deduplicating phase variants...")
    candidates = deduplicate_phase_variants(candidates)
    print(f"  {len(candidates):,} unique patterns after deduplication")

    if not candidates:
        print("  WARNING: no patterns found — try a different file or lower --top threshold.")
        sys.exit(0)

    print("\nMeasuring empirical JSON entry overhead...")
    empirical_overhead = compute_empirical_fixed_overhead()
    print(f"  Fixed metadata cost per entry: ~{empirical_overhead:.1f} bytes")

    print(f"\nSweeping overhead 1 → {args.max_overhead} (step {args.step})...")
    sweep_values = list(range(1, args.max_overhead + 1, args.step))
    # Always include 5 (the default) and the empirical value for comparison
    for extra in [5, round(empirical_overhead)]:
        if extra not in sweep_values and 1 <= extra <= args.max_overhead:
            sweep_values.append(extra)
    sweep_values.sort()

    results = []
    for ov in sweep_values:
        results.append(simulate_overhead(candidates, ov))

    # Find the overhead that maximises net savings
    best = max(results, key=lambda r: r["net_savings"])
    optimal_overhead = best["overhead"]

    report(results, empirical_overhead, optimal_overhead, total_raw_bytes)


if __name__ == "__main__":
    main()
