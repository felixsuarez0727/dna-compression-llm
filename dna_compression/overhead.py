"""Determine a compression dictionary overhead from sequence statistics."""

import argparse
import json
import sys
from collections import Counter


DEFAULT_MAX_OVERHEAD = 200
DEFAULT_TOP_KMERS = 20
TOKEN_BYTES = 3

SAMPLE_JSON_TEMPLATE = {
    "sequence": "",
    "priority": 1,
    "count": 9999,
    "potential_savings": 999999,
}


def load_sequences(path: str) -> list[str]:
    sequences = []
    try:
        with open(path, "r", encoding="utf-8") as input_handle:
            for line in input_handle:
                sequence = line.strip()
                if sequence:
                    sequences.append(sequence)
    except FileNotFoundError:
        print(f"ERROR: file not found: {path}", file=sys.stderr)
        sys.exit(1)
    return sequences


def compute_top_kmers(
    sequences: list[str], k_sizes=(4, 6, 8, 10, 12), top_n: int = DEFAULT_TOP_KMERS
) -> dict:
    kmer_counts = {kmer_size: Counter() for kmer_size in k_sizes}
    for sequence in sequences:
        for kmer_size in k_sizes:
            for index in range(len(sequence) - kmer_size + 1):
                kmer = sequence[index:index + kmer_size]
                if kmer.count("N") <= kmer_size // 4:
                    kmer_counts[kmer_size][kmer] += 1
    return {
        kmer_size: kmer_counts[kmer_size].most_common(top_n)
        for kmer_size in kmer_sizes(kmer_counts)
    }


def kmer_sizes(kmer_counts):
    return kmer_counts.keys()


def validate_and_count(pattern: str, sequences: list[str]) -> int:
    total = 0
    for sequence in sequences:
        start = 0
        while True:
            position = sequence.find(pattern, start)
            if position == -1:
                break
            total += 1
            start = position + 1
    return total


def find_tandem_unit(sequence: str, max_unit: int = 12):
    length = len(sequence)
    for unit_length in range(1, min(max_unit + 1, length // 2 + 1)):
        unit = sequence[:unit_length]
        repeats = length // unit_length
        if unit * repeats == sequence[:unit_length * repeats] and length % unit_length == 0 and repeats >= 2:
            return unit, repeats
    return None, 0


def canonical_unit(unit: str) -> str:
    rotations = [unit[index:] + unit[:index] for index in range(len(unit))]
    return min(rotations)


def build_candidate_patterns(
    sequences: list[str], top_kmers: dict, min_repeats: int = 2, max_repeats: int = 12
) -> dict[str, int]:
    candidates: dict[str, int] = {}
    checked_units: set[str] = set()

    for _, items in top_kmers.items():
        for kmer, _ in items:
            unit, _ = find_tandem_unit(kmer)
            base = unit if unit else kmer

            if base in checked_units:
                continue
            checked_units.add(base)

            for repeats in range(max_repeats, min_repeats - 1, -1):
                candidate = base * repeats
                count = validate_and_count(candidate, sequences)
                if count > 0:
                    if candidate not in candidates or candidates[candidate] < count:
                        candidates[candidate] = count

    return candidates


def deduplicate_phase_variants(patterns: dict[str, int]) -> dict[str, int]:
    groups: dict[str, list] = {}
    for sequence, count in patterns.items():
        unit, _ = find_tandem_unit(sequence)
        canonical = canonical_unit(unit) if unit else sequence
        groups.setdefault(canonical, []).append((sequence, count))

    result: dict[str, int] = {}
    for members in groups.values():
        by_length: dict[int, tuple] = {}
        for sequence, count in members:
            length = len(sequence)
            if length not in by_length or count > by_length[length][1]:
                by_length[length] = (sequence, count)
        for sequence, count in by_length.values():
            result[sequence] = count
    return result


def measure_json_entry_bytes(seq_len: int) -> int:
    entry = dict(SAMPLE_JSON_TEMPLATE)
    entry["sequence"] = "A" * seq_len
    raw = json.dumps({"!": entry}, indent=4)
    inner = raw[2:-2]
    return len(inner.encode("utf-8"))


def compute_empirical_fixed_overhead() -> float:
    lengths = [8, 12, 16, 20, 24, 32, 40, 48, 64, 80, 96, 128]
    differences = []
    for length in lengths:
        total = measure_json_entry_bytes(length)
        differences.append(total - length)
    return sum(differences) / len(differences)


def simulate_overhead(patterns: dict[str, int], overhead: int) -> dict:
    gross_savings = 0
    dictionary_cost = 0
    retained = 0

    for sequence, count in patterns.items():
        sequence_length = len(sequence)
        savings_per_hit = sequence_length - TOKEN_BYTES
        potential = count * savings_per_hit
        threshold = sequence_length + overhead
        if savings_per_hit <= 0:
            continue
        if potential > threshold:
            retained += 1
            gross_savings += potential
            dictionary_cost += measure_json_entry_bytes(sequence_length)

    net_savings = gross_savings - dictionary_cost
    return {
        "overhead": overhead,
        "retained": retained,
        "gross_savings": gross_savings,
        "dict_cost": dictionary_cost,
        "net_savings": net_savings,
    }


def bar(value: int, max_value: int, width: int = 40) -> str:
    if max_value == 0:
        return ""
    filled = round(width * value / max_value)
    return "█" * filled + "░" * (width - filled)


def report(
    results: list[dict], empirical_overhead: float, optimal_overhead: int, total_raw_bytes: int
) -> None:
    separator = "=" * 72
    print(f"\n{separator}")
    print("  OVERHEAD ANALYSIS REPORT")
    print(separator)
    print(f"  Total raw bytes in sequence file : {total_raw_bytes:,}")
    print(f"  Token format                     : <X>  ({TOKEN_BYTES} bytes per substitution)")
    print(f"  Savings formula used             : count × (seq_len - {TOKEN_BYTES})  [corrected]")
    print("  NOTE: dna-compress detect uses count × (seq_len - 1)  ← underestimates by 2×count bytes")
    print(f"  Empirical fixed JSON overhead    : {empirical_overhead:.1f} bytes/entry")
    print(f"  ► Recommended --overhead value   : {optimal_overhead}")
    print("    (default in pattern_detector   : 5)")
    print(separator)

    print(
        f"\n{'OVH':>5}  {'Kept':>5}  {'Gross Savings':>14}  {'Dict Cost':>10}  {'Net Savings':>12}  Chart"
    )
    print("-" * 72)

    max_net = max(result["net_savings"] for result in results) if results else 1
    for result in results:
        marker = " ◄ OPTIMAL" if result["overhead"] == optimal_overhead else ""
        chart = bar(max(result["net_savings"], 0), max_net, width=20)
        print(
            f"{result['overhead']:>5}  "
            f"{result['retained']:>5}  "
            f"{result['gross_savings']:>14,}  "
            f"{result['dict_cost']:>10,}  "
            f"{result['net_savings']:>12,}  "
            f"{chart}{marker}"
        )

    print(separator)
    print("\nINTERPRETATION")
    print("-" * 72)
    print(
        f"  Each pattern substitution writes <TOKEN> ({TOKEN_BYTES} bytes) in the compressed\n"
        f"  file, so the real saving per occurrence is (seq_len - {TOKEN_BYTES}), not (seq_len - 1).\n"
        "\n"
        f"  overhead={optimal_overhead} keeps the patterns whose corrected savings fully\n"
        f"  cover their real JSON dictionary cost (~{empirical_overhead:.0f} bytes fixed metadata\n"
        "  + sequence length), maximising net compression gain.\n"
        "\n"
        "  Use it like this:\n"
        f"    dna-compress detect ... --overhead {optimal_overhead}\n"
        "\n"
        "  ALSO RECOMMENDED: update optimize_patterns() in the detector:\n"
        "    Change:  p['potential_savings'] = p['count'] * (seq_len - 1)\n"
        "    To:      p['potential_savings'] = p['count'] * (seq_len - 3)\n"
    )
    print(separator)


def configure_utf8_stdout():
    reconfigure = getattr(sys.stdout, "reconfigure", None)
    if reconfigure:
        reconfigure(encoding="utf-8")


def run_analysis(file_path, max_overhead=DEFAULT_MAX_OVERHEAD, step=5, top=DEFAULT_TOP_KMERS):
    configure_utf8_stdout()
    print(f"\nLoading sequences from: {file_path}")
    sequences = load_sequences(file_path)
    total_raw_bytes = sum(len(sequence) for sequence in sequences)
    print(f"  {len(sequences):,} sequences loaded  ({total_raw_bytes:,} raw bytes)")

    print("\nComputing k-mer frequencies (k = 4, 6, 8, 10, 12)...")
    top_kmers = compute_top_kmers(sequences, top_n=top)
    for kmer_size in (4, 8, 12):
        top3 = ", ".join(f"{kmer}({count})" for kmer, count in top_kmers[kmer_size][:3])
        print(f"  k={kmer_size}: {top3}")

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

    print(f"\nSweeping overhead 1 → {max_overhead} (step {step})...")
    sweep_values = list(range(1, max_overhead + 1, step))
    for extra in [5, round(empirical_overhead)]:
        if extra not in sweep_values and 1 <= extra <= max_overhead:
            sweep_values.append(extra)
    sweep_values.sort()

    results = [simulate_overhead(candidates, overhead) for overhead in sweep_values]
    best = max(results, key=lambda result: result["net_savings"])
    report(results, empirical_overhead, best["overhead"], total_raw_bytes)


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Find the optimal --overhead value for dna-compress detect",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--file", "-f", required=True, help="Input .seq.txt file (one sequence per line)")
    parser.add_argument(
        "--max-overhead",
        type=int,
        default=DEFAULT_MAX_OVERHEAD,
        help="Upper bound of the overhead sweep range",
    )
    parser.add_argument(
        "--step",
        type=int,
        default=5,
        help="Step size for the overhead sweep (smaller = finer resolution)",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=DEFAULT_TOP_KMERS,
        help="Top-N k-mers per k-size to use as seed patterns",
    )
    return parser.parse_args(argv)


def main(argv=None) -> None:
    args = parse_args(argv)
    run_analysis(args.file, args.max_overhead, args.step, args.top)