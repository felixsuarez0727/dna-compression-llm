"""Pattern detection pipeline preserving the existing provider behavior."""

import argparse
import json
import re
import sys
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from ..config import MICRO_TOKEN_POOL, PROVIDER_DEFAULTS
from ..logging_utils import save_log as write_log
from .local_models import (
    analyze_batch_dnabert2 as run_dnabert2_batch,
    analyze_batch_hyenadna as run_hyenadna_batch,
    analyze_synthesis_dnabert2 as run_dnabert2_synthesis,
    analyze_synthesis_hyenadna as run_hyenadna_synthesis,
)
from .providers import build_client, call_llm


_LOG_FILE = None


def save_log(line):
    write_log(line, filename=_LOG_FILE)


def analyze_batch_dnabert2(client, sequences, batch_num):
    return run_dnabert2_batch(client, sequences, batch_num, save_log)


def analyze_batch_hyenadna(client, sequences, batch_num):
    return run_hyenadna_batch(client, sequences, batch_num, save_log)


def analyze_synthesis_dnabert2(client, all_candidates, top_kmers):
    return run_dnabert2_synthesis(client, all_candidates, top_kmers)


def analyze_synthesis_hyenadna(client, all_candidates, top_kmers):
    return run_hyenadna_synthesis(
        client,
        all_candidates,
        top_kmers,
        find_tandem_unit,
        save_log,
    )


def parse_seq_sequences(file_path):
    try:
        with open(file_path, "r") as input_handle:
            for line in input_handle:
                sequence = line.strip()
                if sequence:
                    yield sequence
    except FileNotFoundError:
        save_log(f"Error: File '{file_path}' not found.")
        sys.exit(1)


def compute_top_kmers(sequences, k_sizes=(4, 6, 8, 10, 12), top_n=20):
    kmer_counts = {kmer_size: Counter() for kmer_size in k_sizes}
    for sequence in sequences:
        for kmer_size in k_sizes:
            for index in range(len(sequence) - kmer_size + 1):
                kmer = sequence[index:index + kmer_size]
                if kmer.count("N") <= kmer_size // 4:
                    kmer_counts[kmer_size][kmer] += 1
    return {
        kmer_size: kmer_counts[kmer_size].most_common(top_n)
        for kmer_size in k_sizes
    }


def format_kmers_for_prompt(top_kmers, max_per_k=10):
    lines = []
    for kmer_size in sorted(top_kmers.keys()):
        items = top_kmers[kmer_size][:max_per_k]
        if items:
            lines.append(f"\n  Length {kmer_size}:")
            for kmer, count in items:
                lines.append(f"    {kmer}  ->  {count} occurrences")
    return "\n".join(lines)


def find_tandem_unit(sequence, max_unit=12):
    length = len(sequence)
    for unit_length in range(1, min(max_unit + 1, length // 2 + 1)):
        unit = sequence[:unit_length]
        repeats = length // unit_length
        remainder = length % unit_length
        if unit * repeats == sequence[:unit_length * repeats] and remainder == 0 and repeats >= 2:
            return unit, repeats
    return None, 0


def canonical_unit(unit):
    rotations = [unit[index:] + unit[:index] for index in range(len(unit))]
    return min(rotations)


def canonicalize_pattern(sequence):
    unit, repeats = find_tandem_unit(sequence)
    if unit:
        canonical = canonical_unit(unit)
        return canonical * repeats
    return sequence


def validate_and_count_pattern(pattern, sequences):
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


def expand_to_tandem_variants(base_unit, sequences, min_repeats=2, max_repeats=12):
    found = {}
    for repeats in range(max_repeats, min_repeats - 1, -1):
        candidate = base_unit * repeats
        count = validate_and_count_pattern(candidate, sequences)
        if count > 0:
            found[candidate] = count
    return found


def parse_model_response(text):
    patterns = {}
    if not text:
        return patterns
    lines = text.strip().split("\n")
    pattern_re = re.compile(r"\b([ACGTN]{4,})\b\D{0,30}?(\d+)", re.IGNORECASE)
    for line in lines:
        if len(line) > 200:
            continue
        match = pattern_re.search(line)
        if match:
            sequence = match.group(1).upper()
            count = int(match.group(2))
            if sequence not in patterns or count > patterns[sequence]:
                patterns[sequence] = count
    return patterns


def analyze_batch_with_context(client, provider, model_id, sequences, batch_num, top_kmers):
    if provider == "dnabert2":
        return analyze_batch_dnabert2(client, sequences, batch_num)
    if provider == "hyenadna":
        return analyze_batch_hyenadna(client, sequences, batch_num)

    sequences_block = "\n".join(sequences)
    kmer_context = format_kmers_for_prompt(top_kmers)

    user_prompt = (
        "You are an expert in analyzing synthetic DNA sequences built from TANDEM REPEATS.\n\n"
        "== CONTEXT ==\n"
        "A tandem repeat is a short unit (motif) repeated consecutively.\n"
        "Example: motif ATCG repeated 4 times = ATCGATCGATCGATCG\n\n"
        "== REAL DATASET STATISTICS (most frequent k-mers) ==\n"
        f"{kmer_context}\n\n"
        "== YOUR TASK ==\n"
        "Analyze the sequences below and identify:\n"
        "1. BASE UNITS (the shortest repeating motif)\n"
        "2. FULL SEQUENCES formed by those repetitions\n"
        "3. Variants with an NNN prefix (e.g. NNNATCGATCG)\n"
        "4. Estimate how many times each pattern appears in the shown sequences\n\n"
        "== RULES ==\n"
        "- Only use characters A, C, G, T, N\n"
        "- Minimum pattern length: 4 characters\n"
        "- Do NOT invent patterns that are not present\n"
        "- Report BOTH the base unit and its extended repeated form\n"
        "  Example: if you see ATCGATCGATCGATCG, report both ATCG and ATCGATCGATCGATCG\n\n"
        "== RESPONSE FORMAT (mandatory, nothing else) ==\n"
        "PATTERN:COUNT\n\n"
        f"Sequences:\n{sequences_block}\n\n"
        "Patterns:"
    )

    try:
        return call_llm(client, provider, model_id, "You are an expert bioinformatics AI.", user_prompt)
    except Exception as error:
        save_log(f"\n  Warning: Error in batch {batch_num}: {error}")
        return None


def analyze_batch_synthesis(client, provider, model_id, all_candidates, top_kmers):
    if provider == "dnabert2":
        return analyze_synthesis_dnabert2(client, all_candidates, top_kmers)
    if provider == "hyenadna":
        return analyze_synthesis_hyenadna(client, all_candidates, top_kmers)

    candidates_block = "\n".join(
        f"{sequence}: {count} times"
        for sequence, count in sorted(all_candidates.items(), key=lambda item: -item[1])[:80]
    )
    kmer_context = format_kmers_for_prompt(top_kmers, max_per_k=5)

    user_prompt = (
        "You are an expert in DNA compression via repetitive pattern substitution.\n\n"
        "== CANDIDATE PATTERNS FOUND ==\n"
        f"{candidates_block}\n\n"
        "== STATISTICAL CONTEXT ==\n"
        f"{kmer_context}\n\n"
        "== SYNTHESIS TASK ==\n"
        "Review the candidates and produce the OPTIMAL list for compression:\n"
        "1. Include the LONGEST patterns\n"
        "2. Include the BASE UNIT of each tandem repeat\n"
        "3. Include NNN-prefixed variants if detected\n"
        "4. Remove redundant duplicates\n"
        "5. Remove any pattern that does not make sense as a tandem repeat\n\n"
        "Priority: length > frequency\n\n"
        "== FORMAT (nothing else) ==\n"
        "PATTERN:COUNT\n\n"
        "Optimized pattern list:"
    )

    try:
        return call_llm(client, provider, model_id, "You are an expert DNA compression AI.", user_prompt)
    except Exception as error:
        save_log(f"\n  Warning: Error in synthesis pass: {error}")
        return None


def token_length_pnn(index):
    return len(f"P{index}")


def best_token(sequence, count, pnn_index, micro_available):
    sequence_length = len(sequence)
    pnn_token = f"P{pnn_index}"
    gain_pnn = (sequence_length - len(pnn_token)) * count
    gain_micro = (sequence_length - 1) * count

    if micro_available and gain_micro > gain_pnn:
        return micro_available[0], gain_micro, "micro"
    if gain_pnn > 0:
        return pnn_token, gain_pnn, "pnn"
    if micro_available and gain_micro > 0:
        return micro_available[0], gain_micro, "micro"
    return None, 0, "discard"


def deduplicate_phase_variants(aggregated_patterns):
    groups = {}

    for sequence, count in aggregated_patterns.items():
        unit, repeats = find_tandem_unit(sequence)
        if unit:
            canonical = canonical_unit(unit)
            groups.setdefault(canonical, []).append((sequence, count))
        else:
            groups.setdefault(sequence, []).append((sequence, count))

    result = {}
    removed = []

    for canonical, members in groups.items():
        if len(members) == 1:
            sequence, count = members[0]
            result[sequence] = count
        else:
            by_length = {}
            for sequence, count in members:
                length = len(sequence)
                if length not in by_length or count > by_length[length][1]:
                    by_length[length] = (sequence, count)

            kept = set()
            for length, (sequence, count) in by_length.items():
                result[sequence] = count
                kept.add(sequence)

            for sequence, count in members:
                if sequence not in kept:
                    removed.append((sequence, count, canonical))

    if removed:
        save_log(f"  Removed {len(removed)} phase-variant duplicates:")
        for sequence, count, canonical in removed[:8]:
            save_log(f"    '{sequence[:40]}' (unit canon='{canonical}', count={count})")
        if len(removed) > 8:
            save_log(f"    ... and {len(removed) - 8} more")

    return result


def build_compressor_patterns(aggregated_patterns):
    pre_sorted = sorted(
        aggregated_patterns.items(),
        key=lambda item: len(item[0]) * item[1],
        reverse=True,
    )

    result = {}
    pnn_index = 1
    micro_pool = list(MICRO_TOKEN_POOL)
    discarded = []
    micro_used = []

    for sequence, count in pre_sorted:
        token, gain, token_type = best_token(sequence, count, pnn_index, micro_pool)

        if token_type == "discard":
            discarded.append((sequence, count))
            continue

        if token_type == "micro":
            micro_char = micro_pool.pop(0)
            micro_used.append((micro_char, sequence))
            result[micro_char] = {
                "sequence": sequence,
                "count": count,
                "token": micro_char,
                "token_type": "micro",
                "bytes_saved": gain,
                "priority": None,
            }
        else:
            result[f"P{pnn_index}"] = {
                "sequence": sequence,
                "count": count,
                "token": f"P{pnn_index}",
                "token_type": "pnn",
                "bytes_saved": gain,
                "priority": None,
            }
            pnn_index += 1

    sorted_by_length = sorted(
        result.items(), key=lambda item: len(item[1]["sequence"]), reverse=True
    )
    for rank, (key, _) in enumerate(sorted_by_length, start=1):
        result[key]["priority"] = rank

    total_saved = sum(value["bytes_saved"] for value in result.values())
    pnn_count = sum(1 for value in result.values() if value["token_type"] == "pnn")
    micro_count = sum(1 for value in result.values() if value["token_type"] == "micro")

    save_log(f"\n  Patterns retained  : {len(result)} (PNN: {pnn_count}, micro: {micro_count})")
    save_log(f"  Micro-token slots  : {micro_count}/{len(MICRO_TOKEN_POOL)} used")
    if micro_used:
        save_log(
            "  Micro assignments  : "
            + ", ".join(f"'{character}'='{sequence}'" for character, sequence in micro_used)
        )
    if discarded:
        save_log("  Discarded (no gain): " + ", ".join(sequence for sequence, _ in discarded))
    save_log(f"  Estimated bytes saved: {total_saved:,}")

    return result


def optimize_patterns(compressor_patterns, dictionary_overhead=5):
    patterns = []
    for key, data in compressor_patterns.items():
        if "token" not in data:
            data["token"] = key
        patterns.append(data)

    for pattern in patterns:
        sequence_length = len(pattern["sequence"])
        pattern["potential_savings"] = pattern["count"] * (sequence_length - 1)

    optimized = [
        pattern
        for pattern in patterns
        if pattern["potential_savings"] > (len(pattern["sequence"]) + dictionary_overhead)
    ]

    optimized.sort(
        key=lambda pattern: (len(pattern["sequence"]), pattern["potential_savings"]),
        reverse=True,
    )

    used_chars = set("ACGT")
    safe_tokens = [
        chr(index)
        for index in range(33, 127)
        if chr(index) not in used_chars and chr(index) not in ("<", ">", '"', "\\")
    ]

    if len(optimized) > len(safe_tokens):
        save_log(f"  Limiting to the best {len(safe_tokens)} patterns due to token capacity.")
        optimized = optimized[:len(safe_tokens)]

    result = {}
    for index, pattern in enumerate(optimized):
        token = safe_tokens[index]
        result[token] = {
            "sequence": pattern["sequence"],
            "priority": index + 1,
            "count": pattern["count"],
            "potential_savings": pattern["potential_savings"],
        }

    save_log(f"  Final patterns after optimization: {len(result)}")
    save_log("  Priority 1 is the longest/most efficient pattern.")

    return result


def save_to_json(data, filename):
    Path(filename).parent.mkdir(parents=True, exist_ok=True)
    with open(filename, "w", encoding="utf-8") as output_handle:
        json.dump(data, output_handle, indent=4, ensure_ascii=False)


def run_detection(args):
    model_id = args.model or PROVIDER_DEFAULTS[args.provider]

    if args.model and args.provider in ("dnabert2", "hyenadna"):
        PROVIDER_DEFAULTS[args.provider] = args.model

    client = build_client(args.provider, args.key, save_log)

    timestamp = time.strftime("%Y%m%d_%H%M%S")
    global _LOG_FILE
    _LOG_FILE = f"{Path(args.output).stem}_log_{timestamp}.txt"

    separator = "=" * 56
    save_log(f"\n{separator}")
    save_log(f"  FASTQ Pattern Detector  [{args.provider}]")
    save_log(f"{separator}")
    save_log(f"  File       : {args.file}")
    save_log(f"  Model      : {model_id}")
    save_log(f"  Batch size : {args.batch_size}")
    save_log(f"  Overhead   : {args.overhead}")
    if args.provider in ("chatgpt", "deepseek"):
        save_log(f"  Threads    : {args.threads}")
    save_log(f"  Log file   : data/logs/{_LOG_FILE}")
    save_log(f"  Output     : {args.output}")
    save_log(f"{separator}\n")

    save_log("Loading sequences...")
    all_sequences = list(parse_seq_sequences(args.file))
    save_log(f"{len(all_sequences)} sequences loaded")

    save_log("\nK-mer pre-analysis...")
    top_kmers = compute_top_kmers(all_sequences, k_sizes=(4, 6, 8, 10, 12))
    for kmer_size in (4, 8, 12):
        if top_kmers.get(kmer_size):
            top3 = ", ".join(
                f"{kmer}({count})" for kmer, count in top_kmers[kmer_size][:3]
            )
            save_log(f"k={kmer_size}: {top3}")

    save_log("\nBatch analysis...")
    aggregated_patterns = {}

    all_batches = []
    for index in range(0, len(all_sequences), args.batch_size):
        if args.max_batches > 0 and len(all_batches) >= args.max_batches:
            save_log(f"Batch limit of {args.max_batches} reached.")
            break
        all_batches.append(
            (len(all_batches) + 1, all_sequences[index:index + args.batch_size])
        )

    if args.provider in ("dnabert2", "hyenadna"):
        for batch_num, batch in all_batches:
            raw = analyze_batch_with_context(
                client, args.provider, model_id, batch, batch_num, top_kmers
            )
            if raw:
                parsed = parse_model_response(raw)
                save_log(f"  Batch #{batch_num} -> {len(parsed)} patterns detected")
                for pattern, count in parsed.items():
                    aggregated_patterns[pattern] = aggregated_patterns.get(pattern, 0) + count
            else:
                save_log(f"  Batch #{batch_num} -> no result")

    elif args.provider in ("chatgpt", "deepseek") and args.threads > 1:
        save_log(f"  Using {args.threads} parallel threads for {len(all_batches)} batches...")

        def run_batch(batch_num, batch):
            thread_client = build_client(args.provider, args.key, save_log)
            raw = analyze_batch_with_context(
                thread_client, args.provider, model_id, batch, batch_num, top_kmers
            )
            return batch_num, raw

        with ThreadPoolExecutor(max_workers=args.threads) as executor:
            futures_map = {
                executor.submit(run_batch, batch_num, batch): batch_num
                for batch_num, batch in all_batches
            }

            for future in as_completed(futures_map):
                batch_num, raw = future.result()
                if raw:
                    parsed = parse_model_response(raw)
                    save_log(f"  Batch #{batch_num} -> {len(parsed)} patterns detected")
                    for pattern, count in parsed.items():
                        aggregated_patterns[pattern] = (
                            aggregated_patterns.get(pattern, 0) + count
                        )
                else:
                    save_log(f"  Batch #{batch_num} -> no result")
    else:
        for batch_num, batch in all_batches:
            save_log(f"Batch #{batch_num} ({len(batch)} sequences)... ")
            raw = analyze_batch_with_context(
                client, args.provider, model_id, batch, batch_num, top_kmers
            )
            if raw:
                parsed = parse_model_response(raw)
                save_log(f"  -> {len(parsed)} patterns detected")
                for pattern, count in parsed.items():
                    aggregated_patterns[pattern] = aggregated_patterns.get(pattern, 0) + count
            else:
                save_log("  -> no result")
            if batch_num < len(all_batches):
                time.sleep(2)

    if aggregated_patterns:
        save_log(f"\nGlobal synthesis pass ({len(aggregated_patterns)} candidates)...")
        raw_syn = analyze_batch_synthesis(
            client, args.provider, model_id, aggregated_patterns, top_kmers
        )
        if raw_syn:
            parsed_syn = parse_model_response(raw_syn)
            save_log(f"Synthesis contributed {len(parsed_syn)} patterns")
            for pattern, count in parsed_syn.items():
                aggregated_patterns[pattern] = aggregated_patterns.get(pattern, 0) + count

    if not args.no_expand:
        save_log("\nExpanding tandem repeat variants...")
        checked_units = set()
        new_from_expansion = {}

        for pattern in list(aggregated_patterns.keys()):
            unit, repeats = find_tandem_unit(pattern)
            if unit and unit not in checked_units:
                checked_units.add(unit)
                variants = expand_to_tandem_variants(unit, all_sequences)
                for variant, count in variants.items():
                    if variant not in aggregated_patterns or aggregated_patterns[variant] < count:
                        new_from_expansion[variant] = count

        aggregated_patterns.update(new_from_expansion)
        save_log(f"Total after expansion: {len(aggregated_patterns)} patterns")

    if not args.no_validate:
        save_log("\nValidating patterns against real sequences...")
        validated = {}
        fake_count = 0

        for pattern, llm_count in aggregated_patterns.items():
            real_count = validate_and_count_pattern(pattern, all_sequences)
            if real_count > 0:
                validated[pattern] = real_count
            else:
                fake_count += 1

        save_log(f"Valid patterns: {len(validated)}")
        save_log(f"Discarded: {fake_count}")
        aggregated_patterns = validated

    save_log("\nDeduplicating phase-variant patterns...")
    before = len(aggregated_patterns)
    aggregated_patterns = deduplicate_phase_variants(aggregated_patterns)
    save_log(f"{before} -> {len(aggregated_patterns)} patterns after deduplication")

    compressor_ready = build_compressor_patterns(aggregated_patterns)

    save_log("\nOptimizing and assigning final tokens...")
    final_patterns = optimize_patterns(compressor_ready, args.overhead)
    save_to_json(final_patterns, args.output)

    save_log(f"\n{separator}")
    save_log("ANALYSIS COMPLETE")
    save_log(f"Final patterns: {len(final_patterns)}")
    save_log(f"Saved to: {args.output}")
    save_log(f"Log saved to: data/logs/{_LOG_FILE}")
    save_log(f"{separator}")


def add_arguments(parser):
    parser.add_argument("--file", "-f", required=True, help="Path to the sequence text file")
    parser.add_argument(
        "--key",
        "-k",
        default="local",
        help="API key for the selected provider (not required for dnabert2/hyenadna)",
    )
    parser.add_argument(
        "--provider",
        "-p",
        required=True,
        choices=["deepseek", "chatgpt", "gemini", "dnabert2", "hyenadna"],
        help="LLM/DNA model provider to use",
    )
    parser.add_argument("--output", "-o", default="patterns.json", help="Output JSON file")
    parser.add_argument("--batch_size", "-b", type=int, default=30, help="Sequences per API call")
    parser.add_argument("--model", "-m", default=None, help="Model ID (defaults to provider default)")
    parser.add_argument("--max_batches", type=int, default=0, help="Max number of LLM batches")
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Parallel threads for chatgpt and deepseek only (default: 4)",
    )
    parser.add_argument(
        "--overhead",
        type=int,
        default=5,
        help="Dictionary entry cost in bytes for final optimization (default: 5)",
    )
    parser.add_argument("--no_validate", action="store_true", help="Disable validation against real sequences")
    parser.add_argument("--no_expand", action="store_true", help="Disable tandem repeat expansion")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="FASTQ Pattern Detector using LLMs")
    add_arguments(parser)
    return parser.parse_args(argv)


def main(argv=None):
    run_detection(parse_args(argv))