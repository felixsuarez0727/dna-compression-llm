import argparse
import time
import sys
import json
import re
import os
from pathlib import Path
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed

MICRO_TOKEN_POOL = list("bdefhijklmoprsvwxyz")

_LOG_FILE = None
_LOG_DIR = "data/logs"


def save_log(line):
    safe_line = line.replace("→", "->")
    print(safe_line)
    if _LOG_FILE:
        os.makedirs(_LOG_DIR, exist_ok=True)
        with open(f"{_LOG_DIR}/{_LOG_FILE}", "a", encoding="utf-8") as f:
            f.write(safe_line + "\n")

PROVIDER_DEFAULTS = {
    "deepseek": "deepseek-chat",
    "chatgpt": "gpt-4o-mini",
    "gemini": "gemini-2.5-flash-lite",
}


def build_client(provider, api_key):
    if provider in ("deepseek", "chatgpt"):
        from openai import OpenAI
        if provider == "deepseek":
            return OpenAI(api_key=api_key, base_url="https://api.deepseek.com")
        return OpenAI(api_key=api_key)
    if provider == "gemini":
        from google import genai
        return genai.Client(api_key=api_key)
    raise ValueError(f"Unknown provider: {provider}")


def call_llm(client, provider, model_id, system_prompt, user_prompt):
    if provider in ("deepseek", "chatgpt"):
        
        response = client.chat.completions.create(
            model=model_id,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_prompt},
            ],
            temperature=0,
        )
        return response.choices[0].message.content
    if provider == "gemini":
        response = client.models.generate_content(
            model=model_id,
            contents=user_prompt,
        )
        return response.text
    raise ValueError(f"Unknown provider: {provider}")


def parse_seq_sequences(file_path):
    try:
        with open(file_path, "r") as f:
            for line in f:
                sequence = line.strip()
                if sequence:
                    yield sequence
    except FileNotFoundError:
        save_log(f"Error: File '{file_path}' not found.")
        sys.exit(1)


def compute_top_kmers(sequences, k_sizes=(4, 6, 8, 10, 12), top_n=20):
    kmer_counts = {k: Counter() for k in k_sizes}
    for seq in sequences:
        for k in k_sizes:
            for i in range(len(seq) - k + 1):
                kmer = seq[i:i + k]
                if kmer.count("N") <= k // 4:
                    kmer_counts[k][kmer] += 1
    return {k: kmer_counts[k].most_common(top_n) for k in k_sizes}


def format_kmers_for_prompt(top_kmers, max_per_k=10):
    lines = []
    for k in sorted(top_kmers.keys()):
        items = top_kmers[k][:max_per_k]
        if items:
            lines.append(f"\n  Length {k}:")
            for kmer, count in items:
                lines.append(f"    {kmer}  ->  {count} occurrences")
    return "\n".join(lines)


def find_tandem_unit(sequence, max_unit=12):
    n = len(sequence)
    for unit_len in range(1, min(max_unit + 1, n // 2 + 1)):
        unit = sequence[:unit_len]
        repeats = n // unit_len
        remainder = n % unit_len
        if unit * repeats == sequence[:unit_len * repeats] and remainder == 0 and repeats >= 2:
            return unit, repeats
    return None, 0


def canonical_unit(unit):
    rotations = [unit[i:] + unit[:i] for i in range(len(unit))]
    return min(rotations)


def canonicalize_pattern(sequence):
    unit, reps = find_tandem_unit(sequence)
    if unit:
        canon = canonical_unit(unit)
        return canon * reps
    return sequence


def validate_and_count_pattern(pattern, sequences):
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


def expand_to_tandem_variants(base_unit, sequences, min_repeats=2, max_repeats=12):
    found = {}
    for r in range(max_repeats, min_repeats - 1, -1):
        candidate = base_unit * r
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
            seq = match.group(1).upper()
            count = int(match.group(2))
            if seq not in patterns or count > patterns[seq]:
                patterns[seq] = count
    return patterns


def is_retriable_error(error):
    msg = str(error).lower()
    retriable_markers = (
        "503",
        "unavailable",
        "rate limit",
        "429",
        "timeout",
        "timed out",
        "connection reset",
        "temporarily unavailable",
    )
    return any(marker in msg for marker in retriable_markers)


def analyze_batch_with_context(
    client,
    provider,
    model_id,
    sequences,
    batch_num,
    top_kmers,
    max_retries=3,
    retry_base_delay=2.0,
):
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

    for attempt in range(max_retries + 1):
        try:
            return call_llm(client, provider, model_id, "You are an expert bioinformatics AI.", user_prompt)
        except Exception as e:
            if attempt < max_retries and is_retriable_error(e):
                delay = retry_base_delay * (2 ** attempt)
                save_log(
                    f"\n  Warning: Retriable error in batch {batch_num}: {e}. "
                    f"Retrying in {delay:.1f}s ({attempt + 1}/{max_retries})"
                )
                time.sleep(delay)
                continue
            save_log(f"\n  Warning: Error in batch {batch_num}: {e}")
            return None


def analyze_batch_synthesis(
    client,
    provider,
    model_id,
    all_candidates,
    top_kmers,
    max_retries=3,
    retry_base_delay=2.0,
):
    candidates_block = "\n".join(
        f"{seq}: {count} times"
        for seq, count in sorted(all_candidates.items(), key=lambda x: -x[1])[:80]
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

    for attempt in range(max_retries + 1):
        try:
            return call_llm(client, provider, model_id, "You are an expert DNA compression AI.", user_prompt)
        except Exception as e:
            if attempt < max_retries and is_retriable_error(e):
                delay = retry_base_delay * (2 ** attempt)
                save_log(
                    f"\n  Warning: Retriable error in synthesis pass: {e}. "
                    f"Retrying in {delay:.1f}s ({attempt + 1}/{max_retries})"
                )
                time.sleep(delay)
                continue
            save_log(f"\n  Warning: Error in synthesis pass: {e}")
            return None


def token_length_pnn(idx):
    return len(f"P{idx}")


def best_token(seq, count, pnn_idx, micro_available):
    slen = len(seq)
    pnn_tok = f"P{pnn_idx}"
    gain_pnn = (slen - len(pnn_tok)) * count
    gain_micro = (slen - 1) * count

    if micro_available and gain_micro > gain_pnn:
        return micro_available[0], gain_micro, "micro"
    elif gain_pnn > 0:
        return pnn_tok, gain_pnn, "pnn"
    elif micro_available and gain_micro > 0:
        return micro_available[0], gain_micro, "micro"
    else:
        return None, 0, "discard"


def deduplicate_phase_variants(aggregated_patterns):
    groups = {}

    for seq, count in aggregated_patterns.items():
        unit, reps = find_tandem_unit(seq)
        if unit:
            canon = canonical_unit(unit)
            groups.setdefault(canon, []).append((seq, count))
        else:
            groups.setdefault(seq, []).append((seq, count))

    result = {}
    removed = []

    for canon, members in groups.items():
        if len(members) == 1:
            seq, count = members[0]
            result[seq] = count
        else:
            by_length = {}
            for seq, count in members:
                length = len(seq)
                if length not in by_length or count > by_length[length][1]:
                    by_length[length] = (seq, count)

            kept = set()
            for length, (seq, count) in by_length.items():
                result[seq] = count
                kept.add(seq)

            for seq, count in members:
                if seq not in kept:
                    removed.append((seq, count, canon))

    if removed:
        save_log(f"  Removed {len(removed)} phase-variant duplicates:")
        for seq, count, canon in removed[:8]:
            save_log(f"    '{seq[:40]}' (unit canon='{canon}', count={count})")
        if len(removed) > 8:
            save_log(f"    ... and {len(removed) - 8} more")

    return result


def prune_shifted_redundant_patterns(
    aggregated_patterns,
    min_length=80,
    max_length_gap=24,
    count_tolerance=0.20,
):
    sorted_items = sorted(
        aggregated_patterns.items(),
        key=lambda x: (len(x[0]), x[1]),
        reverse=True,
    )

    kept = []
    pruned = []
    result = {}

    for seq, count in sorted_items:
        should_prune = False

        if len(seq) >= min_length:
            for super_seq, super_count in kept:
                if len(super_seq) <= len(seq):
                    continue
                if len(super_seq) - len(seq) > max_length_gap:
                    continue
                if seq not in super_seq:
                    continue

                max_count = max(super_count, count, 1)
                rel_diff = abs(super_count - count) / max_count
                if rel_diff <= count_tolerance:
                    should_prune = True
                    pruned.append((seq, count, super_seq[:40]))
                    break

        if should_prune:
            continue

        result[seq] = count
        kept.append((seq, count))

    if pruned:
        save_log(f"  Pruned {len(pruned)} shifted redundant patterns:")
        for seq, count, super_preview in pruned[:8]:
            save_log(f"    '{seq[:40]}' (count={count}) shadowed by '{super_preview}...'")
        if len(pruned) > 8:
            save_log(f"    ... and {len(pruned) - 8} more")

    return result


def build_compressor_patterns(aggregated_patterns):
    pre_sorted = sorted(
        aggregated_patterns.items(),
        key=lambda x: len(x[0]) * x[1],
        reverse=True,
    )

    result = {}
    pnn_idx = 1
    micro_pool = list(MICRO_TOKEN_POOL)
    discarded = []
    micro_used = []

    for seq, count in pre_sorted:
        token, gain, ttype = best_token(seq, count, pnn_idx, micro_pool)

        if ttype == "discard":
            discarded.append((seq, count))
            continue

        if ttype == "micro":
            micro_char = micro_pool.pop(0)
            micro_used.append((micro_char, seq))
            result[micro_char] = {
                "sequence": seq,
                "count": count,
                "token": micro_char,
                "token_type": "micro",
                "bytes_saved": gain,
                "priority": None,
            }
        else:
            result[f"P{pnn_idx}"] = {
                "sequence": seq,
                "count": count,
                "token": f"P{pnn_idx}",
                "token_type": "pnn",
                "bytes_saved": gain,
                "priority": None,
            }
            pnn_idx += 1

    sorted_by_len = sorted(result.items(), key=lambda x: len(x[1]["sequence"]), reverse=True)
    for rank, (key, _) in enumerate(sorted_by_len, start=1):
        result[key]["priority"] = rank

    total_saved = sum(v["bytes_saved"] for v in result.values())
    pnn_count = sum(1 for v in result.values() if v["token_type"] == "pnn")
    micro_count = sum(1 for v in result.values() if v["token_type"] == "micro")

    save_log(f"\n  Patterns retained  : {len(result)} (PNN: {pnn_count}, micro: {micro_count})")
    save_log(f"  Micro-token slots  : {micro_count}/{len(MICRO_TOKEN_POOL)} used")
    if micro_used:
        save_log("  Micro assignments  : " + ", ".join(f"'{c}'='{s}'" for c, s in micro_used))
    if discarded:
        save_log("  Discarded (no gain): " + ", ".join(f"'{s}'" for s, _ in discarded))
    save_log(f"  Estimated bytes saved: {total_saved:,}")

    return result


def optimize_patterns(compressor_patterns, dictionary_overhead=5):
    patterns = []
    for key, data in compressor_patterns.items():
        if "token" not in data:
            data["token"] = key
        patterns.append(data)

    for p in patterns:
        seq_len = len(p["sequence"])
        p["potential_savings"] = p["count"] * (seq_len - 1)

    optimized = [
        p for p in patterns
        if p["potential_savings"] > (len(p["sequence"]) + dictionary_overhead)
    ]

    optimized.sort(
        key=lambda x: (len(x["sequence"]), x["potential_savings"]),
        reverse=True,
    )

    used_chars = set("ACGT")
    safe_tokens = [
        chr(i) for i in range(33, 127)
        if chr(i) not in used_chars and chr(i) not in ("<", ">", '"', "\\")
    ]

    if len(optimized) > len(safe_tokens):
        save_log(f"  Limiting to the best {len(safe_tokens)} patterns due to token capacity.")
        optimized = optimized[:len(safe_tokens)]

    result = {}
    for i, p in enumerate(optimized):
        token = safe_tokens[i]
        result[token] = {
            "sequence": p["sequence"],
            "priority": i + 1,
            "count": p["count"],
            "potential_savings": p["potential_savings"],
        }

    save_log(f"  Final patterns after optimization: {len(result)}")
    save_log("  Priority 1 is the longest/most efficient pattern.")

    return result


def save_to_json(data, filename):
    Path(filename).parent.mkdir(parents=True, exist_ok=True)
    with open(filename, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=4, ensure_ascii=False)


def main():
    parser = argparse.ArgumentParser(description="FASTQ Pattern Detector using LLMs")
    parser.add_argument("--file", "-f", required=True, help="Path to the sequence text file")
    parser.add_argument("--key", "-k", required=True, help="API key for the selected provider")
    parser.add_argument(
        "--provider", "-p",
        required=True,
        choices=["deepseek", "chatgpt", "gemini"],
        help="LLM provider to use",
    )
    parser.add_argument("--output", "-o", default="patterns.json", help="Output JSON file")
    parser.add_argument("--batch_size", "-b", type=int, default=30, help="Sequences per API call")
    parser.add_argument("--model", "-m", default=None, help="Model ID (defaults to provider default)")
    parser.add_argument("--max_batches", type=int, default=0, help="Max number of LLM batches")
    parser.add_argument("--threads", type=int, default=4, help="Parallel threads for batch analysis (default: 4)")
    parser.add_argument("--overhead", type=int, default=5, help="Dictionary entry cost in bytes for final optimization (default: 5)")
    parser.add_argument("--max_retries", type=int, default=3, help="Retry attempts for transient API errors (default: 3)")
    parser.add_argument("--retry_base_delay", type=float, default=2.0, help="Base delay in seconds for exponential backoff (default: 2.0)")
    parser.add_argument("--no_validate", action="store_true", help="Disable validation against real sequences")
    parser.add_argument("--no_expand", action="store_true", help="Disable tandem repeat expansion")
    parser.add_argument("--no_shift_prune", action="store_true", help="Disable shifted redundancy pruning")
    parser.add_argument("--shift_prune_min_len", type=int, default=80, help="Minimum pattern length to apply shifted pruning (default: 80)")
    parser.add_argument("--shift_prune_max_gap", type=int, default=24, help="Max length gap between overlapping shifted patterns (default: 24)")
    parser.add_argument("--shift_prune_count_tol", type=float, default=0.20, help="Relative count tolerance for shifted pruning (default: 0.20)")
    parser.add_argument("--log-dir", default=None, help="Directory for detector logs (default: data/logs)")

    args = parser.parse_args()

    model_id = args.model or PROVIDER_DEFAULTS[args.provider]
    client = build_client(args.provider, args.key)

    ts = time.strftime("%Y%m%d_%H%M%S")
    global _LOG_FILE, _LOG_DIR
    if args.log_dir:
        _LOG_DIR = args.log_dir
    _LOG_FILE = f"{Path(args.output).stem}_log_{ts}.txt"

    sep = "=" * 56
    save_log(f"\n{sep}")
    save_log(f"  FASTQ Pattern Detector  [{args.provider}]")
    save_log(f"{sep}")
    save_log(f"  File       : {args.file}")
    save_log(f"  Model      : {model_id}")
    save_log(f"  Batch size : {args.batch_size}")
    save_log(f"  Overhead   : {args.overhead}")
    if args.threads > 1:
        save_log(f"  Threads    : {args.threads}")
    save_log(f"  Log file   : {_LOG_DIR}/{_LOG_FILE}")
    save_log(f"  Output     : {args.output}")
    save_log(f"{sep}\n")

    save_log("Loading sequences...")
    all_sequences = list(parse_seq_sequences(args.file))
    save_log(f"{len(all_sequences)} sequences loaded")

    save_log("\nK-mer pre-analysis...")
    top_kmers = compute_top_kmers(all_sequences, k_sizes=(4, 6, 8, 10, 12))
    for k in (4, 8, 12):
        if top_kmers.get(k):
            top3 = ", ".join(f"{km}({c})" for km, c in top_kmers[k][:3])
            save_log(f"k={k}: {top3}")

    save_log("\nBatch analysis...")
    aggregated_patterns = {}

    all_batches = []
    for i in range(0, len(all_sequences), args.batch_size):
        if args.max_batches > 0 and len(all_batches) >= args.max_batches:
            save_log(f"Batch limit of {args.max_batches} reached.")
            break
        all_batches.append((len(all_batches) + 1, all_sequences[i:i + args.batch_size]))

    if args.threads > 1:
        save_log(f"  Using {args.threads} parallel threads for {len(all_batches)} batches...")

        def run_batch(batch_num, batch):
            thread_client = build_client(args.provider, args.key)
            raw = analyze_batch_with_context(
                thread_client,
                args.provider,
                model_id,
                batch,
                batch_num,
                top_kmers,
                args.max_retries,
                args.retry_base_delay,
            )
            return batch_num, raw

        with ThreadPoolExecutor(max_workers=args.threads) as executor:
            futures_map = {executor.submit(run_batch, batch_num, batch): batch_num
                          for batch_num, batch in all_batches}

            for future in as_completed(futures_map):
                batch_num, raw = future.result()
                if raw:
                    parsed = parse_model_response(raw)
                    save_log(f"  Batch #{batch_num} -> {len(parsed)} patterns detected")
                    for pat, count in parsed.items():
                        aggregated_patterns[pat] = aggregated_patterns.get(pat, 0) + count
                else:
                    save_log(f"  Batch #{batch_num} -> no result")
    else:
        for batch_num, batch in all_batches:
            save_log(f"Batch #{batch_num} ({len(batch)} sequences)... ")
            raw = analyze_batch_with_context(
                client,
                args.provider,
                model_id,
                batch,
                batch_num,
                top_kmers,
                args.max_retries,
                args.retry_base_delay,
            )
            if raw:
                parsed = parse_model_response(raw)
                save_log(f"  -> {len(parsed)} patterns detected")
                for pat, count in parsed.items():
                    aggregated_patterns[pat] = aggregated_patterns.get(pat, 0) + count
            else:
                save_log("  -> no result")
            if batch_num < len(all_batches):
                time.sleep(2)

    if aggregated_patterns:
        save_log(f"\nGlobal synthesis pass ({len(aggregated_patterns)} candidates)...")
        raw_syn = analyze_batch_synthesis(
            client,
            args.provider,
            model_id,
            aggregated_patterns,
            top_kmers,
            args.max_retries,
            args.retry_base_delay,
        )
        if raw_syn:
            parsed_syn = parse_model_response(raw_syn)
            save_log(f"Synthesis contributed {len(parsed_syn)} patterns")
            for pat, count in parsed_syn.items():
                aggregated_patterns[pat] = aggregated_patterns.get(pat, 0) + count

    if not args.no_expand:
        save_log("\nExpanding tandem repeat variants...")
        checked_units = set()
        new_from_expansion = {}

        for pattern in list(aggregated_patterns.keys()):
            unit, reps = find_tandem_unit(pattern)
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

    if not args.no_shift_prune:
        save_log("\nPruning shifted redundant patterns...")
        before_shift = len(aggregated_patterns)
        aggregated_patterns = prune_shifted_redundant_patterns(
            aggregated_patterns,
            min_length=args.shift_prune_min_len,
            max_length_gap=args.shift_prune_max_gap,
            count_tolerance=args.shift_prune_count_tol,
        )
        save_log(f"{before_shift} -> {len(aggregated_patterns)} patterns after shift pruning")

    compressor_ready = build_compressor_patterns(aggregated_patterns)

    save_log("\nOptimizing and assigning final tokens...")
    final_patterns = optimize_patterns(compressor_ready, args.overhead)
    save_to_json(final_patterns, args.output)

    save_log(f"\n{sep}")
    save_log("ANALYSIS COMPLETE")
    save_log(f"Final patterns: {len(final_patterns)}")
    save_log(f"Saved to: {args.output}")
    save_log(f"Log saved to: {_LOG_DIR}/{_LOG_FILE}")
    save_log(f"{sep}")


if __name__ == "__main__":
    main()
