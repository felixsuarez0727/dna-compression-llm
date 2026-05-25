import argparse
import time
import sys
import json
import re
import os
import hashlib
from pathlib import Path
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

MICRO_TOKEN_POOL = list("bdefhijklmoprsvwxyz")

_LOG_FILE = None


def save_log(line):
    safe_line = line.replace("→", "->")
    print(safe_line)
    if _LOG_FILE:
        os.makedirs("data/logs", exist_ok=True)
        with open(f"data/logs/{_LOG_FILE}", "a", encoding="utf-8") as f:
            f.write(safe_line + "\n")

PROVIDER_DEFAULTS = {
    "deepseek": "deepseek-chat",
    "chatgpt": "gpt-4o-mini",
    "gemini": "gemini-2.5-flash-lite",
    "dnabert2": "zhihan1996/DNABERT-2-117M",
    "hyenadna": "LongSafari/hyenadna-large-1m-seqlen-hf",
}




# ─────────────────────────────────────────────────────────────
# DNABERT-2 helpers
# ─────────────────────────────────────────────────────────────

def _load_dnabert2(model_id):
    try:
        from transformers import AutoTokenizer, AutoModel
        import torch
    except ImportError:
        save_log(
            "Error: 'transformers' and 'torch' are required for DNABERT-2.\n"
            "  pip install transformers torch"
        )
        sys.exit(1)

    save_log(f"  Loading DNABERT-2 from '{model_id}' (may download on first run)...")
    tokenizer = AutoTokenizer.from_pretrained(model_id, trust_remote_code=True)
    model = AutoModel.from_pretrained(model_id, trust_remote_code=True)
    model.eval()

    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = model.to(device)
    save_log(f"  DNABERT-2 loaded on {device.upper()}.")

    return {"tokenizer": tokenizer, "model": model, "device": device, "torch": torch}


def _dnabert2_tokenize_sequences(client, sequences):
    tokenizer = client["tokenizer"]
    SPECIAL = {"[CLS]", "[SEP]", "[PAD]", "[MASK]", "[UNK]"}
    all_tokens = []
    for seq in sequences:
        ids = tokenizer(seq, return_tensors="pt", padding=False)["input_ids"][0]
        tokens = tokenizer.convert_ids_to_tokens(ids)
        dna_tokens = [t for t in tokens if t not in SPECIAL]
        all_tokens.append(dna_tokens)
    return all_tokens


def _find_tandem_runs_in_tokens(token_list, min_repeats=2):
    runs = []
    i = 0
    while i < len(token_list):
        j = i + 1
        while j < len(token_list) and token_list[j] == token_list[i]:
            j += 1
        run_len = j - i
        if run_len >= min_repeats:
            runs.append((token_list[i].upper(), run_len, i))
        i = j
    return runs


def _build_attention_patterns(client, sequences, max_seq_len=256):
    tokenizer = client["tokenizer"]
    model     = client["model"]
    device    = client["device"]
    torch     = client["torch"]

    pattern_counter = Counter()
    SPECIAL = {"[CLS]", "[SEP]", "[PAD]", "[MASK]"}

    for seq in sequences:
        seq = seq[:max_seq_len]
        inputs = tokenizer(seq, return_tensors="pt", padding=True).to(device)
        try:
            with torch.no_grad():
                outputs = model(**inputs, output_attentions=True)

            attentions = None
            if hasattr(outputs, "attentions") and outputs.attentions is not None:
                attentions = outputs.attentions
            elif isinstance(outputs, tuple):
                for item in outputs:
                    if (
                        isinstance(item, tuple)
                        and len(item) > 0
                        and isinstance(item[0], torch.Tensor)
                        and item[0].dim() == 4
                    ):
                        attentions = item
                        break

            tokens = tokenizer.convert_ids_to_tokens(inputs["input_ids"][0])
            dna_tokens = [(i, t) for i, t in enumerate(tokens) if t not in SPECIAL]

            if attentions is not None:
                attn_layers = torch.stack(list(attentions))
                avg_attn = attn_layers.mean(dim=(0, 1, 2))
                col_sum = avg_attn.sum(dim=0)
            else:
                hidden = None
                if isinstance(outputs, tuple):
                    for item in outputs:
                        if isinstance(item, torch.Tensor) and item.dim() == 3:
                            hidden = item
                            break
                elif hasattr(outputs, "last_hidden_state"):
                    hidden = outputs.last_hidden_state
                if hidden is None:
                    continue
                col_sum = hidden[0].norm(dim=-1)

            top_k = min(10, len(dna_tokens))
            if top_k == 0:
                continue

            _, top_indices = torch.topk(col_sum[:len(tokens)], top_k)
            top_set = set(top_indices.tolist())

            for idx, (pos, tok) in enumerate(dna_tokens):
                if pos in top_set:
                    for window in range(1, 4):
                        end = idx + window
                        if end <= len(dna_tokens):
                            candidate = "".join(t for _, t in dna_tokens[idx:end]).upper()
                            if len(candidate) >= 4 and re.fullmatch(r"[ACGTN]+", candidate):
                                pattern_counter[candidate] += 1

        except Exception as e:
            save_log(f"  [DNABERT-2] Atencion no disponible para una secuencia: {e}")
            continue

    return pattern_counter


def analyze_batch_dnabert2(client, sequences, batch_num):
    save_log(f"  [DNABERT-2] Batch #{batch_num}: tokenizing {len(sequences)} sequences...")

    all_token_lists = _dnabert2_tokenize_sequences(client, sequences)
    pattern_counter = Counter()

    for token_list in all_token_lists:
        for tok in token_list:
            tok_upper = tok.upper()
            if len(tok_upper) >= 4 and re.fullmatch(r"[ACGTN]+", tok_upper):
                pattern_counter[tok_upper] += 1

        runs = _find_tandem_runs_in_tokens(token_list, min_repeats=2)
        for unit, reps, _ in runs:
            if len(unit) >= 4 and re.fullmatch(r"[ACGTN]+", unit):
                for r in range(2, reps + 1):
                    pattern_counter[unit * r] += 1
                pattern_counter[unit] += reps

    save_log(f"  [DNABERT-2] Batch #{batch_num}: extracting attention motifs...")
    attn_patterns = _build_attention_patterns(client, sequences)
    pattern_counter.update(attn_patterns)

    lines = [
        f"{pat}:{cnt}"
        for pat, cnt in pattern_counter.most_common(100)
        if len(pat) >= 4 and re.fullmatch(r"[ACGTN]+", pat)
    ]
    return "\n".join(lines) if lines else None


def analyze_synthesis_dnabert2(client, all_candidates, top_kmers):
    tokenizer = client["tokenizer"]
    SPECIAL = {"[CLS]", "[SEP]", "[PAD]", "[MASK]"}

    boosted = {}
    for seq, count in all_candidates.items():
        ids = tokenizer(seq, return_tensors="pt", padding=False)["input_ids"][0]
        tokens = [t for t in tokenizer.convert_ids_to_tokens(ids) if t not in SPECIAL]
        boost = 1.5 if len(tokens) == 1 else 1.0
        boosted[seq] = int(count * boost)

    lines = [
        f"{pat}:{cnt}"
        for pat, cnt in sorted(boosted.items(), key=lambda x: -x[1])[:80]
        if re.fullmatch(r"[ACGTN]+", pat)
    ]
    return "\n".join(lines)


# ─────────────────────────────────────────────────────────────
# HyenaDNA helpers
# ─────────────────────────────────────────────────────────────

def _load_hyenadna(model_id):
    try:
        from transformers import AutoTokenizer, AutoModel
        import torch
    except ImportError:
        save_log(
            "Error: 'transformers' and 'torch' are required for HyenaDNA.\n"
            "  pip install transformers torch"
        )
        sys.exit(1)

    save_log(f"  Loading HyenaDNA from '{model_id}' (may download on first run)...")
    tokenizer = AutoTokenizer.from_pretrained(model_id, trust_remote_code=True)
    model = AutoModel.from_pretrained(model_id, trust_remote_code=True)
    model.eval()

    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = model.to(device)
    save_log(f"  HyenaDNA loaded on {device.upper()}.")

    return {"tokenizer": tokenizer, "model": model, "device": device, "torch": torch}


def _mine_long_repeats(sequences, min_len=24, max_len=512, stride=4,
                       min_occurrences=3, max_candidates=20000):
    save_log("  [HyenaDNA] Mining long repeated substrings...")
    candidates = defaultdict(int)

    lengths = [24, 32, 48, 64, 96, 128, 160, 192, 256, 384, 512]
    lengths = [x for x in lengths if x <= max_len]

    for seq in sequences:
        seq_len = len(seq)
        for k in lengths:
            if seq_len < k:
                continue
            for i in range(0, seq_len - k + 1, stride):
                chunk = seq[i:i + k]
                if "N" in chunk and chunk.count("N") > (k // 8):
                    continue
                candidates[chunk] += 1

        for unit_len in range(4, 65):
            for start in range(0, seq_len - unit_len):
                unit = seq[start:start + unit_len]
                repeats = 1
                pos = start + unit_len
                while pos + unit_len <= seq_len:
                    if seq[pos:pos + unit_len] == unit:
                        repeats += 1
                        pos += unit_len
                    else:
                        break
                if repeats >= 3:
                    repeated = unit * repeats
                    if len(repeated) >= min_len:
                        candidates[repeated] += repeats * 3

    filtered = {}
    for pattern, count in candidates.items():
        if count < min_occurrences:
            continue
        gain = (len(pattern) - 1) * count
        if gain <= len(pattern):
            continue
        filtered[pattern] = {"count": count, "gain": gain, "length": len(pattern)}

    sorted_patterns = sorted(
        filtered.items(),
        key=lambda x: (x[1]["gain"], x[1]["length"], x[1]["count"]),
        reverse=True,
    )[:max_candidates]

    save_log(f"  [HyenaDNA] Long-repeat candidates retained: {len(sorted_patterns)}")
    return sorted_patterns


def _hyenadna_embeddings(client, sequences, max_length=1024):
    tokenizer = client["tokenizer"]
    model     = client["model"]
    device    = client["device"]
    torch     = client["torch"]

    embeddings = []
    for seq in sequences:
        seq = seq[:max_length]
        inputs = tokenizer(
            seq, return_tensors="pt", truncation=True, max_length=max_length
        ).to(device)
        with torch.no_grad():
            outputs = model(**inputs)
        hidden = outputs.last_hidden_state[0]
        embeddings.append(hidden.mean(dim=0).cpu().numpy())

    return embeddings


def _deduplicate_compression_patterns(sorted_patterns):
    save_log("  [HyenaDNA] Deduplicating overlapping patterns...")
    kept = {}
    fingerprints = set()

    for pattern, meta in sorted_patterns:
        fp = hashlib.md5(pattern[:64].encode()).hexdigest()
        if fp in fingerprints:
            continue
        redundant = any(
            len(existing) > len(pattern) and pattern in existing
            for existing in kept
        )
        if redundant:
            continue
        kept[pattern] = meta
        fingerprints.add(fp)

    save_log(f"  [HyenaDNA] Final deduplicated patterns: {len(kept)}")
    return kept


def analyze_batch_hyenadna(client, sequences, batch_num):
    save_log(f"  [HyenaDNA] Batch #{batch_num}: compression-oriented analysis...")

    try:
        _ = _hyenadna_embeddings(client, sequences[:8])
    except Exception as e:
        save_log(f"  [HyenaDNA] embedding warning: {e}")

    candidates = _mine_long_repeats(
        sequences, min_len=24, max_len=512, stride=8,
        min_occurrences=2, max_candidates=5000,
    )
    dedup = _deduplicate_compression_patterns(candidates)

    lines = []
    for pattern, meta in sorted(
        dedup.items(), key=lambda x: (x[1]["gain"], x[1]["length"]), reverse=True
    )[:200]:
        lines.append(f"{pattern}:{meta['count']}")

    save_log(f"  [HyenaDNA] Batch #{batch_num}: {len(lines)} high-value patterns")
    return "\n".join(lines)


def analyze_synthesis_hyenadna(client, all_candidates, top_kmers):
    save_log("  [HyenaDNA] Compression-oriented synthesis...")

    rescored = {}
    for seq, count in all_candidates.items():
        gain = (len(seq) - 1) * count
        tandem_bonus = 1.0
        unit, reps = find_tandem_unit(seq)
        if unit:
            tandem_bonus += min(reps * 0.25, 4.0)
        length_bonus = 1.0 + (len(seq) / 128)
        rescored[seq] = int(gain * tandem_bonus * length_bonus)

    lines = []
    for seq, score in sorted(rescored.items(), key=lambda x: (x[1], len(x[0])), reverse=True)[:400]:
        if len(seq) >= 16:
            lines.append(f"{seq}:{score}")

    save_log(f"  [HyenaDNA] synthesis retained {len(lines)} patterns")
    return "\n".join(lines)



def build_client(provider, api_key):
    if provider in ("deepseek", "chatgpt"):
        from openai import OpenAI
        if provider == "deepseek":
            return OpenAI(api_key=api_key, base_url="https://api.deepseek.com")
        return OpenAI(api_key=api_key)
    if provider == "gemini":
        from google import genai
        return genai.Client(api_key=api_key)
    if provider == "dnabert2":
        return _load_dnabert2(PROVIDER_DEFAULTS["dnabert2"])
    if provider == "hyenadna":
        return _load_hyenadna(PROVIDER_DEFAULTS["hyenadna"])
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
    except Exception as e:
        save_log(f"\n  Warning: Error in batch {batch_num}: {e}")
        return None


def analyze_batch_synthesis(client, provider, model_id, all_candidates, top_kmers):
    if provider == "dnabert2":
        return analyze_synthesis_dnabert2(client, all_candidates, top_kmers)
    if provider == "hyenadna":
        return analyze_synthesis_hyenadna(client, all_candidates, top_kmers)

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

    try:
        return call_llm(client, provider, model_id, "You are an expert DNA compression AI.", user_prompt)
    except Exception as e:
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
    parser.add_argument("--key", "-k", default="local", help="API key for the selected provider (not required for dnabert2/hyenadna)")
    parser.add_argument(
        "--provider", "-p",
        required=True,
        choices=["deepseek", "chatgpt", "gemini", "dnabert2", "hyenadna"],
        help="LLM/DNA model provider to use",
    )
    parser.add_argument("--output", "-o", default="patterns.json", help="Output JSON file")
    parser.add_argument("--batch_size", "-b", type=int, default=30, help="Sequences per API call")
    parser.add_argument("--model", "-m", default=None, help="Model ID (defaults to provider default)")
    parser.add_argument("--max_batches", type=int, default=0, help="Max number of LLM batches")
    parser.add_argument("--threads", type=int, default=4, help="Parallel threads for batch analysis (chatgpt and deepseek only, default: 4)")
    parser.add_argument("--overhead", type=int, default=5, help="Dictionary entry cost in bytes for final optimization (default: 5)")
    parser.add_argument("--no_validate", action="store_true", help="Disable validation against real sequences")
    parser.add_argument("--no_expand", action="store_true", help="Disable tandem repeat expansion")

    args = parser.parse_args()

    model_id = args.model or PROVIDER_DEFAULTS[args.provider]

    if args.model and args.provider in ("dnabert2", "hyenadna"):
        PROVIDER_DEFAULTS[args.provider] = args.model

    client = build_client(args.provider, args.key)

    ts = time.strftime("%Y%m%d_%H%M%S")
    global _LOG_FILE
    _LOG_FILE = f"{Path(args.output).stem}_log_{ts}.txt"

    sep = "=" * 56
    save_log(f"\n{sep}")
    save_log(f"  FASTQ Pattern Detector  [{args.provider}]")
    save_log(f"{sep}")
    save_log(f"  File       : {args.file}")
    save_log(f"  Model      : {model_id}")
    save_log(f"  Batch size : {args.batch_size}")
    save_log(f"  Overhead   : {args.overhead}")
    if args.provider in ("chatgpt", "deepseek"):
        save_log(f"  Threads    : {args.threads}")
    save_log(f"  Log file   : data/logs/{_LOG_FILE}")
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

    if args.provider in ("dnabert2", "hyenadna"):
        # Local GPU/CPU-bound models: always sequential
        for batch_num, batch in all_batches:
            raw = analyze_batch_with_context(client, args.provider, model_id, batch, batch_num, top_kmers)
            if raw:
                parsed = parse_model_response(raw)
                save_log(f"  Batch #{batch_num} -> {len(parsed)} patterns detected")
                for pat, count in parsed.items():
                    aggregated_patterns[pat] = aggregated_patterns.get(pat, 0) + count
            else:
                save_log(f"  Batch #{batch_num} -> no result")

    elif args.provider in ("chatgpt", "deepseek") and args.threads > 1:
        save_log(f"  Using {args.threads} parallel threads for {len(all_batches)} batches...")

        def run_batch(batch_num, batch):
            thread_client = build_client(args.provider, args.key)
            raw = analyze_batch_with_context(thread_client, args.provider, model_id, batch, batch_num, top_kmers)
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
            raw = analyze_batch_with_context(client, args.provider, model_id, batch, batch_num, top_kmers)
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
        raw_syn = analyze_batch_synthesis(client, args.provider, model_id, aggregated_patterns, top_kmers)
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

    compressor_ready = build_compressor_patterns(aggregated_patterns)

    save_log("\nOptimizing and assigning final tokens...")
    final_patterns = optimize_patterns(compressor_ready, args.overhead)
    save_to_json(final_patterns, args.output)

    save_log(f"\n{sep}")
    save_log("ANALYSIS COMPLETE")
    save_log(f"Final patterns: {len(final_patterns)}")
    save_log(f"Saved to: {args.output}")
    save_log(f"Log saved to: data/logs/{_LOG_FILE}")
    save_log(f"{sep}")


if __name__ == "__main__":
    main()