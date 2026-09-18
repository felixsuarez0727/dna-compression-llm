"""Lazy DNABERT-2 and HyenaDNA model helpers."""

import hashlib
import re
import sys
from collections import Counter, defaultdict


def _place_model_on_available_device(model, torch, save_log, model_name):
    if torch.cuda.is_available():
        try:
            model = model.to("cuda")
        except (AssertionError, RuntimeError) as error:
            save_log(f"  {model_name}: CUDA placement failed ({error}); falling back to CPU.")
        else:
            save_log(f"  {model_name} loaded on CUDA.")
            return model, "cuda"

    model = model.to("cpu")
    save_log(f"  {model_name} loaded on CPU.")
    return model, "cpu"


def load_dnabert2(model_id, save_log):
    try:
        from transformers import AutoModel, AutoTokenizer
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

    model, device = _place_model_on_available_device(
        model, torch, save_log, "DNABERT-2"
    )

    return {"tokenizer": tokenizer, "model": model, "device": device, "torch": torch}


def _dnabert2_tokenize_sequences(client, sequences):
    tokenizer = client["tokenizer"]
    special_tokens = {"[CLS]", "[SEP]", "[PAD]", "[MASK]", "[UNK]"}
    all_tokens = []
    for sequence in sequences:
        ids = tokenizer(sequence, return_tensors="pt", padding=False)["input_ids"][0]
        tokens = tokenizer.convert_ids_to_tokens(ids)
        dna_tokens = [token for token in tokens if token not in special_tokens]
        all_tokens.append(dna_tokens)
    return all_tokens


def _find_tandem_runs_in_tokens(token_list, min_repeats=2):
    runs = []
    index = 0
    while index < len(token_list):
        end = index + 1
        while end < len(token_list) and token_list[end] == token_list[index]:
            end += 1
        run_length = end - index
        if run_length >= min_repeats:
            runs.append((token_list[index].upper(), run_length, index))
        index = end
    return runs


def _build_attention_patterns(client, sequences, save_log, max_seq_len=256):
    tokenizer = client["tokenizer"]
    model = client["model"]
    device = client["device"]
    torch = client["torch"]

    pattern_counter = Counter()
    special_tokens = {"[CLS]", "[SEP]", "[PAD]", "[MASK]"}

    for sequence in sequences:
        sequence = sequence[:max_seq_len]
        inputs = tokenizer(sequence, return_tensors="pt", padding=True).to(device)
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
            dna_tokens = [
                (index, token)
                for index, token in enumerate(tokens)
                if token not in special_tokens
            ]

            if attentions is not None:
                attention_layers = torch.stack(list(attentions))
                average_attention = attention_layers.mean(dim=(0, 1, 2))
                column_sum = average_attention.sum(dim=0)
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
                column_sum = hidden[0].norm(dim=-1)

            top_k = min(10, len(dna_tokens))
            if top_k == 0:
                continue

            _, top_indices = torch.topk(column_sum[:len(tokens)], top_k)
            top_set = set(top_indices.tolist())

            for index, (position, token) in enumerate(dna_tokens):
                if position in top_set:
                    for window in range(1, 4):
                        end = index + window
                        if end <= len(dna_tokens):
                            candidate = "".join(
                                item for _, item in dna_tokens[index:end]
                            ).upper()
                            if len(candidate) >= 4 and re.fullmatch(r"[ACGTN]+", candidate):
                                pattern_counter[candidate] += 1

        except Exception as error:
            save_log(f"  [DNABERT-2] Atencion no disponible para una secuencia: {error}")
            continue

    return pattern_counter


def analyze_batch_dnabert2(client, sequences, batch_num, save_log):
    save_log(f"  [DNABERT-2] Batch #{batch_num}: tokenizing {len(sequences)} sequences...")

    all_token_lists = _dnabert2_tokenize_sequences(client, sequences)
    pattern_counter = Counter()

    for token_list in all_token_lists:
        for token in token_list:
            token_upper = token.upper()
            if len(token_upper) >= 4 and re.fullmatch(r"[ACGTN]+", token_upper):
                pattern_counter[token_upper] += 1

        runs = _find_tandem_runs_in_tokens(token_list, min_repeats=2)
        for unit, repeats, _ in runs:
            if len(unit) >= 4 and re.fullmatch(r"[ACGTN]+", unit):
                for repeat_count in range(2, repeats + 1):
                    pattern_counter[unit * repeat_count] += 1
                pattern_counter[unit] += repeats

    save_log(f"  [DNABERT-2] Batch #{batch_num}: extracting attention motifs...")
    attention_patterns = _build_attention_patterns(client, sequences, save_log)
    pattern_counter.update(attention_patterns)

    lines = [
        f"{pattern}:{count}"
        for pattern, count in pattern_counter.most_common(100)
        if len(pattern) >= 4 and re.fullmatch(r"[ACGTN]+", pattern)
    ]
    return "\n".join(lines) if lines else None


def analyze_synthesis_dnabert2(client, all_candidates, top_kmers):
    tokenizer = client["tokenizer"]
    special_tokens = {"[CLS]", "[SEP]", "[PAD]", "[MASK]"}

    boosted = {}
    for sequence, count in all_candidates.items():
        ids = tokenizer(sequence, return_tensors="pt", padding=False)["input_ids"][0]
        tokens = [
            token for token in tokenizer.convert_ids_to_tokens(ids) if token not in special_tokens
        ]
        boost = 1.5 if len(tokens) == 1 else 1.0
        boosted[sequence] = int(count * boost)

    lines = [
        f"{pattern}:{count}"
        for pattern, count in sorted(boosted.items(), key=lambda item: -item[1])[:80]
        if re.fullmatch(r"[ACGTN]+", pattern)
    ]
    return "\n".join(lines)


def load_hyenadna(model_id, save_log):
    try:
        from transformers import AutoModel, AutoTokenizer
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

    model, device = _place_model_on_available_device(
        model, torch, save_log, "HyenaDNA"
    )

    return {"tokenizer": tokenizer, "model": model, "device": device, "torch": torch}


def _mine_long_repeats(
    sequences, save_log, min_len=24, max_len=512, stride=4, min_occurrences=3, max_candidates=20000
):
    save_log("  [HyenaDNA] Mining long repeated substrings...")
    candidates = defaultdict(int)

    lengths = [24, 32, 48, 64, 96, 128, 160, 192, 256, 384, 512]
    lengths = [length for length in lengths if length <= max_len]

    for sequence in sequences:
        sequence_length = len(sequence)
        for length in lengths:
            if sequence_length < length:
                continue
            for index in range(0, sequence_length - length + 1, stride):
                chunk = sequence[index:index + length]
                if "N" in chunk and chunk.count("N") > (length // 8):
                    continue
                candidates[chunk] += 1

        for unit_length in range(4, 65):
            for start in range(0, sequence_length - unit_length):
                unit = sequence[start:start + unit_length]
                repeats = 1
                position = start + unit_length
                while position + unit_length <= sequence_length:
                    if sequence[position:position + unit_length] == unit:
                        repeats += 1
                        position += unit_length
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
        key=lambda item: (item[1]["gain"], item[1]["length"], item[1]["count"]),
        reverse=True,
    )[:max_candidates]

    save_log(f"  [HyenaDNA] Long-repeat candidates retained: {len(sorted_patterns)}")
    return sorted_patterns


def _hyenadna_embeddings(client, sequences, max_length=1024):
    tokenizer = client["tokenizer"]
    model = client["model"]
    device = client["device"]
    torch = client["torch"]

    embeddings = []
    for sequence in sequences:
        sequence = sequence[:max_length]
        inputs = tokenizer(
            sequence, return_tensors="pt", truncation=True, max_length=max_length
        ).to(device)
        with torch.no_grad():
            outputs = model(**inputs)
        hidden = outputs.last_hidden_state[0]
        embeddings.append(hidden.mean(dim=0).cpu().numpy())

    return embeddings


def _deduplicate_compression_patterns(sorted_patterns, save_log):
    save_log("  [HyenaDNA] Deduplicating overlapping patterns...")
    kept = {}
    fingerprints = set()

    for pattern, metadata in sorted_patterns:
        fingerprint = hashlib.md5(pattern[:64].encode()).hexdigest()
        if fingerprint in fingerprints:
            continue
        redundant = any(
            len(existing) > len(pattern) and pattern in existing for existing in kept
        )
        if redundant:
            continue
        kept[pattern] = metadata
        fingerprints.add(fingerprint)

    save_log(f"  [HyenaDNA] Final deduplicated patterns: {len(kept)}")
    return kept


def analyze_batch_hyenadna(client, sequences, batch_num, save_log):
    save_log(f"  [HyenaDNA] Batch #{batch_num}: compression-oriented analysis...")

    try:
        _hyenadna_embeddings(client, sequences[:8])
    except Exception as error:
        save_log(f"  [HyenaDNA] embedding warning: {error}")

    candidates = _mine_long_repeats(
        sequences,
        save_log,
        min_len=24,
        max_len=512,
        stride=8,
        min_occurrences=2,
        max_candidates=5000,
    )
    deduplicated = _deduplicate_compression_patterns(candidates, save_log)

    lines = []
    for pattern, metadata in sorted(
        deduplicated.items(), key=lambda item: (item[1]["gain"], item[1]["length"]), reverse=True
    )[:200]:
        lines.append(f"{pattern}:{metadata['count']}")

    save_log(f"  [HyenaDNA] Batch #{batch_num}: {len(lines)} high-value patterns")
    return "\n".join(lines)


def analyze_synthesis_hyenadna(client, all_candidates, top_kmers, find_tandem_unit, save_log):
    save_log("  [HyenaDNA] Compression-oriented synthesis...")

    rescored = {}
    for sequence, count in all_candidates.items():
        gain = (len(sequence) - 1) * count
        tandem_bonus = 1.0
        unit, repeats = find_tandem_unit(sequence)
        if unit:
            tandem_bonus += min(repeats * 0.25, 4.0)
        length_bonus = 1.0 + (len(sequence) / 128)
        rescored[sequence] = int(gain * tandem_bonus * length_bonus)

    lines = []
    for sequence, score in sorted(
        rescored.items(), key=lambda item: (item[1], len(item[0])), reverse=True
    )[:400]:
        if len(sequence) >= 16:
            lines.append(f"{sequence}:{score}")

    save_log(f"  [HyenaDNA] synthesis retained {len(lines)} patterns")
    return "\n".join(lines)