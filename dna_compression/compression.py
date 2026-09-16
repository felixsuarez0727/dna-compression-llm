"""Lossless token substitution and expansion for sequence files."""

from pathlib import Path

from .logging_utils import save_log
from .pattern_store import load_sorted_patterns, load_token_map


def load_compression_patterns(pattern_file, log_file):
    sorted_patterns = load_sorted_patterns(pattern_file)
    save_log(f"Loaded {len(sorted_patterns)} patterns", filename=log_file)
    return sorted_patterns


def load_decompression_patterns(pattern_file):
    token_map = load_token_map(pattern_file)
    print(f"Loaded {len(token_map)} patterns")
    return token_map


def tokenize_sequence(sequence, sorted_patterns):
    for token, data in sorted_patterns:
        pattern = data["sequence"]
        sequence = sequence.replace(pattern, f"<{token}>")
    return sequence


def expand_sequence(tokenized_seq, token_map):
    result = ""
    index = 0
    length = len(tokenized_seq)

    while index < length:
        if tokenized_seq[index] == "<":
            token_end = tokenized_seq.find(">", index)
            if token_end == -1:
                raise ValueError("Malformed token")

            token = tokenized_seq[index + 1:token_end]
            if token not in token_map:
                raise ValueError(f"Unknown token: {token}")

            result += token_map[token]
            index = token_end + 1
        else:
            result += tokenized_seq[index]
            index += 1

    return result


def line_ending(line):
    if line.endswith("\r\n"):
        return "\r\n"
    if line.endswith(("\n", "\r")):
        return line[-1]
    return ""


def compress_seq(input_file, pattern_file, output_file, log_file):
    sorted_patterns = load_compression_patterns(pattern_file, log_file)

    save_log(f"Compressing: {input_file}", filename=log_file)
    total = 0
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)

    with open(input_file, "r", encoding="utf-8", newline="") as input_handle, open(
        output_file, "w", encoding="utf-8", newline=""
    ) as output_handle:
        for line in input_handle:
            sequence = line.strip()
            if not sequence:
                continue

            tokenized = tokenize_sequence(sequence, sorted_patterns)
            output_handle.write(tokenized + line_ending(line))
            total += 1

    original_size = Path(input_file).stat().st_size
    compressed_size = Path(output_file).stat().st_size
    ratio = 100 * (1 - compressed_size / original_size) if original_size > 0 else 0

    save_log("\nCompression completed", filename=log_file)
    save_log(f"Sequences processed : {total}", filename=log_file)
    save_log(f"Compression ratio   : {ratio:.2f}%", filename=log_file)
    save_log(f"Output file         : {output_file}", filename=log_file)


def decompress_seq(input_file, pattern_file, output_file, log_file):
    token_map = load_decompression_patterns(pattern_file)

    save_log(f"Decompressing: {input_file}", filename=log_file)
    total = 0
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)

    with open(input_file, "r", encoding="utf-8", newline="") as input_handle, open(
        output_file, "w", encoding="utf-8", newline=""
    ) as output_handle:
        for line in input_handle:
            tokenized = line.strip()
            if not tokenized:
                continue

            restored = expand_sequence(tokenized, token_map)
            output_handle.write(restored + line_ending(line))
            total += 1

    save_log("\nDecompression completed", filename=log_file)
    save_log(f"Sequences restored : {total}", filename=log_file)
    save_log(f"Output file        : {output_file}", filename=log_file)