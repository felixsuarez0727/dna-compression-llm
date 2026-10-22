"""Lossless token substitution and expansion for sequence files."""

from pathlib import Path

from .logging_utils import save_log
from .pattern_store import TOKEN_CLOSE, TOKEN_OPEN, load_sorted_patterns, load_token_map


def load_compression_patterns(pattern_file, log_file):
    sorted_patterns = load_sorted_patterns(pattern_file)
    save_log(f"Loaded {len(sorted_patterns)} patterns", filename=log_file)
    return sorted_patterns


def load_decompression_patterns(pattern_file, log_file=None):
    token_map = load_token_map(pattern_file)
    save_log(f"Loaded {len(token_map)} patterns", filename=log_file)
    return token_map


def tokenize_sequence(sequence, sorted_patterns):
    segments = [(False, sequence)]
    for token, data in sorted_patterns:
        pattern = data["sequence"]
        marker = f"{TOKEN_OPEN}{token}{TOKEN_CLOSE}"
        updated = []
        for is_marker, text in segments:
            if is_marker or pattern not in text:
                updated.append((is_marker, text))
                continue
            parts = text.split(pattern)
            for index, part in enumerate(parts):
                if index:
                    updated.append((True, marker))
                if part:
                    updated.append((False, part))
        segments = updated
    return "".join(text for _, text in segments)


def expand_sequence(tokenized_seq, token_map):
    parts = []
    index = 0
    length = len(tokenized_seq)

    while index < length:
        if tokenized_seq[index] == TOKEN_OPEN:
            token_end = tokenized_seq.find(TOKEN_CLOSE, index)
            if token_end == -1:
                raise ValueError("Malformed token")

            token = tokenized_seq[index + 1:token_end]
            if token not in token_map:
                raise ValueError(f"Unknown token: {token}")

            parts.append(token_map[token])
            index = token_end + 1
        else:
            next_token = tokenized_seq.find(TOKEN_OPEN, index)
            if next_token == -1:
                next_token = length
            parts.append(tokenized_seq[index:next_token])
            index = next_token

    return "".join(parts)


def line_ending(line):
    if line.endswith("\r\n"):
        return "\r\n"
    if line.endswith(("\n", "\r")):
        return line[-1]
    return ""


def split_line(line):
    ending = line_ending(line)
    content = line[: len(line) - len(ending)] if ending else line
    return content, ending


def compress_seq(input_file, pattern_file, output_file, log_file):
    sorted_patterns = load_compression_patterns(pattern_file, log_file)

    save_log(f"Compressing: {input_file}", filename=log_file)
    total = 0
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)

    with open(input_file, "r", encoding="utf-8", newline="") as input_handle, open(
        output_file, "w", encoding="utf-8", newline=""
    ) as output_handle:
        for line_number, line in enumerate(input_handle, start=1):
            content, ending = split_line(line)
            if TOKEN_OPEN in content or TOKEN_CLOSE in content:
                raise ValueError(
                    f"{input_file}, line {line_number}: input contains the reserved "
                    f"characters '{TOKEN_OPEN}' or '{TOKEN_CLOSE}' and cannot be "
                    "compressed losslessly"
                )

            output_handle.write(tokenize_sequence(content, sorted_patterns) + ending)
            if content.strip():
                total += 1

    original_size = Path(input_file).stat().st_size
    compressed_size = Path(output_file).stat().st_size
    ratio = 100 * (1 - compressed_size / original_size) if original_size > 0 else 0

    save_log("\nCompression completed", filename=log_file)
    save_log(f"Sequences processed : {total}", filename=log_file)
    save_log(f"Compression ratio   : {ratio:.2f}%", filename=log_file)
    save_log(f"Output file         : {output_file}", filename=log_file)


def decompress_seq(input_file, pattern_file, output_file, log_file):
    token_map = load_decompression_patterns(pattern_file, log_file)

    save_log(f"Decompressing: {input_file}", filename=log_file)
    total = 0
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)

    with open(input_file, "r", encoding="utf-8", newline="") as input_handle, open(
        output_file, "w", encoding="utf-8", newline=""
    ) as output_handle:
        for line in input_handle:
            content, ending = split_line(line)
            output_handle.write(expand_sequence(content, token_map) + ending)
            if content.strip():
                total += 1

    save_log("\nDecompression completed", filename=log_file)
    save_log(f"Sequences restored : {total}", filename=log_file)
    save_log(f"Output file        : {output_file}", filename=log_file)
