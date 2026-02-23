import json
import argparse
import os
from pathlib import Path
import time

def load_patterns(pattern_file):
    with open(pattern_file, "r", encoding="utf-8") as f:
        patterns = json.load(f)

    token_map = {
        token: data["sequence"]
        for token, data in patterns.items()
    }

    print(f"Loaded {len(token_map)} patterns")
    return token_map

def expand_sequence(tokenized_seq, token_map):
    result = ""
    i = 0
    n = len(tokenized_seq)

    while i < n:
        if tokenized_seq[i] == "<":
            j = tokenized_seq.find(">", i)
            if j == -1:
                raise ValueError("Malformed token")

            token = tokenized_seq[i+1:j]

            if token not in token_map:
                raise ValueError(f"Unknown token: {token}")

            result += token_map[token]
            i = j + 1
        else:
            result += tokenized_seq[i]
            i += 1

    return result

def decompress_seq(input_file, pattern_file, output_file, log_file):
    token_map = load_patterns(pattern_file)

    save_log(f"Decompressing: {input_file}", filename=log_file)

    total = 0

    with open(input_file, "r", encoding="utf-8") as inp, \
         open(output_file, "w", encoding="utf-8") as out:

        for line in inp:
            tokenized = line.strip()
            if not tokenized:
                continue

            restored = expand_sequence(tokenized, token_map)
            out.write(restored + "\n")
            total += 1

    save_log(f"\nDecompression completed", filename=log_file)
    save_log(f"Sequences restored : {total}", filename=log_file)
    save_log(f"Output file        : {output_file}", filename=log_file)

def save_log(line, filename="compression_stats.log"):
    safe_line = line.replace("→", "->")
    os.makedirs(f"data/", exist_ok=True)
    with open(f"data/{filename}", "a", encoding="utf-8") as log_file:
        print(safe_line)
        log_file.write(safe_line + "\n")

def main():
    parser = argparse.ArgumentParser(description="Sequence file decompressor")
    parser.add_argument("-i", "--input", required=True,
                        help="Input compressed sequence file")
    parser.add_argument("-p", "--patterns", required=True,
                        help="Pattern JSON file")
    parser.add_argument("-o", "--output", required=True,
                        help="Output restored sequence file")

    args = parser.parse_args()

    ts = time.strftime("%Y%m%d_%H%M%S")
    log_file = f"{Path(args.output).stem}.decompress_log_{ts}.txt"

    decompress_seq(args.input, args.patterns, args.output, log_file)

if __name__ == "__main__":
    main()
