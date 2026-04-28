import json
import argparse
import os
from pathlib import Path
import time

def load_patterns(pattern_file, log_file):
    with open(pattern_file, 'r', encoding='utf-8') as f:
        patterns = json.load(f)

    sorted_patterns = sorted(
        patterns.items(),
        key=lambda x: x[1]["priority"]
    )

    save_log(f"Loaded {len(patterns)} patterns", filename=log_file)
    return sorted_patterns

def tokenize_sequence(sequence, sorted_patterns):
    for token, data in sorted_patterns:
        pattern = data["sequence"]
        sequence = sequence.replace(pattern, f"<{token}>")
    return sequence

def compress_seq(input_file, pattern_file, output_file, log_file):
    sorted_patterns = load_patterns(pattern_file, log_file)

    save_log(f"Compressing: {input_file}", filename=log_file)
    total = 0
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)

    # Se agrega encoding='utf-8' para soportar los tokens de 1 byte (0-255)
    with open(input_file, "r", encoding='utf-8') as inp, \
         open(output_file, "w", encoding='utf-8') as out:

        for line in inp:
            seq = line.strip()
            if not seq:
                continue

            tokenized = tokenize_sequence(seq, sorted_patterns)
            out.write(tokenized + "\n")
            total += 1

    original_size = Path(input_file).stat().st_size
    compressed_size = Path(output_file).stat().st_size
    
    ratio = 100 * (1 - compressed_size / original_size) if original_size > 0 else 0

    save_log("\nCompression completed", filename=log_file)
    save_log(f"Sequences processed : {total}", filename=log_file)
    save_log(f"Compression ratio   : {ratio:.2f}%", filename=log_file)
    save_log(f"Output file         : {output_file}", filename=log_file)
    
  

def save_log(line, filename="compression_stats.log"):
    safe_line = line.replace("→", "->")
    os.makedirs("data/logs", exist_ok=True)
    with open(f"data/logs/{filename}", "a", encoding="utf-8") as log_file:
        print(safe_line)
        log_file.write(safe_line + "\n")

def main():
    parser = argparse.ArgumentParser(description="Sequence file compressor")
    parser.add_argument("-i", "--input", required=True, help="Input sequence file")
    parser.add_argument("-p", "--patterns", required=True, help="Pattern JSON file")
    parser.add_argument("-o", "--output", required=True, help="Output compressed file")

    args = parser.parse_args()

    ts = time.strftime("%Y%m%d_%H%M%S")
    log_file = f"{Path(args.output).stem}.compress_log_{ts}.txt"

    compress_seq(args.input, args.patterns, args.output, log_file)

if __name__ == "__main__":
    main()
