"""FASTQ extraction and synthetic FASTQ generation helpers."""

import gzip
import os
import random
from pathlib import Path


COMMON_PATTERNS = [
    "ATCGATCGATCGATCG",
    "GGCCTTAAAGGCCTTAAA",
    "CCCCGGGGCCCCGGGG",
    "NNNATCGATCGATCGATCG",
    "AAAATTTTAAAATTTT",
    "GCTAGCTAGCTAGCTA",
    "TATATATATATATATA",
    "CGCGCGCGCGCGCGCG",
]


def extract_sequences(fastq_file, output_file, n):
    input_path = Path(fastq_file)
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    opener = gzip.open if input_path.suffix == ".gz" else open

    with opener(input_path, "rt", encoding="utf-8") as input_handle, open(
        output_path, "w", encoding="utf-8"
    ) as output_handle:
        count = 0
        while count < n:
            header = input_handle.readline()
            if not header:
                break
            sequence = input_handle.readline().strip()
            input_handle.readline()
            input_handle.readline()
            output_handle.write(sequence + "\n")
            count += 1
    print(f"Process Ended: Extracted {count} sequences to {output_file}")


def generate_quality_string(length):
    quality_chars = "IIIIIIIIIHHHHHHHHHJJJJJJJJJ"
    return "".join(random.choice(quality_chars) for _ in range(length))


def generate_sequence(seq_id, target_length=100):
    sequence = ""

    patterns_to_use = random.randint(3, 6)
    for _ in range(patterns_to_use):
        sequence += random.choice(COMMON_PATTERNS)

    while len(sequence) < target_length:
        sequence += random.choice("ATCG")

    sequence = sequence[:target_length]
    quality = generate_quality_string(len(sequence))
    return f"@SEQ_ID_{seq_id}\n{sequence}\n+\n{quality}\n"


def generate_fastq_file(output_file, num_sequences=1000, seq_length=100):
    print("Generating synthetic FASTQ file...")
    print(f"   Sequences: {num_sequences}")
    print(f"   Length: {seq_length} nucleotides")

    Path(output_file).parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, "w") as output_handle:
        for sequence_number in range(1, num_sequences + 1):
            fastq_entry = generate_sequence(sequence_number, seq_length)
            output_handle.write(fastq_entry)

            if sequence_number % 100 == 0:
                print(f"   Generated {sequence_number}/{num_sequences} sequences...")

    print(f"\nFile generated: {output_file}")
    file_size = os.path.getsize(output_file)
    print(f"   Size: {file_size:,} bytes ({file_size / 1024:.1f} KB)")