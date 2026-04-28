import random
import sys
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

def generate_quality_string(length):
    """Generates a random quality string"""
    quality_chars = "IIIIIIIIIHHHHHHHHHJJJJJJJJJ"
    return ''.join(random.choice(quality_chars) for _ in range(length))

def generate_sequence(seq_id, target_length=100):
    """
    Generates a FASTQ sequence with repeated patterns
    """
    sequence = ""
    
    patterns_to_use = random.randint(3, 6)
    for _ in range(patterns_to_use):
        pattern = random.choice(COMMON_PATTERNS)
        sequence += pattern
    
    while len(sequence) < target_length:
        sequence += random.choice("ATCG")
    
    sequence = sequence[:target_length]
    
    quality = generate_quality_string(len(sequence))
    
    return f"@SEQ_ID_{seq_id}\n{sequence}\n+\n{quality}\n"

def generate_fastq_file(output_file, num_sequences=1000, seq_length=100):
    """
    Generates a synthetic FASTQ file
    
    Args:
        output_file: Output file name
        num_sequences: Number of sequences to generate
        seq_length: Length of each sequence
    """
    print("Generating synthetic FASTQ file...")
    print(f"   Sequences: {num_sequences}")
    print(f"   Length: {seq_length} nucleotides")
    
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, 'w') as f:
        for i in range(1, num_sequences + 1):
            fastq_entry = generate_sequence(i, seq_length)
            f.write(fastq_entry)
            
            if i % 100 == 0:
                print(f"   Generated {i}/{num_sequences} sequences...")
    
    print(f"\nFile generated: {output_file}")
    
    import os
    file_size = os.path.getsize(output_file)
    print(f"   Size: {file_size:,} bytes ({file_size/1024:.1f} KB)")

def main():
    if len(sys.argv) < 2:
        print("Usage: python generate_fastq.py <num_sequences> [output.fastq]")
        print("\nExamples:")
        print("  python generate_fastq.py 1000")
        print("  python generate_fastq.py 5000 large_data.fastq")
        print("  python generate_fastq.py 10000 xl_data.fastq")
        sys.exit(1)
    
    num_sequences = int(sys.argv[1])
    output_file = sys.argv[2] if len(sys.argv) > 2 else f"synthetic_{num_sequences}.fastq"
    
    generate_fastq_file(output_file, num_sequences, seq_length=100)

if __name__ == "__main__":
    main()
