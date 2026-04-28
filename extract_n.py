import argparse
import gzip
from pathlib import Path

def extract_sequences(fastq_file, output_file, n):
    input_path = Path(fastq_file)
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    opener = gzip.open if input_path.suffix == ".gz" else open

    with opener(input_path, "rt", encoding="utf-8") as f, open(output_path, "w", encoding="utf-8") as out:
        count = 0
        while count < n:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()  # "+" line
            f.readline()  # quality line
            out.write(seq + "\n")
            count += 1
    print(f"Process Ended: Extracted {count} sequences to {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description="Extract N sequences (sequence lines only) from a FASTQ file."
    )
    
    parser.add_argument("--input", required=True, help="Input FASTQ file")
    parser.add_argument("--output", required=True, help="Output file containing sequences")
    parser.add_argument("--n", required=True, type=int, help="Number of sequences to extract")
    
    args = parser.parse_args()
    
    extract_sequences(args.input, args.output, args.n)

if __name__ == "__main__":
    main()
