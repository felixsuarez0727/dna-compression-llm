import argparse

def extract_sequences(fastq_file, output_file, n):
    with open(fastq_file, "r") as f, open(output_file, "w") as out:
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
    
    parser.add_argument("--input", help="Input FASTQ file")
    parser.add_argument("--output", help="Output file containing sequences")
    parser.add_argument("--n", type=int, help="Number of sequences to extract")
    
    args = parser.parse_args()
    
    extract_sequences(args.input, args.output, args.n)

if __name__ == "__main__":
    main()
