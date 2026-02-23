import argparse

def sanitized_files_are_equal(file1, file2):
    with open(file1, 'r', encoding='utf-8') as f1, open(file2, 'r', encoding='utf-8') as f2:
        lines1 = [line.strip() for line in f1 if line.strip()]
        lines2 = [line.strip() for line in f2 if line.strip()]
        return lines1 == lines2

def main():
    parser = argparse.ArgumentParser(description="Compare two files ignoring spaces and empty lines.")
    parser.add_argument("file1", help="Path to the first file")
    parser.add_argument("file2", help="Path to the second file")
    args = parser.parse_args()

    if sanitized_files_are_equal(args.file1, args.file2):
        print("The files are equal (ignoring spaces and empty lines)")
    else:
        print("The files are different")

if __name__ == "__main__":
    main()
