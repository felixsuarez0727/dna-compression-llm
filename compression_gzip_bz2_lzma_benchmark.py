"""
benchmark_compressors.py
Compresses a sequence file with gzip, bzip2 and lzma and reports ratios.

Usage:
    python benchmark_compressors.py <input_file>

Example:
    python benchmark_compressors.py data/ERR15993673_5000.seq.txt
"""

import gzip
import bz2
import lzma
import os
import sys
import time

COMPRESSORS = [
    ("gzip",  lambda d: gzip.compress(d, compresslevel=9)),
    ("bzip2", lambda d: bz2.compress(d, compresslevel=9)),
    ("lzma",  lambda d: lzma.compress(d, preset=9)),
]

def human(n):
    for unit in ["B", "KB", "MB", "GB"]:
        if n < 1024:
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} TB"

def main():
    if len(sys.argv) < 2:
        print("Usage: python benchmark_compressors.py <input_file>")
        sys.exit(1)

    path = sys.argv[1]

    if not os.path.exists(path):
        print(f"File not found: {path}")
        sys.exit(1)

    with open(path, "rb") as f:
        data = f.read()

    original = len(data)
    print(f"\nInput file : {path}")
    print(f"Original   : {human(original)} ({original:,} bytes)")
    print(f"\n{'Method':<10} {'Output':>10} {'Ratio':>10} {'Time':>8}")
    print("-" * 44)

    for name, fn in COMPRESSORS:
        t0 = time.time()
        compressed = fn(data)
        elapsed = time.time() - t0
        size = len(compressed)
        ratio = 100 * (1 - size / original)
        print(f"{name:<10} {human(size):>10} {ratio:>9.2f}%  {elapsed:>6.2f}s")

    print()

if __name__ == "__main__":
    main()