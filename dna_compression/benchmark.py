"""Benchmark standard compressors against a sequence file."""

import bz2
import gzip
import lzma
import os
import time


COMPRESSORS = [
    ("gzip", lambda data: gzip.compress(data, compresslevel=9)),
    ("bzip2", lambda data: bz2.compress(data, compresslevel=9)),
    ("lzma", lambda data: lzma.compress(data, preset=9)),
]


def human(n):
    for unit in ["B", "KB", "MB", "GB"]:
        if n < 1024:
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} TB"


def benchmark_file(path):
    if not os.path.exists(path):
        print(f"File not found: {path}")
        return False

    with open(path, "rb") as input_handle:
        data = input_handle.read()

    original = len(data)
    print(f"\nInput file : {path}")
    print(f"Original   : {human(original)} ({original:,} bytes)")
    print(f"\n{'Method':<10} {'Output':>10} {'Ratio':>10} {'Time':>8}")
    print("-" * 44)

    for name, compressor in COMPRESSORS:
        start_time = time.time()
        compressed = compressor(data)
        elapsed = time.time() - start_time
        size = len(compressed)
        ratio = 100 * (1 - size / original)
        print(f"{name:<10} {human(size):>10} {ratio:>9.2f}%  {elapsed:>6.2f}s")

    print()
    return True