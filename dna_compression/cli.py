"""Command-line entry point for DNA compression workflows."""

import argparse
import time
from pathlib import Path

from .benchmark import benchmark_file
from .compression import compress_seq, decompress_seq
from .detection.pipeline import add_arguments as add_detect_arguments
from .detection.pipeline import run_detection
from .fastq import extract_sequences, generate_fastq_file
from .integrity import sanitized_files_are_equal
from .motif_tools import TOOLS as MOTIF_TOOLS
from .motif_tools import run_comparison
from .overhead import DEFAULT_MAX_OVERHEAD, DEFAULT_TOP_KMERS, run_analysis


def _add_compress_arguments(parser):
    parser.add_argument("-i", "--input", required=True, help="Input sequence file")
    parser.add_argument("-p", "--patterns", required=True, help="Pattern JSON file")
    parser.add_argument("-o", "--output", required=True, help="Output compressed file")


def _add_decompress_arguments(parser):
    parser.add_argument(
        "-i", "--input", required=True, help="Input compressed sequence file"
    )
    parser.add_argument("-p", "--patterns", required=True, help="Pattern JSON file")
    parser.add_argument("-o", "--output", required=True, help="Output restored sequence file")


def _add_extract_arguments(parser):
    parser.add_argument("--input", required=True, help="Input FASTQ file")
    parser.add_argument("--output", required=True, help="Output file containing sequences")
    parser.add_argument("--n", required=True, type=int, help="Number of sequences to extract")


def _add_generate_arguments(parser):
    parser.add_argument("num_sequences", type=int, help="Number of sequences to generate")
    parser.add_argument("output_file", nargs="?", help="Output FASTQ file")


def _add_compare_arguments(parser):
    parser.add_argument("file1", help="Path to the first file")
    parser.add_argument("file2", help="Path to the second file")


def _add_benchmark_arguments(parser):
    parser.add_argument("file", help="Input sequence file")


def _add_analyze_arguments(parser):
    parser.add_argument("--file", "-f", required=True, help="Input .seq.txt file (one sequence per line)")
    parser.add_argument(
        "--max-overhead",
        type=int,
        default=DEFAULT_MAX_OVERHEAD,
        help="Upper bound of the overhead sweep range",
    )
    parser.add_argument(
        "--step",
        type=int,
        default=5,
        help="Step size for the overhead sweep (smaller = finer resolution)",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=DEFAULT_TOP_KMERS,
        help="Top-N k-mers per k-size to use as seed patterns",
    )


def _add_compare_motifs_arguments(parser):
    parser.add_argument("--file", "-f", required=True, help="Input .seq.txt file (one sequence per line)")
    parser.add_argument("--outdir", "-o", default="data/outputs/motif_comparison", help="Output directory")
    parser.add_argument(
        "--tools",
        nargs="+",
        default=list(MOTIF_TOOLS),
        choices=MOTIF_TOOLS,
        help="Motif tools to run (default: all)",
    )
    parser.add_argument("--overhead", type=int, default=101, help="Dictionary overhead used to filter patterns")
    parser.add_argument("--bin-dir", default=None, help="Directory containing the tool binaries")
    parser.add_argument(
        "--llm",
        nargs="*",
        default=[],
        metavar="LABEL=PATTERNS.json",
        help="Existing LLM pattern dictionaries to include in the table",
    )


def _run_compare_motifs(args):
    extra = []
    for item in args.llm:
        label, _, path = item.partition("=")
        if not path:
            raise SystemExit(f"--llm expects LABEL=path, got: {item}")
        extra.append((label, path))
    run_comparison(args.file, args.outdir, args.tools, args.overhead, args.bin_dir, extra)


def _run_compress(args):
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    log_file = f"{Path(args.output).stem}.compress_log_{timestamp}.txt"
    compress_seq(args.input, args.patterns, args.output, log_file)


def _run_decompress(args):
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    log_file = f"{Path(args.output).stem}.decompress_log_{timestamp}.txt"
    decompress_seq(args.input, args.patterns, args.output, log_file)


def _run_extract(args):
    extract_sequences(args.input, args.output, args.n)


def _run_generate(args):
    output_file = args.output_file or f"synthetic_{args.num_sequences}.fastq"
    generate_fastq_file(output_file, args.num_sequences, seq_length=100)


def _run_compare(args):
    if sanitized_files_are_equal(args.file1, args.file2):
        print("The files are equal (ignoring spaces and empty lines)")
    else:
        print("The files are different")


def _run_benchmark(args):
    if not benchmark_file(args.file):
        raise SystemExit(1)


def _run_analyze(args):
    run_analysis(args.file, args.max_overhead, args.step, args.top)


def _run_detect(args):
    run_detection(args)


COMMANDS = {
    "detect": ("Detect DNA patterns", add_detect_arguments, _run_detect),
    "compress": ("Compress a sequence file", _add_compress_arguments, _run_compress),
    "decompress": (
        "Restore a compressed sequence file",
        _add_decompress_arguments,
        _run_decompress,
    ),
    "extract": ("Extract sequences from a FASTQ file", _add_extract_arguments, _run_extract),
    "generate": ("Generate a synthetic FASTQ file", _add_generate_arguments, _run_generate),
    "compare": ("Compare sequence files", _add_compare_arguments, _run_compare),
    "benchmark": ("Benchmark standard compressors", _add_benchmark_arguments, _run_benchmark),
    "analyze": ("Analyze pattern dictionary overhead", _add_analyze_arguments, _run_analyze),
    "compare-motifs": (
        "Compare external motif tools with LLM dictionaries",
        _add_compare_motifs_arguments,
        _run_compare_motifs,
    ),
}


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="dna-compress",
        description="DNA sequence compression workflows",
    )
    subparsers = parser.add_subparsers(dest="command")

    for command, (description, add_arguments, handler) in COMMANDS.items():
        command_parser = subparsers.add_parser(command, help=description)
        add_arguments(command_parser)
        command_parser.set_defaults(handler=handler)

    args = parser.parse_args(argv)
    if not hasattr(args, "handler"):
        parser.print_help()
        return
    args.handler(args)


if __name__ == "__main__":
    main()