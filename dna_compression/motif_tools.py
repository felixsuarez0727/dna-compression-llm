"""Compare the pattern dictionaries of external motif tools with LLM-guided ones.

Each tool proposes candidate motifs. Every tool then goes through the same
validation, phase deduplication, token assignment and overhead filter as the
LLM pipeline, and through the same compressor, so results are comparable.
"""

import itertools
import json
import os
import random
import re
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

from .compression import compress_seq, decompress_seq
from .detection.pipeline import (
    build_compressor_patterns,
    deduplicate_phase_variants,
    expand_to_tandem_variants,
    find_tandem_unit,
    optimize_patterns,
    validate_and_count_pattern,
)
from .integrity import sanitized_files_are_equal
from .logging_utils import save_log

TOOLS = ("jellyfish", "trf", "mreps", "streme", "homer")
BIN_DIR_ENV = "MOTIF_TOOLS_BIN"

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "N": "ACGT",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG",
}
MAX_IUPAC_EXPANSION = 16
DNA_ONLY = re.compile(r"^[ACGT]+$")


def find_binary(name, bin_dir=None):
    executable = "homer2" if name == "homer" else name
    directory = bin_dir or os.environ.get(BIN_DIR_ENV)
    if directory:
        candidate = Path(directory) / executable
        if candidate.exists():
            return str(candidate)
    return shutil.which(executable)


def write_fasta(sequences, path):
    with open(path, "w", encoding="utf-8") as handle:
        for index, sequence in enumerate(sequences, start=1):
            handle.write(f">r{index}\n{sequence}\n")


def expand_iupac(consensus):
    choices = [IUPAC.get(char) for char in consensus.upper()]
    if any(choice is None for choice in choices):
        return []
    total = 1
    for choice in choices:
        total *= len(choice)
    if total > MAX_IUPAC_EXPANSION:
        return []
    return ["".join(combo) for combo in itertools.product(*choices)]


def _run(command, cwd=None, stdout=None):
    subprocess.run(
        command,
        cwd=cwd,
        stdout=stdout if stdout is not None else subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True,
    )


def candidates_jellyfish(binary, sequences, workdir, top=200, k_sizes=(4, 6, 8, 10, 12)):
    fasta = workdir / "reads.fa"
    write_fasta(sequences, fasta)
    candidates = {}
    for k in k_sizes:
        db = workdir / f"k{k}.jf"
        _run([binary, "count", "-m", str(k), "-s", "50M", "-t", "4", "-o", str(db), str(fasta)])
        dump = subprocess.run(
            [binary, "dump", "-c", str(db)], capture_output=True, text=True, check=True
        )
        counts = []
        for line in dump.stdout.splitlines():
            kmer, count = line.split()
            if DNA_ONLY.match(kmer):
                counts.append((int(count), kmer))
        counts.sort(reverse=True)
        for count, kmer in counts[:top]:
            candidates[kmer] = count
    return candidates


def candidates_trf(binary, sequences, workdir):
    fasta = workdir / "reads.fa"
    write_fasta(sequences, fasta)
    with open(workdir / "trf.out", "w") as handle:
        subprocess.run(
            [binary, str(fasta), "2", "7", "7", "80", "10", "30", "500", "-h", "-ngs"],
            cwd=workdir,
            stdout=handle,
            stderr=subprocess.DEVNULL,
            check=False,
        )
    candidates = {}
    for line in (workdir / "trf.out").read_text().splitlines():
        fields = line.split()
        if line.startswith("@") or len(fields) < 15:
            continue
        for sequence in (fields[13], fields[14]):
            sequence = sequence.upper()
            if len(sequence) >= 4 and DNA_ONLY.match(sequence):
                candidates[sequence] = candidates.get(sequence, 0) + 1
    return candidates


def candidates_mreps(binary, sequences, workdir):
    fasta = workdir / "reads.fa"
    write_fasta(sequences, fasta)
    with open(workdir / "mreps.out", "w") as handle:
        subprocess.run(
            [binary, "-fasta", "-res", "3", "-minsize", "8", str(fasta)],
            stdout=handle,
            stderr=subprocess.DEVNULL,
            check=False,
        )
    row = re.compile(r"^\s*\d+\s*->\s*\d+\s*:\s*\d+\s*<\d+>\s*\[[\d.]+\]\s*[\d.]+\s*(.+)$")
    candidates = {}
    for line in (workdir / "mreps.out").read_text().splitlines():
        match = row.match(line)
        if not match:
            continue
        sequence = match.group(1).replace(" ", "").upper()
        if len(sequence) >= 4 and DNA_ONLY.match(sequence):
            candidates[sequence] = candidates.get(sequence, 0) + 1
    return candidates


def candidates_streme(binary, sequences, workdir, max_motifs=20):
    fasta = workdir / "reads.fa"
    write_fasta(sequences, fasta)
    out = workdir / "streme_out"
    _run(
        [
            binary, "--dna", "--p", str(fasta), "--oc", str(out),
            "--minw", "4", "--maxw", "30", "--nmotifs", str(max_motifs),
        ]
    )
    candidates = {}
    for line in (out / "streme.txt").read_text().splitlines():
        if not line.startswith("MOTIF"):
            continue
        motif_id = line.split()[1]
        consensus = motif_id.split("-", 1)[-1]
        for sequence in expand_iupac(consensus):
            candidates[sequence] = candidates.get(sequence, 0) + 1
    return candidates


def candidates_homer(binary, sequences, workdir, seed=1):
    fasta = workdir / "reads.fa"
    write_fasta(sequences, fasta)
    rng = random.Random(seed)
    shuffled = []
    for sequence in sequences:
        letters = list(sequence)
        rng.shuffle(letters)
        shuffled.append("".join(letters))
    background = workdir / "background.fa"
    write_fasta(shuffled, background)
    motifs = workdir / "homer.motifs"
    _run(
        [
            binary, "denovo", "-i", str(fasta), "-b", str(background),
            "-len", "6,8,10,12", "-S", "10", "-o", str(motifs),
        ]
    )
    candidates = {}
    for line in motifs.read_text().splitlines():
        if not line.startswith(">"):
            continue
        consensus = line[1:].split("\t")[0]
        for sequence in expand_iupac(consensus):
            candidates[sequence] = candidates.get(sequence, 0) + 1
    return candidates


CANDIDATE_FUNCTIONS = {
    "jellyfish": candidates_jellyfish,
    "trf": candidates_trf,
    "mreps": candidates_mreps,
    "streme": candidates_streme,
    "homer": candidates_homer,
}


def build_dictionary(candidates, sequences, overhead, expand=True):
    aggregated = dict(candidates)
    if expand:
        checked = set()
        for pattern in list(aggregated):
            unit, _ = find_tandem_unit(pattern)
            if unit and unit not in checked:
                checked.add(unit)
                for variant, count in expand_to_tandem_variants(unit, sequences).items():
                    aggregated[variant] = max(aggregated.get(variant, 0), count)

    validated = {}
    for pattern in aggregated:
        real_count = validate_and_count_pattern(pattern, sequences)
        if real_count > 0:
            validated[pattern] = real_count

    validated = deduplicate_phase_variants(validated)
    compressor_ready = build_compressor_patterns(validated)
    return optimize_patterns(compressor_ready, overhead)


def base_coverage(sequences, dictionary):
    total = sum(len(sequence) for sequence in sequences)
    if total == 0:
        return 0.0
    covered = 0
    ordered = sorted(dictionary.values(), key=lambda entry: entry["priority"])
    for sequence in sequences:
        marker = [False] * len(sequence)
        for entry in ordered:
            pattern = entry["sequence"]
            start = sequence.find(pattern)
            while start != -1:
                if not any(marker[start:start + len(pattern)]):
                    marker[start:start + len(pattern)] = [True] * len(pattern)
                start = sequence.find(pattern, start + len(pattern))
        covered += sum(marker)
    return 100 * covered / total


def evaluate(name, dictionary, sequences, input_file, outdir, elapsed):
    patterns_file = outdir / f"{name}_patterns.json"
    compressed_file = outdir / f"{name}.seq.compress"
    restored_file = outdir / f"{name}.seq.restored"
    with open(patterns_file, "w", encoding="utf-8") as handle:
        json.dump(dictionary, handle, indent=4, ensure_ascii=False)

    log = f"motif_compare_{name}.txt"
    if dictionary:
        compress_seq(str(input_file), str(patterns_file), str(compressed_file), log)
        decompress_seq(str(compressed_file), str(patterns_file), str(restored_file), log)
        original_size = Path(input_file).stat().st_size
        compressed_size = compressed_file.stat().st_size
        ratio = 100 * (1 - compressed_size / original_size)
        lossless = sanitized_files_are_equal(str(input_file), str(restored_file))
    else:
        compressed_size = Path(input_file).stat().st_size
        ratio = 0.0
        lossless = True

    lengths = [len(entry["sequence"]) for entry in dictionary.values()]
    return {
        "tool": name,
        "patterns": len(dictionary),
        "mean_length": sum(lengths) / len(lengths) if lengths else 0.0,
        "coverage_pct": base_coverage(sequences, dictionary) if dictionary else 0.0,
        "ratio_pct": ratio,
        "compressed_kb": compressed_size / 1024,
        "dictionary_kb": patterns_file.stat().st_size / 1024,
        "seconds": elapsed,
        "lossless": lossless,
    }


def format_table(rows):
    header = (
        f"{'Tool':<12}{'Patterns':>9}{'Mean L':>9}{'Cover %':>9}"
        f"{'R %':>8}{'Compr. KB':>11}{'Dict KB':>9}{'Time s':>9}{'Lossless':>10}"
    )
    lines = [header, "-" * len(header)]
    for row in rows:
        lines.append(
            f"{row['tool']:<12}{row['patterns']:>9}{row['mean_length']:>9.1f}"
            f"{row['coverage_pct']:>9.2f}{row['ratio_pct']:>8.2f}"
            f"{row['compressed_kb']:>11.2f}{row['dictionary_kb']:>9.2f}"
            f"{row['seconds']:>9.1f}{'yes' if row['lossless'] else 'NO':>10}"
        )
    return "\n".join(lines)


def format_latex(rows, caption):
    lines = [
        "\\begin{table}[htbp]",
        "\\centering",
        f"\\caption{{{caption}}}",
        "\\label{tab:motif_tools}",
        "\\small",
        "\\setlength{\\tabcolsep}{3pt}",
        "\\resizebox{\\columnwidth}{!}{%",
        "\\begin{tabular}{lcccccc}",
        "\\toprule",
        "\\textbf{Tool} & \\textbf{Patterns} & \\textbf{Mean $L$ (bp)} & "
        "\\textbf{Coverage (\\%)} & \\textbf{Ratio $R$ (\\%)} & "
        "\\textbf{Time (s)} & \\textbf{Lossless} \\\\",
        "\\midrule",
    ]
    for row in rows:
        mark = "\\checkmark" if row["lossless"] else "no"
        lines.append(
            f"{row['tool']} & {row['patterns']} & {row['mean_length']:.1f} & "
            f"{row['coverage_pct']:.2f} & {row['ratio_pct']:.2f} & "
            f"{row['seconds']:.1f} & {mark} \\\\"
        )
    lines += ["\\bottomrule", "\\end{tabular}}", "\\end{table}"]
    return "\n".join(lines)


def run_comparison(input_file, outdir, tools, overhead, bin_dir=None, extra_dictionaries=()):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    with open(input_file, "r", encoding="utf-8") as handle:
        sequences = [line.strip() for line in handle if line.strip()]
    save_log(f"{len(sequences)} sequences loaded from {input_file}")

    rows = []
    for name in tools:
        binary = find_binary(name, bin_dir)
        if not binary:
            save_log(f"[{name}] binary not found; skipping (set {BIN_DIR_ENV} or --bin-dir)")
            continue
        save_log(f"[{name}] running...")
        start = time.time()
        with tempfile.TemporaryDirectory() as tmp:
            candidates = CANDIDATE_FUNCTIONS[name](binary, sequences, Path(tmp))
        save_log(f"[{name}] {len(candidates)} candidate motifs")
        dictionary = build_dictionary(candidates, sequences, overhead)
        elapsed = time.time() - start
        rows.append(evaluate(name, dictionary, sequences, input_file, outdir, elapsed))

    for label, path in extra_dictionaries:
        with open(path, "r", encoding="utf-8") as handle:
            dictionary = json.load(handle)
        rows.append(evaluate(label, dictionary, sequences, input_file, outdir, 0.0))

    table = format_table(rows)
    print("\n" + table + "\n")
    (outdir / "comparison.txt").write_text(table + "\n", encoding="utf-8")
    (outdir / "comparison.json").write_text(json.dumps(rows, indent=2), encoding="utf-8")
    caption = f"Motif tools vs.\\ LLM dictionaries ($N={len(sequences)}$, $\\delta={overhead}$)"
    (outdir / "comparison.tex").write_text(format_latex(rows, caption) + "\n", encoding="utf-8")
    return rows
