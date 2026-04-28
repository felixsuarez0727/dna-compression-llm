# DNA Compression With LLM
## Table of Contents

- [How It Works](#how-it-works)
- [How the Pattern Detector Works](#how-the-pattern-detector-works)
  - [Key Features](#key-features)
  - [How the Pipeline Works](#how-the-pipeline-works)
    - [1. Statistical Pre-Analysis](#1-statistical-pre-analysis)
    - [2. Batch Processing & AI Analysis](#2-batch-processing--ai-analysis)
    - [3. Expansion & Deduplication](#3-expansion--deduplication)
    - [4. Mathematical Optimization](#4-mathematical-optimization)
    - [5. Final Dictionary Generation](#5-final-dictionary-generation)
  - [Technical Specifications](#technical-specifications)
- [Installation](#installation)
- [Files](#files)
- [Step 1 - Download FASTQ File](#step-1---download-fastq-file)
- [Step 2 - Setting Python Environment](#step-2---setting-python-environment)
- [Step 3 - Dependencies](#step-3---dependencies)
- [Step 4 - Prepare the Sequence File](#step-4---prepare-the-sequence-file)
- [Step 5 - Analyze Overhead](#step-5---analyze-overhead)
  - [How It Works](#how-it-works-1)
    - [1. Statistical "Ground Truth" Analysis](#1-statistical-ground-truth-analysis)
    - [2. Empirical Measurement of JSON "Cost"](#2-empirical-measurement-of-json-cost)
    - [3. Mathematical Simulation (The Sweep)](#3-mathematical-simulation-the-sweep)
    - [4. Optimization Discovery](#4-optimization-discovery)
  - [Key Takeaway for the User](#key-takeaway-for-the-user)
- [Step 6 - Detect Patterns](#step-6---detect-patterns)
- [Step 7 - Compress](#step-7---compress)
- [Step 8 - Decompress](#step-8---decompress)
- [Step 9 - Verify Integrity](#step-9---verify-integrity)
- [About compression_gzip_bz2_lzma_benchmark.py](#about-compression_gzip_bz2_lzma_benchmarkpy)
- [Conclusions](#conclusions)
- [Notes](#notes)
- [Citation](#citation)

A pipeline for compressing DNA sequence files by detecting repetitive patterns using Large Language Models. The LLM analyzes batches of sequences, identifies tandem repeat structures, and produces a pattern dictionary that a compressor uses to substitute long repetitive motifs with short tokens. A paired decompressor fully restores the original sequences.

---

## How It Works

1. Sequences are extracted from a FASTQ file into a plain text file (one sequence per line).
2. The pattern detector sends batches of sequences to an LLM, which identifies tandem repeat units and their extended forms.
3. The detected patterns are validated against the real sequences, deduplicated, and optimized into a JSON dictionary.
4. The compressor reads that dictionary and substitutes each matching pattern with a single-character token.
5. The decompressor reads the same dictionary and reverses every substitution to restore the original sequences.
6. The integrity checker confirms the restored file is byte-for-byte identical to the original.

---

## How the Pattern Detector Works?
The **FASTQ Pattern Detector** is a hybrid bioinformatics pipeline designed to identify, validate, and optimize DNA motifs for sequence compression. It leverages the pattern-recognition strengths of Large Language Models (LLMs) and anchors them with deterministic statistical validation to ensure 100% data integrity.

---

###  Key Features

* **Hybrid Intelligence:** Combines traditional k-mer frequency analysis with advanced LLM reasoning (GPT-4, DeepSeek, or Gemini).
* **Multi-Model Integration:** Out-of-the-box support for OpenAI, DeepSeek, and Google Gemini APIs.
* **Automated Validation:** Every pattern suggested by the AI is cross-verified against the raw FASTQ data to eliminate "hallucinations."
* **Tandem Repeat Expansion:** Programmatically identifies varying lengths of the same repeating motif (e.g., 2x, 4x, or 8x repeats).
* **Token Optimization:** Assigns single-character tokens to patterns based on their potential to reduce file size.

---

###  How the Pipeline Works

#### 1. Statistical Pre-Analysis
The script performs a preliminary scan of the sequences to find the most frequent **k-mers** (sub-sequences) of lengths 4, 6, 8, 10, and 12. This "Ground Truth" data is passed to the LLM to provide context for its analysis.

#### 2. Batch Processing & AI Analysis
Sequences are processed in batches (default: 30) and analyzed by the LLM. The AI identifies:
* **Base Units:** The fundamental repeating motif.
* **Extended Repeats:** The full tandem sequence found in the data.
* **NNN Prefixes:** Specific variants commonly found in synthetic DNA.

#### 3. Expansion & Deduplication
The pipeline performs a **Phase-Variant Deduplication**. It recognizes that shifted patterns like `ATCGATCG` and `TCGATCGA` represent the same biological repeat and merges them into the most efficient form. It also programmatically searches for higher-order repetitions of found motifs.

#### 4. Mathematical Optimization
Patterns are ranked by their **Potential Savings**. The script calculates the efficiency of each pattern using the following formula:

$$Savings = Count 	imes (SequenceLength - 1)$$

Only patterns that provide a net gain after accounting for "Dictionary Overhead" (the cost of storing the pattern in the JSON file) are retained.

#### 5. Final Dictionary Generation
The result is a structured `patterns.json` file. Each entry maps a unique single-character token to a high-frequency DNA sequence, ready to be used by a compression engine.

---

### Technical Specifications

| Component | Detail |
| :--- | :--- |
| **Language** | Python 3.x |
| **Concurrency** | `ThreadPoolExecutor` (for OpenAI/DeepSeek) |
| **Input Format** | `.fastq` or raw sequence text |
| **Output Format** | Structured `.json` |
| **Default Model** | `gpt-4o-mini` / `deepseek-chat` / `gemini-2.0-flash` |


## Installation

Python 3.9 or later is required.

```
pip install -r requirements.txt
```

The requirements file installs the two API client libraries:

```
google-genai
openai
```

DeepSeek uses the OpenAI-compatible client pointed at `https://api.deepseek.com`. No additional library is needed for it.

---

## Files

| File | Purpose |
|---|---|
| `extract_n.py` | Extracts N sequences from a FASTQ file into a plain text file |
| `pattern_detector.py` | Sends sequences to an LLM and produces a pattern JSON dictionary |
| `fastq_compressor.py` | Compresses a sequence file using a pattern dictionary |
| `fastq_decompressor.py` | Restores a compressed file to the original sequences |
| `check_files.py` | Verifies that two files are identical (ignoring whitespace) |

---
## Step 1 - Download FASTQ File
1. Download the testing file from: https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR15993673

2. Then place the FASTQ file alone in this folder: `/fastq_files`

## Step 2 - Setting Python Environment
```bash
python -m venv venv

CMD: .\venv\Scripts\activate.bat 
PowerShell: .\venv\Scripts\Activate.ps1
Bash: source venv/bin/activate  
```

---

## Step 3 - Dependencies

Make sure to install:

```bash
pip install -r requirements.txt
```


## Step 4 - Prepare the Sequence File

Extract a number of sequences from a FASTQ file into a plain text file. Each output line contains one raw sequence.

```
python extract_n.py --input ./fast_files/ERR15993673.fastq --output ./data/ERR15993673_5000.seq.txt --n 5000
```

| Argument | Description |
|---|---|
| `--input` | Input FASTQ file |
| `--output` | Output plain text file |
| `--n` | Number of sequences to extract |

---
## Step 5 - Analyze Overhead

This script is a diagnostic utility designed to fine-tune the compression efficiency of the DNA pipeline.

Its primary purpose is to determine the optimal `--overhead` value — a threshold used to decide if a DNA pattern is frequent enough to justify the "cost" of storing it in the JSON dictionary.

---

## How It Works

The script performs a simulation that mimics the compression logic without requiring an expensive LLM call.

### 1. Statistical "Ground Truth" Analysis

- **K-mer Scan**  
  It scans the `.seq.txt` file for frequent sub-sequences (lengths 4–12) to create a list of candidate patterns.

- **Deterministic Expansion**  
  It identifies *tandem repeats* (e.g., seeing `ATCG` and finding `ATCGATCGATCG`) and removes redundant shifted variants to ensure the pattern list is clean.

---

### 2. Empirical Measurement of JSON "Cost"

- The script calculates the real byte size of a JSON dictionary entry.
- It measures the **Fixed Metadata Cost** — the bytes taken up by JSON syntax (brackets, quotes, labels like `"sequence"`, `"count"`).
- This cost is identified as approximately **45–50 bytes per entry**.

---

### 3. Mathematical Simulation (The Sweep)

The script runs a simulation across a range of overhead values (1 to 200).

For each value, it applies a corrected savings formula:

$$
Net\ Savings = (Count \times (SequenceLength - 3)) - RealEntrySize
$$

- **Correction**  
  Unlike the main pipeline (which assumes a 1-character cost), this script accounts for the 3-byte `<X>` token format used in the compressed file.

---

### 4. Optimization Discovery

- It compares:
  - **Gross Savings** (bytes removed from the sequence)
  - **Dictionary Cost** (bytes added by the JSON file)

- It identifies the **"Sweet Spot"** where:
Net Savings = Gross Savings - Cost is maximized.

---

## Key Takeaway for the User

The script reveals that the default `--overhead 5` in the main pipeline is likely too aggressive  
(i.e., it keeps patterns that actually increase the total file size).

It typically recommends a higher overhead (e.g., **~45–50**) to ensure that every pattern stored in the dictionary provides a definitive net gain in compression.



## Step 6 - Detect Patterns

Run the pattern detector against the sequence file. The script sends batches to the chosen LLM, aggregates detected patterns across all batches, runs a synthesis pass, optionally expands tandem repeat variants, validates every candidate against the real sequences, deduplicates phase variants, and finally produces an optimized token dictionary saved as JSON.

### Using DeepSeek

```
python pattern_detector.py --threads 8 -f .\data\ERR15993673_5000.seq.txt -o .\data\ERR15993673_5000_deepseek_patterns.json -b 80 -p deepseek -m deepseek-chat -k [KEY]
```

### Using ChatGPT

```
python pattern_detector.py --threads 8 -f .\data\ERR15993673_5000.seq.txt -o .\data\ERR15993673_5000_chatgpt_patterns.json -b 80 -p chatgpt -m gpt-4o-mini -k [KEY]
```

### Using Gemini

```
python pattern_detector.py -f .\data\ERR15993673_5000.seq.txt -o .\data\ERR15993673_5000_gemini_patterns.json -b 80 -p gemini -m gemini-2.0-flash -k [KEY]
```

### All Arguments

| Argument | Short | Required | Default | Description |
|---|---|---|---|---|
| `--file` | `-f` | Yes | | Input sequence text file |
| `--key` | `-k` | Yes | | API key for the selected provider |
| `--provider` | `-p` | Yes | | LLM provider: `deepseek`, `chatgpt`, or `gemini` |
| `--output` | `-o` | No | `patterns.json` | Output JSON pattern file |
| `--batch_size` | `-b` | No | `30` | Number of sequences per API call |
| `--model` | `-m` | No | Provider default | Model ID. Defaults: `deepseek-chat`, `gpt-4o-mini`, `gemini-2.0-flash` |
| `--max_batches` | | No | `0` (all) | Stop after this many batches |
| `--threads` | | No | `4` | Parallel threads. Applies to `chatgpt` and `deepseek` only |
| `--overhead` | | No | `5` | Dictionary entry cost in bytes used during final optimization |
| `--no_validate` | | No | | Skip validation of patterns against real sequences |
| `--no_expand` | | No | | Skip tandem repeat variant expansion |

Threading is not available for Gemini. The `--threads` argument is silently ignored when `--provider gemini` is used.

---

## Step 7 - Compress

### Using DeepSeek
```
python fastq_compressor.py -i .\data\ERR15993673_5000.seq.txt -o .\data\ERR15993673_5000_deepseek.seq.compress -p .\data\ERR15993673_5000_deepseek_patterns.json
```

### Using ChatGPT
```
python fastq_compressor.py -i .\data\ERR15993673_5000.seq.txt -o .\data\ERR15993673_5000_chatgpt.seq.compress -p .\data\ERR15993673_5000_chatgpt_patterns.json
```

### Using Gemini
```
python fastq_compressor.py -i .\data\ERR15993673_5000.seq.txt -o .\data\ERR15993673_5000_gemini.seq.compress -p .\data\ERR15993673_5000_gemini_patterns.json
```

### All Arguments
| Argument | Short | Description |
|---|---|---|
| `--input` | `-i` | Input sequence file |
| `--patterns` | `-p` | Pattern JSON file produced by `pattern_detector.py` |
| `--output` | `-o` | Output compressed file |

The script reports the number of sequences processed and the compression ratio on completion.

---

## Step 8 - Decompress

### Using DeepSeek
```
python fastq_decompressor.py -i .\data\ERR15993673_5000_deepseek.seq.compress -o .\data\ERR15993673_5000_deepseek.seq.restored -p .\data\ERR15993673_5000_deepseek_patterns.json
```

### Using ChatGPT
```
python fastq_decompressor.py -i .\data\ERR15993673_5000_chatgpt.seq.compress -o .\data\ERR15993673_5000_chatgpt.seq.restored -p .\data\ERR15993673_5000_chatgpt_patterns.json
```

### Using Gemini
```
python fastq_decompressor.py -i .\data\ERR15993673_5000_gemini.seq.compress -o .\data\ERR15993673_5000_gemini.seq.restored -p .\data\ERR15993673_5000_gemini_patterns.json
```

| Argument | Short | Description |
|---|---|---|
| `--input` | `-i` | Input compressed file |
| `--patterns` | `-p` | Pattern JSON file used during compression |
| `--output` | `-o` | Output restored sequence file |

The same pattern JSON file used during compression must be supplied here. Without it the tokens cannot be resolved and decompression will fail.

---

## Step 9 - Verify Integrity

Confirm the restored file is identical to the original.

### Using DeepSeek
```
python check_files.py .\data\ERR15993673_5000.seq.txt .\data\ERR15993673_5000_deepseek.seq.restored
```
### Using ChatGPT
```
python check_files.py .\data\ERR15993673_5000.seq.txt .\data\ERR15993673_5000_chatgpt.seq.restored
```

### Using Gemini
```
python check_files.py .\data\ERR15993673_5000.seq.txt .\data\ERR15993673_5000_gemini.seq.restored
```

The script compares both files line by line, ignoring leading and trailing whitespace and empty lines. It prints whether the files are equal or different.
---

## About compression_gzip_bz2_lzma_benchmark.py
This code was built for testing how much the common compression algorithms perform over the example used in this project.
This is the list of algorithms:

- Gzip: Uses the DEFLATE algorithm, balancing speed and compression.

- Bzip2: Uses the Burrows-Wheeler transform; typically slower than gzip but achieves higher compression ratios.

- LZMA: The algorithm behind 7-Zip; it usually provides the highest compression ratio at the cost of significant memory and time.


---
## Conclusions
Using the exact same DNA file and the same 5,000 sequences, the only thing that changed in your experiment was the model used to detect patterns. Everything else in the pipeline stayed identical. Even so, the compressed file sizes were quite different. Gemini produced the smallest file, while ChatGPT and DeepSeek resulted in noticeably larger ones, with Gemini achieving about 30% better compression overall.

In simple terms, this means Gemini was better at spotting longer and more useful repeating DNA patterns, which allowed the compressor to replace more data with fewer tokens. ChatGPT performed reasonably well but was more conservative in the patterns it selected, and DeepSeek likely focused on shorter or less efficient repeats. Since the input data and processing steps were the same, these differences come directly from how each model understands and abstracts repetitive structure in the DNA sequences.

Please, consult the `data` folder to find the files that where generated during running the files: 
- `extract_n.py`
- `pattern_detector.py`
- `fastq_compressor.py`
- `fastq_decompressor.py`

Also consider testing with other FASTQ Files in order to explore behaviour in further DNA Sequences.
Just, for your notice, `generate_fastq.py` was a file intended to generate randomly a Fastq File, this may be executed like this: `python generate_fastq.py 10000 xl_data.fastq`

## Notes

The quality of compression depends directly on the patterns the LLM identifies. Results are not fully deterministic across separate runs because LLM responses can vary between sessions even at temperature zero, and because hosted model versions change silently over time. To ensure reproducibility, keep the JSON pattern file generated from a successful run. The compressor and decompressor are fully deterministic once the pattern file is fixed.

The `--threads` option accelerates the batch analysis phase by sending multiple API requests in parallel. Each thread creates its own isolated API client. For providers with strict rate limits, reduce thread count if batches return `no result` errors.

---

## Citation

Test data used during development:
```
ERR15993673 - NCBI Trace Archive
https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR15993673&display=analysis
```
This file is a raw metagenomic sequencing data from a human vaginal sample.