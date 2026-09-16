# DNA Compression With LLM

A lossless DNA sequence compression workflow. It extracts sequence lines from
FASTQ files, finds repetitive motifs with hosted or local models, stores the
selected patterns in JSON, replaces those patterns with tokens, and restores
the original sequence file from the same dictionary.

The project ships as an installable package with one command-line interface:
`dna-compress`.

## Requirements

- Python 3.10 or 3.11
- [uv](https://docs.astral.sh/uv/)

## Install

The core utilities use only the standard library:

```powershell
uv sync
uv run dna-compress --help
```

Install optional dependencies only for the detection providers you use:

```powershell
uv sync --extra llm
uv sync --extra local-models
uv sync --all-extras
```

| Extra | Includes |
| --- | --- |
| `llm` | OpenAI-compatible clients, DeepSeek, and Google Gemini |
| `local-models` | CPU PyTorch, Transformers, DNABERT-2, and HyenaDNA dependencies |

`pyproject.toml` and `uv.lock` are the project's dependency source of truth.
There is no maintained `requirements.txt` file.

## Commands

Run all commands as `uv run dna-compress <command>`:

| Command | Purpose |
| --- | --- |
| `detect` | Create a pattern dictionary with a hosted or local model |
| `compress` | Replace dictionary sequences with tokens |
| `decompress` | Restore a compressed sequence file |
| `analyze` | Sweep dictionary-overhead estimates without calling a model |
| `extract` | Extract sequence lines from FASTQ or FASTQ.GZ |
| `generate` | Generate synthetic FASTQ data |
| `compare` | Compare normalized sequence files |
| `benchmark` | Compare gzip, bzip2, and LZMA on an input file |

Use `uv run dna-compress <command> --help` for the full argument list.

## Typical Workflow

Download a FASTQ file from the
[NCBI Trace Archive](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR15993673)
and place it under `fastq_files/input/`.

Extract sequence lines:

```powershell
uv run dna-compress extract --input .\fastq_files\input\ERR15993673.fastq.gz --output .\data\input\ERR15993673_5000.seq.txt --n 5000
```

Estimate an appropriate detector overhead without an API call:

```powershell
uv run dna-compress analyze --file .\data\input\ERR15993673_5000.seq.txt
```

Detect patterns with a hosted provider. Install the `llm` extra first and pass
the provider key through `--key`:

```powershell
uv run dna-compress detect -f .\data\input\ERR15993673_5000.seq.txt -o .\data\outputs\patterns.json -p gemini -b 80 -k $env:GEMINI_API_KEY
```

Compress and restore using the generated JSON dictionary:

```powershell
uv run dna-compress compress -i .\data\input\ERR15993673_5000.seq.txt -p .\data\outputs\patterns.json -o .\data\outputs\sequences.compress
uv run dna-compress decompress -i .\data\outputs\sequences.compress -p .\data\outputs\patterns.json -o .\data\outputs\sequences.restored
```

`compare` intentionally keeps the original whitespace-tolerant behavior. For
a strict lossless check in PowerShell, compare hashes:

```powershell
(Get-FileHash .\data\input\ERR15993673_5000.seq.txt -Algorithm SHA256).Hash -eq (Get-FileHash .\data\outputs\sequences.restored -Algorithm SHA256).Hash
```

## Detection Providers

| Provider | Extra | API key | Default model |
| --- | --- | --- | --- |
| `chatgpt` | `llm` | Required | `gpt-4o-mini` |
| `deepseek` | `llm` | Required | `deepseek-chat` |
| `gemini` | `llm` | Required | `gemini-2.5-flash-lite` |
| `dnabert2` | `local-models` | Not required | `zhihan1996/DNABERT-2-117M` |
| `hyenadna` | `local-models` | Not required | `LongSafari/hyenadna-large-1m-seqlen-hf` |

Local models download their weights the first time they run. DNABERT-2 and
HyenaDNA retain their existing `trust_remote_code=True` behavior, so use the
locked optional environment supplied by this repository.

## Docker

The Docker image installs every locked optional dependency, pre-downloads
DNABERT-2, and starts at `dna-compress detect` by default.

```powershell
docker build -t dna-patterns .
docker run --rm -v "${PWD}\data:/data" dna-patterns --help
```

For a DNABERT-2 detection run:

```powershell
docker run --rm -v "${PWD}\data:/data" dna-patterns -f /data/input/ERR15993673_5000.seq.txt -o /data/outputs/dnabert2_patterns.json -p dnabert2 -b 80 --overhead 101
```

The container writes logs under the mounted `/data/logs` directory. Set
`DNA_COMPRESSION_LOG_DIR` for a different log directory outside Docker.

## Project Layout

| Location | Purpose |
| --- | --- |
| `dna_compression/` | Reusable library modules and command dispatcher |
| `dna_compression/detection/` | Provider clients, local-model helpers, and detection pipeline |
| `data/input/` | Extracted sequence input files |
| `data/outputs/` | Pattern JSON, compressed files, and restored files |
| `data/logs/` | Command logs |
| `docs/architecture.md` | Module ownership and deliberate compatibility boundaries |

See [docs/architecture.md](docs/architecture.md) for the package structure and
known scope boundaries.

## Notes

- The pattern JSON produced by `detect` must be retained for decompression.
- Hosted model output can vary over time, even when a temperature setting is
  fixed. Once a pattern JSON is saved, compression and decompression are
  deterministic.
- The detector's historical savings estimate and the overhead diagnostic use
  different token-cost assumptions. This organizational refactor preserves
  that existing behavior.
- The bundled test data is from `ERR15993673`, a raw metagenomic sequencing
  data set from a human vaginal sample in the NCBI Trace Archive.
