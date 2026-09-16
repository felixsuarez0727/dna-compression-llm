# Architecture

## Goal

The project is an installable Python package with a unified CLI. This
organization preserves the existing compression behavior, prompts, token
syntax, and pattern JSON schema.

## Module Ownership

```text
dna_compression/
  cli.py                  Command parsing and command dispatcher
  config.py               Provider defaults and log directory configuration
  logging_utils.py        Shared UTF-8 log writer
  pattern_store.py        Pattern dictionary readers
  compression.py          Token replacement and expansion
  fastq.py                FASTQ extraction and synthetic data generation
  integrity.py            Whitespace-tolerant comparison helper
  benchmark.py            gzip, bzip2, and LZMA benchmark
  overhead.py             Deterministic overhead analysis
  detection/
    providers.py          Hosted provider clients
    local_models.py       DNABERT-2 and HyenaDNA helpers
    pipeline.py           Detection orchestration and optimization
```

`cli.py` parses command-line arguments and delegates to the owning module.
Library modules do not depend on the CLI. `dna-compress` is the only supported
command-line entry point.

## Dependencies

The core package uses only the Python standard library. Optional capabilities
are isolated in uv extras:

| Extra | Purpose |
| --- | --- |
| `llm` | OpenAI-compatible, DeepSeek, and Gemini hosted providers |
| `local-models` | CPU PyTorch and Transformer dependencies for local models |

Heavy imports are lazy. A core installation can import the package and use
offline commands without installing client or model libraries.

The `local-models` extra intentionally retains packages that may be imported
by remote DNABERT-2 and HyenaDNA code through `trust_remote_code=True`.

## Tests and CI

The `test` dependency group installs pytest. Tests cover the offline CLI,
compression round trips, FASTQ utilities, overhead analysis, provider request
adapters, and a simulated detection run. No API key, model download, or Docker
daemon is required.

GitHub Actions runs `uv lock --check`, synchronizes the test group, executes
pytest, and builds the package on Python 3.10 and 3.11 for Linux and Windows.

## Artifacts

Pattern dictionaries are JSON mappings from a token to metadata. The final
format remains compatible with prior outputs:

```json
{
  "!": {
    "sequence": "ATCGATCG",
    "priority": 1,
    "count": 42,
    "potential_savings": 294
  }
}
```

Logs default to `data/logs`. `DNA_COMPRESSION_LOG_DIR` changes that location;
the Docker image sets it to `/data/logs` so logs persist in the mounted volume.

## Deliberately Preserved Behavior

- Pattern detection uses its existing batching, prompts, validation, phase
  deduplication, and token allocation logic.
- `compare` ignores blank lines and leading/trailing whitespace, as did the
  original checker. It is not a byte-for-byte validator.
- The detector estimates savings with `seq_len - 1`; overhead analysis uses
  `seq_len - 3` for its `<X>` token model. Resolving that discrepancy is a
  separate algorithmic change.
- The project is CPU-only by default. CUDA support, linting, type checks, PyPI
  publication, and automatic `.env` loading are intentionally out of scope.