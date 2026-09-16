FROM python:3.10-slim

ENV DEBIAN_FRONTEND=noninteractive \
    DNA_COMPRESSION_LOG_DIR=/data/logs \
    PYTHONUNBUFFERED=1 \
    UV_NO_DEV=1

RUN apt-get update && apt-get install -y --no-install-recommends \
        git \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY --from=ghcr.io/astral-sh/uv:0.12.15 /uv /uvx /bin/

WORKDIR /app

COPY pyproject.toml uv.lock README.md ./
RUN uv sync --locked --all-extras --no-install-project

COPY dna_compression ./dna_compression
RUN uv sync --locked --all-extras

RUN .venv/bin/python -c "from huggingface_hub import snapshot_download; \
p = snapshot_download(repo_id='zhihan1996/DNABERT-2-117M'); \
print('Model Downloaded:', p)"

ENV PATH="/app/.venv/bin:$PATH"

VOLUME ["/data"]

ENTRYPOINT ["dna-compress", "detect"]
CMD ["--help"]