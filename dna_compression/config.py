"""Shared configuration for DNA compression commands."""

import os
from pathlib import Path


MICRO_TOKEN_POOL = list("bdefhijklmoprsvwxyz")

PROVIDER_DEFAULTS = {
    "deepseek": "deepseek-chat",
    "chatgpt": "gpt-4o-mini",
    "gemini": "gemini-2.5-flash-lite",
    "dnabert2": "zhihan1996/DNABERT-2-117M",
    "hyenadna": "LongSafari/hyenadna-large-1m-seqlen-hf",
}


def get_logs_dir():
    return Path(os.environ.get("DNA_COMPRESSION_LOG_DIR", "data/logs"))