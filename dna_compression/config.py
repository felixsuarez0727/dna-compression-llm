"""Shared configuration for DNA compression commands."""

import os
from pathlib import Path


MICRO_TOKEN_POOL = list("bdefhijklmoprsvwxyz")

TOKEN_BYTES = 3

LOCAL_PROVIDERS = ("dnabert2", "hyenadna")

API_KEY_ENV_VARS = {
    "chatgpt": "OPENAI_API_KEY",
    "deepseek": "DEEPSEEK_API_KEY",
    "gemini": "GEMINI_API_KEY",
}

PROVIDER_DEFAULTS = {
    "deepseek": "deepseek-chat",
    "chatgpt": "gpt-4o-mini",
    "gemini": "gemini-2.5-flash-lite",
    "dnabert2": "zhihan1996/DNABERT-2-117M",
    "hyenadna": "LongSafari/hyenadna-large-1m-seqlen-hf",
}


def get_logs_dir():
    return Path(os.environ.get("DNA_COMPRESSION_LOG_DIR", "data/logs"))

def resolve_api_key(provider, explicit_key=None):
    if provider in LOCAL_PROVIDERS:
        return None
    if explicit_key:
        return explicit_key
    env_var = API_KEY_ENV_VARS.get(provider)
    key = os.environ.get(env_var) if env_var else None
    if not key:
        raise ValueError(
            f"Provider '{provider}' needs an API key: set {env_var} or pass --key"
        )
    return key
