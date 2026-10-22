"""Read pattern dictionaries used by compression commands."""

import json

TOKEN_OPEN = "<"
TOKEN_CLOSE = ">"


def _validate_patterns(patterns, pattern_file):
    if not isinstance(patterns, dict):
        raise ValueError(f"{pattern_file}: pattern dictionary must be a JSON object")

    for token, data in patterns.items():
        if not token or TOKEN_OPEN in token or TOKEN_CLOSE in token:
            raise ValueError(
                f"{pattern_file}: invalid token {token!r} (must be non-empty and "
                f"must not contain '{TOKEN_OPEN}' or '{TOKEN_CLOSE}')"
            )
        if not isinstance(data, dict):
            raise ValueError(f"{pattern_file}: entry for token {token!r} must be an object")
        sequence = data.get("sequence")
        if not isinstance(sequence, str) or not sequence:
            raise ValueError(f"{pattern_file}: token {token!r} has no 'sequence'")
        if TOKEN_OPEN in sequence or TOKEN_CLOSE in sequence:
            raise ValueError(
                f"{pattern_file}: sequence for token {token!r} contains a reserved character"
            )
    return patterns


def load_pattern_dictionary(pattern_file):
    with open(pattern_file, "r", encoding="utf-8") as pattern_handle:
        patterns = json.load(pattern_handle)
    return _validate_patterns(patterns, pattern_file)


def _priority_key(item):
    _, data = item
    priority = data.get("priority")
    if isinstance(priority, (int, float)):
        return (0, priority, 0)
    return (1, 0, -len(data["sequence"]))


def load_sorted_patterns(pattern_file):
    patterns = load_pattern_dictionary(pattern_file)
    return sorted(patterns.items(), key=_priority_key)


def load_token_map(pattern_file):
    patterns = load_pattern_dictionary(pattern_file)
    return {
        token: data["sequence"]
        for token, data in patterns.items()
    }
