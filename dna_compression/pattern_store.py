"""Read pattern dictionaries used by compression commands."""

import json


def load_pattern_dictionary(pattern_file):
    with open(pattern_file, "r", encoding="utf-8") as pattern_handle:
        return json.load(pattern_handle)


def load_sorted_patterns(pattern_file):
    patterns = load_pattern_dictionary(pattern_file)
    return sorted(patterns.items(), key=lambda item: item[1]["priority"])


def load_token_map(pattern_file):
    patterns = load_pattern_dictionary(pattern_file)
    return {
        token: data["sequence"]
        for token, data in patterns.items()
    }