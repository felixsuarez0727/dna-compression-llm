import json

import pytest

from dna_compression.compression import (
    compress_seq,
    decompress_seq,
    expand_sequence,
    line_ending,
)
from dna_compression.pattern_store import load_sorted_patterns, load_token_map


def write_patterns(path):
    path.write_text(
        json.dumps(
            {
                "!": {"sequence": "ATCGATCG", "priority": 1},
                "#": {"sequence": "GG", "priority": 2},
            }
        ),
        encoding="utf-8",
    )


def test_round_trip_preserves_tokens_and_line_endings(tmp_path, monkeypatch):
    monkeypatch.setenv("DNA_COMPRESSION_LOG_DIR", str(tmp_path / "logs"))
    input_file = tmp_path / "input.seq.txt"
    patterns_file = tmp_path / "patterns.json"
    compressed_file = tmp_path / "output.compress"
    restored_file = tmp_path / "output.restored"
    input_file.write_bytes(b"ATCGATCG\r\nGGGG\r\nATCGATCG")
    write_patterns(patterns_file)

    compress_seq(input_file, patterns_file, compressed_file, "compress.log")
    decompress_seq(compressed_file, patterns_file, restored_file, "decompress.log")

    assert compressed_file.read_bytes() == b"<!>\r\n<#><#>\r\n<!>"
    assert restored_file.read_bytes() == input_file.read_bytes()
    assert (tmp_path / "logs" / "compress.log").exists()
    assert (tmp_path / "logs" / "decompress.log").exists()


def test_pattern_store_sorts_by_priority_and_builds_token_map(tmp_path):
    patterns_file = tmp_path / "patterns.json"
    write_patterns(patterns_file)

    assert [token for token, _ in load_sorted_patterns(patterns_file)] == ["!", "#"]
    assert load_token_map(patterns_file) == {"!": "ATCGATCG", "#": "GG"}


@pytest.mark.parametrize(
    ("line", "expected"),
    [
        ("ATCG\r\n", "\r\n"),
        ("ATCG\n", "\n"),
        ("ATCG\r", "\r"),
        ("ATCG", ""),
    ],
)
def test_line_ending_preserves_supported_terminators(line, expected):
    assert line_ending(line) == expected


@pytest.mark.parametrize(
    ("tokenized_sequence", "token_map", "message"),
    [
        ("ATCG<", {}, "Malformed token"),
        ("ATCG<?>", {}, "Unknown token"),
    ],
)
def test_expand_sequence_rejects_invalid_tokens(tokenized_sequence, token_map, message):
    with pytest.raises(ValueError, match=message):
        expand_sequence(tokenized_sequence, token_map)