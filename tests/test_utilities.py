from dna_compression.benchmark import benchmark_file, human
from dna_compression.integrity import sanitized_files_are_equal
from dna_compression.overhead import (
    compute_top_kmers,
    find_tandem_unit,
    measure_json_entry_bytes,
    simulate_overhead,
    validate_and_count,
)


def test_sanitized_files_are_equal_ignores_blank_lines_and_outer_whitespace(tmp_path):
    first_file = tmp_path / "first.seq.txt"
    second_file = tmp_path / "second.seq.txt"
    first_file.write_text("  ATCG  \n\nTGCA\n", encoding="utf-8")
    second_file.write_text("ATCG\nTGCA\n", encoding="utf-8")

    assert sanitized_files_are_equal(first_file, second_file)

    second_file.write_text("ATCG\nAAAA\n", encoding="utf-8")
    assert not sanitized_files_are_equal(first_file, second_file)


def test_benchmark_reports_available_compressors_and_missing_input(tmp_path, capsys):
    input_file = tmp_path / "sequences.seq.txt"
    input_file.write_bytes(b"ATCG" * 100)

    assert benchmark_file(input_file)
    assert "gzip" in capsys.readouterr().out
    assert not benchmark_file(tmp_path / "missing.seq.txt")


def test_human_formats_binary_units():
    assert human(1023) == "1023.0 B"
    assert human(1024) == "1.0 KB"
    assert human(1024 * 1024) == "1.0 MB"


def test_overhead_kmers_and_overlapping_pattern_count():
    top_kmers = compute_top_kmers(["ATAT", "ATAT"], k_sizes=(2,), top_n=2)

    assert top_kmers[2] == [("AT", 4), ("TA", 2)]
    assert validate_and_count("ATA", ["ATATA"]) == 2
    assert find_tandem_unit("ATCGATCG") == ("ATCG", 2)


def test_simulate_overhead_uses_token_cost_and_json_cost():
    result = simulate_overhead({"ATCGATCG": 3}, overhead=1)

    assert result["retained"] == 1
    assert result["gross_savings"] == 15
    assert result["dict_cost"] == measure_json_entry_bytes(8)
    assert result["net_savings"] == result["gross_savings"] - result["dict_cost"]