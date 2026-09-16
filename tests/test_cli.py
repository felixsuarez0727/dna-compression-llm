import json

from dna_compression import cli


def write_patterns(path):
    path.write_text(
        json.dumps({"!": {"sequence": "ATCG", "priority": 1}}),
        encoding="utf-8",
    )


def test_main_without_a_command_prints_available_commands(capsys):
    cli.main([])

    output = capsys.readouterr().out
    assert "DNA sequence compression workflows" in output
    assert "detect" in output
    assert "decompress" in output


def test_detect_command_parses_arguments_and_delegates(tmp_path, monkeypatch):
    captured = {}
    input_file = tmp_path / "sequences.seq.txt"

    monkeypatch.setattr(cli, "run_detection", lambda args: captured.update(vars(args)))

    cli.main(
        [
            "detect",
            "--file",
            str(input_file),
            "--provider",
            "gemini",
            "--key",
            "test-key",
            "--max_batches",
            "2",
        ]
    )

    assert captured["file"] == str(input_file)
    assert captured["provider"] == "gemini"
    assert captured["key"] == "test-key"
    assert captured["max_batches"] == 2


def test_compress_and_decompress_commands_restore_input(tmp_path, monkeypatch):
    monkeypatch.setenv("DNA_COMPRESSION_LOG_DIR", str(tmp_path / "logs"))
    input_file = tmp_path / "input.seq.txt"
    patterns_file = tmp_path / "patterns.json"
    compressed_file = tmp_path / "output.compress"
    restored_file = tmp_path / "output.restored"
    input_file.write_bytes(b"ATCGATCG\r\nGATTACA")
    write_patterns(patterns_file)

    cli.main(
        [
            "compress",
            "--input",
            str(input_file),
            "--patterns",
            str(patterns_file),
            "--output",
            str(compressed_file),
        ]
    )
    cli.main(
        [
            "decompress",
            "--input",
            str(compressed_file),
            "--patterns",
            str(patterns_file),
            "--output",
            str(restored_file),
        ]
    )

    assert restored_file.read_bytes() == input_file.read_bytes()


def test_compare_command_reports_normalized_equality(tmp_path, capsys):
    first_file = tmp_path / "first.seq.txt"
    second_file = tmp_path / "second.seq.txt"
    first_file.write_text("  ATCG  \n\nTGCA\n", encoding="utf-8")
    second_file.write_text("ATCG\nTGCA\n", encoding="utf-8")

    cli.main(["compare", str(first_file), str(second_file)])

    assert "The files are equal" in capsys.readouterr().out