from dna_compression.motif_tools import (
    build_dictionary,
    expand_iupac,
    format_latex,
    run_comparison,
)


def test_expand_iupac():
    assert sorted(expand_iupac("GCM")) == ["GCA", "GCC"]
    assert expand_iupac("GCX") == []
    assert expand_iupac("N" * 4) == []


def test_build_dictionary_drops_unseen_candidates():
    sequences = ["ACGTACGTACGTACGT" * 4] * 20
    dictionary = build_dictionary({"ACGTACGT": 1, "TTTTTTTT": 5}, sequences, overhead=5)
    found = {entry["sequence"] for entry in dictionary.values()}
    assert "TTTTTTTT" not in found
    assert any(sequence.startswith("ACGT") for sequence in found)


def test_run_comparison_with_existing_dictionary(tmp_path):
    input_file = tmp_path / "reads.seq.txt"
    input_file.write_text("ACGTACGTACGTACGT\n" * 30, encoding="utf-8")
    dictionary = build_dictionary(
        {"ACGTACGTACGTACGT": 30}, input_file.read_text().split(), overhead=5
    )
    import json

    patterns = tmp_path / "p.json"
    patterns.write_text(json.dumps(dictionary), encoding="utf-8")
    rows = run_comparison(
        input_file, tmp_path / "out", tools=[], overhead=5,
        extra_dictionaries=[("llm", patterns)],
    )
    assert rows[0]["lossless"] and rows[0]["ratio_pct"] > 0
    assert "\\begin{table}" in format_latex(rows, "cap")
