import json
from types import SimpleNamespace

import pytest

from dna_compression.detection import pipeline


def test_parse_model_response_normalizes_and_keeps_highest_count():
    response = "\n".join(
        [
            "ATCGATCG: 2",
            "atcgatcg appears 5 times",
            "NNNN: 3",
            "invalid: 99",
            "A" * 201,
        ]
    )

    assert pipeline.parse_model_response(response) == {
        "ATCGATCG": 5,
        "NNNN": 3,
    }


def test_deduplicate_phase_variants_keeps_the_highest_count_per_length(monkeypatch):
    monkeypatch.setattr(pipeline, "save_log", lambda line: None)

    deduplicated = pipeline.deduplicate_phase_variants(
        {
            "ATCGATCG": 2,
            "TCGATCGA": 5,
            "GCGCGCGC": 3,
        }
    )

    assert deduplicated == {
        "TCGATCGA": 5,
        "GCGCGCGC": 3,
    }


def test_run_detection_writes_a_valid_dictionary_with_a_fake_provider(tmp_path, monkeypatch):
    input_file = tmp_path / "sequences.seq.txt"
    output_file = tmp_path / "patterns.json"
    input_file.write_text(
        "ATCGATCGATCGATCG\nATCGATCGATCGATCG\nGGGGGGGGGGGGGGGG\n",
        encoding="utf-8",
    )

    monkeypatch.setattr(pipeline, "save_log", lambda line: None)
    monkeypatch.setattr(
        pipeline,
        "build_client",
        lambda provider, api_key, save_log, model_id=None: object(),
    )
    monkeypatch.setattr(
        pipeline,
        "analyze_batch_with_context",
        lambda client, provider, model_id, sequences, batch_num, top_kmers: (
            "ATCGATCGATCGATCG:2\nGGGGGGGGGGGGGGGG:1"
        ),
    )
    monkeypatch.setattr(
        pipeline,
        "analyze_batch_synthesis",
        lambda client, provider, model_id, candidates, top_kmers: "ATCGATCGATCGATCG:2",
    )

    pipeline.run_detection(
        SimpleNamespace(
            file=str(input_file),
            provider="gemini",
            key="test-key",
            output=str(output_file),
            batch_size=3,
            model=None,
            max_batches=1,
            threads=1,
            overhead=5,
            no_validate=False,
            no_expand=False,
        )
    )

    patterns = json.loads(output_file.read_text(encoding="utf-8"))
    assert patterns
    assert all(
        {"sequence", "priority", "count", "potential_savings"} <= set(pattern)
        for pattern in patterns.values()
    )
    assert any(pattern["sequence"] == "ATCGATCGATCGATCG" for pattern in patterns.values())

def test_optimize_patterns_drops_patterns_that_would_grow_the_file(monkeypatch):
    monkeypatch.setattr(pipeline, "save_log", lambda line: None)

    result = pipeline.optimize_patterns(
        {
            "a": {"sequence": "GG", "count": 500},
            "b": {"sequence": "ATCGATCGATCG", "count": 50},
        },
        dictionary_overhead=0,
    )

    sequences = [entry["sequence"] for entry in result.values()]
    assert sequences == ["ATCGATCGATCG"]
    assert result[next(iter(result))]["potential_savings"] == 50 * (12 - 3)


def test_validate_and_count_pattern_counts_non_overlapping_hits():
    assert pipeline.validate_and_count_pattern("ATAT", ["ATATATATAT"]) == 2


def test_hosted_provider_requires_an_api_key(tmp_path, monkeypatch):
    from dna_compression.config import resolve_api_key

    monkeypatch.delenv("GEMINI_API_KEY", raising=False)
    with pytest.raises(ValueError, match="GEMINI_API_KEY"):
        resolve_api_key("gemini", None)

    monkeypatch.setenv("GEMINI_API_KEY", "from-env")
    assert resolve_api_key("gemini", None) == "from-env"
    assert resolve_api_key("gemini", "explicit") == "explicit"
    assert resolve_api_key("dnabert2", None) is None


def test_custom_local_model_is_passed_without_mutating_defaults(tmp_path, monkeypatch):
    from dna_compression.config import PROVIDER_DEFAULTS

    captured = {}
    input_file = tmp_path / "sequences.seq.txt"
    input_file.write_text("ATCGATCGATCGATCG\n", encoding="utf-8")
    default_before = PROVIDER_DEFAULTS["dnabert2"]

    def fake_build_client(provider, api_key, save_log, model_id=None):
        captured["model_id"] = model_id
        return object()

    monkeypatch.setattr(pipeline, "save_log", lambda line: None)
    monkeypatch.setattr(pipeline, "build_client", fake_build_client)
    monkeypatch.setattr(pipeline, "analyze_batch_with_context", lambda *args: None)

    pipeline.run_detection(
        SimpleNamespace(
            file=str(input_file),
            provider="dnabert2",
            key=None,
            output=str(tmp_path / "patterns.json"),
            batch_size=10,
            model="my-org/custom-model",
            max_batches=0,
            threads=1,
            overhead=5,
            no_validate=False,
            no_expand=False,
        )
    )

    assert captured["model_id"] == "my-org/custom-model"
    assert PROVIDER_DEFAULTS["dnabert2"] == default_before
