import json
from types import SimpleNamespace

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
        lambda provider, api_key, save_log: object(),
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