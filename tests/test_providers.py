from types import SimpleNamespace

import pytest

from dna_compression.detection.providers import build_client, call_llm


def test_call_llm_uses_openai_compatible_request_shape():
    captured = {}

    def create(**kwargs):
        captured.update(kwargs)
        return SimpleNamespace(
            choices=[SimpleNamespace(message=SimpleNamespace(content="ATCG:2"))]
        )

    client = SimpleNamespace(
        chat=SimpleNamespace(completions=SimpleNamespace(create=create))
    )

    response = call_llm(client, "chatgpt", "test-model", "system", "user")

    assert response == "ATCG:2"
    assert captured["model"] == "test-model"
    assert captured["temperature"] == 0
    assert captured["messages"] == [
        {"role": "system", "content": "system"},
        {"role": "user", "content": "user"},
    ]


def test_call_llm_uses_gemini_request_shape():
    captured = {}

    def generate_content(**kwargs):
        captured.update(kwargs)
        return SimpleNamespace(text="GGGG:4")

    client = SimpleNamespace(models=SimpleNamespace(generate_content=generate_content))

    response = call_llm(client, "gemini", "gemini-test", "unused", "prompt")

    assert response == "GGGG:4"
    assert captured == {"model": "gemini-test", "contents": "prompt"}


@pytest.mark.parametrize("function", [build_client, call_llm])
def test_unknown_provider_is_rejected(function):
    if function is build_client:
        arguments = ("unknown", "key", lambda line: None)
    else:
        arguments = (object(), "unknown", "model", "system", "user")

    with pytest.raises(ValueError, match="Unknown provider: unknown"):
        function(*arguments)