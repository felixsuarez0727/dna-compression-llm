"""Lazy clients for hosted and local pattern-detection providers."""

from ..config import PROVIDER_DEFAULTS


def build_client(provider, api_key, save_log):
    if provider in ("deepseek", "chatgpt"):
        from openai import OpenAI

        if provider == "deepseek":
            return OpenAI(api_key=api_key, base_url="https://api.deepseek.com")
        return OpenAI(api_key=api_key)
    if provider == "gemini":
        from google import genai

        return genai.Client(api_key=api_key)
    if provider == "dnabert2":
        from .local_models import load_dnabert2

        return load_dnabert2(PROVIDER_DEFAULTS["dnabert2"], save_log)
    if provider == "hyenadna":
        from .local_models import load_hyenadna

        return load_hyenadna(PROVIDER_DEFAULTS["hyenadna"], save_log)
    raise ValueError(f"Unknown provider: {provider}")


def call_llm(client, provider, model_id, system_prompt, user_prompt):
    if provider in ("deepseek", "chatgpt"):
        response = client.chat.completions.create(
            model=model_id,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_prompt},
            ],
            temperature=0,
        )
        return response.choices[0].message.content
    if provider == "gemini":
        response = client.models.generate_content(
            model=model_id,
            contents=user_prompt,
        )
        return response.text
    raise ValueError(f"Unknown provider: {provider}")