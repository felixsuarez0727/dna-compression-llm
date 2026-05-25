FROM python:3.10-slim

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        git \
        curl \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

RUN pip install --no-cache-dir \
    torch==2.2.2 \
    --index-url https://download.pytorch.org/whl/cpu

RUN pip install --no-cache-dir "numpy==1.26.4"

RUN pip install --no-cache-dir \
    "transformers==4.29.2" \
    "huggingface-hub==0.16.4" \
    "tokenizers==0.13.3" \
    "sentencepiece==0.2.0" \
    "einops==0.8.0" \
    "accelerate==0.27.2"

RUN pip install --no-cache-dir \
    "openai==1.30.1" \
    "google-genai==0.7.0"

RUN python -c "from huggingface_hub import snapshot_download; \
p = snapshot_download(repo_id='zhihan1996/DNABERT-2-117M'); \
print('Model Downloaded:', p)"

COPY pattern_detector.py .

VOLUME ["/data"]

ENTRYPOINT ["python", "pattern_detector.py"]
CMD ["--help"]