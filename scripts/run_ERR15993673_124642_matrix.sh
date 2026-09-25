#!/usr/bin/env bash
# Full experiment matrix: ERR15993673 large dataset (124,642 reads)
# Providers: gemini, chatgpt, deepseek | Overhead: 5 and 101
#
# Required env vars (do NOT commit keys):
#   export GEMINI_API_KEY="..."
#   export OPENAI_API_KEY="..."
#   export DEEPSEEK_API_KEY="..."

set -euo pipefail
cd "$(dirname "$0")/.."

RUN_DIR="data/runs/ERR15993673_124642_20260522"
INPUT="data/input/ERR15993673_249284.seq.txt"
OUT_DIR="${RUN_DIR}/outputs"
LOG_DIR="${RUN_DIR}/logs"
MANIFEST="${RUN_DIR}/run_manifest.txt"

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

if [[ ! -f "${INPUT}" ]]; then
  echo "Missing input: ${INPUT}"
  exit 1
fi

for v in GEMINI_API_KEY OPENAI_API_KEY DEEPSEEK_API_KEY; do
  if [[ -z "${!v:-}" ]]; then
    echo "Missing env var: ${v}"
    exit 1
  fi
done

source venv/bin/activate

log_step() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${MANIFEST}"
}

run_detector() {
  local provider="$1"
  local overhead="$2"
  local key_var="$3"
  local model="$4"
  local threads="$5"
  local tag="delta${overhead}"
  local out_json="${OUT_DIR}/${provider}_${tag}_patterns.json"

  log_step "START detector ${provider} overhead=${overhead}"
  PYTHONUNBUFFERED=1 python pattern_detector.py \
    -f "${INPUT}" \
    -o "${out_json}" \
    -b 80 \
    -p "${provider}" \
    -m "${model}" \
    --threads "${threads}" \
    --overhead "${overhead}" \
    --max_retries 4 \
    --retry_base_delay 2.0 \
    --log-dir "${LOG_DIR}" \
    -k "${!key_var}"
  log_step "DONE detector ${provider} overhead=${overhead}"
}

run_compress_pipeline() {
  local provider="$1"
  local overhead="$2"
  local tag="delta${overhead}"
  local prefix="ERR15993673_124642_${provider}_${tag}"
  local patterns="${OUT_DIR}/${provider}_${tag}_patterns.json"
  local compressed="${OUT_DIR}/${prefix}.seq.compress"
  local restored="${OUT_DIR}/${prefix}.seq.restored"

  log_step "START compress ${provider} overhead=${overhead}"
  python fastq_compressor.py \
    -i "${INPUT}" \
    -o "${compressed}" \
    -p "${patterns}"
  python fastq_decompressor.py \
    -i "${compressed}" \
    -o "${restored}" \
    -p "${patterns}"
  python check_files.py "${INPUT}" "${restored}" | tee -a "${MANIFEST}"
  log_step "DONE compress ${provider} overhead=${overhead}"
}

log_step "RUN_DIR=${RUN_DIR}"
log_step "INPUT=${INPUT} ($(wc -l < "${INPUT}") lines)"

for overhead in 5 101; do
  run_detector "deepseek" "${overhead}" "DEEPSEEK_API_KEY" "deepseek-chat" 4
  run_compress_pipeline "deepseek" "${overhead}"

  run_detector "chatgpt" "${overhead}" "OPENAI_API_KEY" "gpt-4o-mini" 4
  run_compress_pipeline "chatgpt" "${overhead}"

  run_detector "gemini" "${overhead}" "GEMINI_API_KEY" "gemini-2.5-flash-lite" 2
  run_compress_pipeline "gemini" "${overhead}"
done

log_step "ALL EXPERIMENTS COMPLETED"
