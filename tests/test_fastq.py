import gzip
import random

import pytest

from dna_compression.fastq import (
    extract_sequences,
    generate_fastq_file,
    generate_sequence,
)


def test_generate_sequence_creates_a_valid_fastq_record():
    random.seed(7)

    record = generate_sequence(3, target_length=31)
    header, sequence, separator, quality = record.splitlines()

    assert header == "@SEQ_ID_3"
    assert set(sequence) <= set("ACGTN")
    assert len(sequence) == 31
    assert separator == "+"
    assert len(quality) == len(sequence)


def test_generate_fastq_file_creates_requested_records(tmp_path):
    random.seed(11)
    output_file = tmp_path / "synthetic.fastq"

    generate_fastq_file(output_file, num_sequences=3, seq_length=20)

    lines = output_file.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 12
    assert lines[0] == "@SEQ_ID_1"
    assert lines[4] == "@SEQ_ID_2"
    assert all(len(lines[index]) == 20 for index in (1, 5, 9))


@pytest.mark.parametrize("compressed", [False, True])
def test_extract_sequences_supports_plain_and_gzip_fastq(tmp_path, compressed):
    input_file = tmp_path / ("reads.fastq.gz" if compressed else "reads.fastq")
    output_file = tmp_path / "sequences.seq.txt"
    content = "@read-1\nACGT\n+\nIIII\n@read-2\nTGCA\n+\nJJJJ\n"

    if compressed:
        with gzip.open(input_file, "wt", encoding="utf-8") as output_handle:
            output_handle.write(content)
    else:
        input_file.write_text(content, encoding="utf-8")

    extract_sequences(input_file, output_file, n=1)

    assert output_file.read_text(encoding="utf-8") == "ACGT\n"