"""Compatibility comparison for restored sequence files."""


def sanitized_files_are_equal(file1, file2):
    with open(file1, "r", encoding="utf-8") as first_handle, open(
        file2, "r", encoding="utf-8"
    ) as second_handle:
        lines1 = [line.strip() for line in first_handle if line.strip()]
        lines2 = [line.strip() for line in second_handle if line.strip()]
        return lines1 == lines2