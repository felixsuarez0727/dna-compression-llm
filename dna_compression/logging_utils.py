"""Logging helpers shared by command modules."""

from .config import get_logs_dir


def save_log(line, filename=None):
    safe_line = line.replace("→", "->")
    print(safe_line)

    if filename:
        logs_dir = get_logs_dir()
        logs_dir.mkdir(parents=True, exist_ok=True)
        with open(logs_dir / filename, "a", encoding="utf-8") as log_file:
            log_file.write(safe_line + "\n")