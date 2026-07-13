from __future__ import annotations

from pathlib import Path
import re

import apple2


def test_public_version_matches_project_metadata() -> None:
    pyproject = (Path(__file__).resolve().parents[1] / "pyproject.toml").read_text(
        encoding="utf-8"
    )
    match = re.search(r'^version\s*=\s*"([^"]+)"', pyproject, flags=re.MULTILINE)
    assert match is not None
    assert apple2.__version__ == match.group(1)
