from __future__ import annotations

import hashlib
import os
from pathlib import Path
import shutil
import subprocess

import pytest


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_failed_validation_still_finalizes_stable_external_hashes(tmp_path: Path) -> None:
    if shutil.which("shasum") is None:
        pytest.skip("validation script requires the macOS-compatible shasum command")

    root = Path(__file__).resolve().parents[1]
    script = root / "scripts" / "validate_macos_radia.sh"
    dummy_python = tmp_path / "python-without-pytest"
    dummy_python.write_text(
        """#!/usr/bin/env bash
if [[ ${1:-} == --version ]]; then
  printf 'Python dummy\\n'
  exit 0
fi
if [[ ${1:-} == -c ]]; then
  exit 1
fi
exit 1
""",
        encoding="utf-8",
    )
    dummy_python.chmod(0o755)
    output = tmp_path / "validation-output"

    # Preserve the host command-search path. The failure-path test is intended to
    # replace only the selected Python interpreter; reducing PATH to the shasum
    # directory prevents /usr/bin/env from locating bash on macOS (/bin/bash).
    env = os.environ.copy()
    env["PYTHON_BIN"] = str(dummy_python)
    env["AUTO_INSTALL_TEST_DEPS"] = "0"

    result = subprocess.run(
        [str(script), str(output)],
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 2
    log = output / "validation.log"
    log_sidecar = output / "validation.log.sha256"
    manifest = output / "SHA256SUMS.txt"
    manifest_sidecar = output / "SHA256SUMS.txt.sha256"
    assert all(path.is_file() for path in (log, log_sidecar, manifest, manifest_sidecar))
    assert "Validation stopped with exit status 2." in log.read_text(encoding="utf-8")

    expected_log_sidecar = f"{_sha256(log)}  validation.log\n"
    assert log_sidecar.read_text(encoding="utf-8") == expected_log_sidecar
    expected_manifest_sidecar = f"{_sha256(manifest)}  SHA256SUMS.txt\n"
    assert manifest_sidecar.read_text(encoding="utf-8") == expected_manifest_sidecar

    entries = {}
    for line in manifest.read_text(encoding="utf-8").splitlines():
        digest, size, name = line.split(maxsplit=2)
        entries[name] = (digest, int(size))
    assert entries["validation.log"] == (_sha256(log), log.stat().st_size)
    assert entries["validation.log.sha256"] == (
        _sha256(log_sidecar),
        log_sidecar.stat().st_size,
    )
