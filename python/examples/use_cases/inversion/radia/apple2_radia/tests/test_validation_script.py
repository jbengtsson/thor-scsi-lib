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
    sizes = output / "BYTES.txt"
    manifest = output / "SHA256SUMS.txt"
    manifest_sidecar = output / "SHA256SUMS.txt.sha256"
    assert all(path.is_file() for path in (log, log_sidecar, sizes, manifest, manifest_sidecar))
    assert "Validation stopped with exit status 2." in log.read_text(encoding="utf-8")

    expected_log_sidecar = f"{_sha256(log)}  validation.log\n"
    assert log_sidecar.read_text(encoding="utf-8") == expected_log_sidecar
    expected_manifest_sidecar = f"{_sha256(manifest)}  SHA256SUMS.txt\n"
    assert manifest_sidecar.read_text(encoding="utf-8") == expected_manifest_sidecar

    # SHA256SUMS.txt must be directly consumable by macOS shasum -c.
    check = subprocess.run(
        ["shasum", "-a", "256", "-c", manifest.name],
        cwd=output,
        text=True,
        capture_output=True,
        check=False,
    )
    assert check.returncode == 0, check.stdout + check.stderr

    hashes = {}
    for line in manifest.read_text(encoding="utf-8").splitlines():
        digest, name = line.split("  ", maxsplit=1)
        assert len(digest) == 64
        hashes[name] = digest
    assert hashes["validation.log"] == _sha256(log)
    assert hashes["validation.log.sha256"] == _sha256(log_sidecar)
    assert hashes["BYTES.txt"] == _sha256(sizes)

    byte_lengths = {}
    for line in sizes.read_text(encoding="utf-8").splitlines():
        size, name = line.split(maxsplit=1)
        byte_lengths[name] = int(size)
    assert byte_lengths["validation.log"] == log.stat().st_size
    assert byte_lengths["validation.log.sha256"] == log_sidecar.stat().st_size
