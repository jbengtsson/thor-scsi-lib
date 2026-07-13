from __future__ import annotations

import importlib
import math
from types import ModuleType
from typing import Iterable


def require_radia(radia_module: ModuleType | object | None = None):
    """Return an injected RADIA-like module or import the real extension lazily."""
    if radia_module is not None:
        return radia_module
    try:
        return importlib.import_module("radia")
    except ImportError as exc:
        raise ImportError(
            "RADIA is required for geometry or field evaluation; verify "
            "`python -c \"import radia; print(radia.__file__)\"`."
        ) from exc


def finite_number(name: str, value: float) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def positive_number(name: str, value: float) -> float:
    value = finite_number(name, value)
    if value <= 0:
        raise ValueError(f"{name} must be > 0")
    return value


def nonnegative_number(name: str, value: float) -> float:
    value = finite_number(name, value)
    if value < 0:
        raise ValueError(f"{name} must be >= 0")
    return value


def finite_sequence(name: str, values: Iterable[float]) -> tuple[float, ...]:
    result = tuple(float(value) for value in values)
    if not all(math.isfinite(value) for value in result):
        raise ValueError(f"{name} must contain only finite values")
    return result
