import numpy as np
import pytest

from apple_ii.field import field_integrals


def test_constant_field_integrals():
    z = np.linspace(0.0, 10.0, 101)
    rows = np.column_stack((z, np.full_like(z, 2.0), np.zeros_like(z), np.zeros_like(z)))
    result = field_integrals(rows)
    assert result.first_t_mm[0] == pytest.approx(20.0)
    assert result.second_t_mm2[0] == pytest.approx(100.0)


def test_field_integrals_require_increasing_finite_coordinates():
    with pytest.raises(ValueError, match="strictly increasing"):
        field_integrals([[1, 0, 0, 0], [0, 0, 0, 0]])
    with pytest.raises(ValueError, match="finite"):
        field_integrals([[0, 0, 0, 0], [1, float("nan"), 0, 0]])
