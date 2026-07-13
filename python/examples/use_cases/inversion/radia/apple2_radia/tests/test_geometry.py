import pytest

from apple2.geometry import Apple2Parameters, build_device
from fake_radia import FakeRadia


def test_end_blocks_are_contiguous_by_default():
    radia = FakeRadia()
    parameters = Apple2Parameters(n_periods=1, end_block_fraction=0.5, end_clearance_mm=0.0)
    device = build_device(parameters, radia_module=radia)
    first_row_id = device.rows["UL"].object_id
    block_ids = radia.containers[first_row_id]
    centers = [radia.blocks[block_id]["center"][2] for block_id in block_ids]
    dimensions = [radia.blocks[block_id]["dimensions"][2] for block_id in block_ids]
    assert len(block_ids) == 6
    assert centers[-2] + dimensions[-2] / 2 == pytest.approx(centers[0] - dimensions[0] / 2)
    assert centers[-1] - dimensions[-1] / 2 == pytest.approx(centers[3] + dimensions[3] / 2)


def test_end_clearance_is_explicit():
    radia = FakeRadia()
    parameters = Apple2Parameters(n_periods=1, end_clearance_mm=1.25)
    device = build_device(parameters, radia_module=radia)
    block_ids = radia.containers[device.rows["UL"].object_id]
    lower_end = radia.blocks[block_ids[-2]]
    first = radia.blocks[block_ids[0]]
    gap = first["center"][2] - first["dimensions"][2] / 2 - (
        lower_end["center"][2] + lower_end["dimensions"][2] / 2
    )
    assert gap == pytest.approx(1.25)


def test_phase_changes_translate_existing_rows_without_rebuild():
    radia = FakeRadia()
    device = build_device(Apple2Parameters(n_periods=1), radia_module=radia)
    object_count = radia.next_id
    row_ids = device.row_ids()
    device.set_phase(10.0)
    device.set_phase(12.0)
    assert radia.next_id == object_count
    assert device.row_ids() == row_ids
    ur_translations = [
        transform[1][2]
        for object_id, transform in radia.transforms
        if object_id == row_ids["UR"]
    ]
    assert ur_translations == pytest.approx([10.0, 2.0])
    assert device.phase_mm == pytest.approx(12.0)
    assert device.parameters.phase_mm == pytest.approx(12.0)


def test_parameter_validation_rejects_nonfinite_and_invalid_values():
    with pytest.raises(ValueError, match="period_mm"):
        Apple2Parameters(period_mm=float("nan")).validate()
    with pytest.raises(ValueError, match="blocks_per_period"):
        Apple2Parameters(blocks_per_period=8).validate()
    with pytest.raises(ValueError, match="end_block_fraction"):
        Apple2Parameters(end_block_fraction=1.1).validate()
    with pytest.raises(ValueError, match="end_magnetization_fraction"):
        Apple2Parameters(end_magnetization_fraction=1.1).validate()
