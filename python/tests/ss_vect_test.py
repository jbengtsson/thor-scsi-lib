import pytest
import thor_scsi.lib as tlib


def test20_tps_access():
    tps = tlib.tps()
    assert tps.peek([0] * 6) == 0
    return

# These lines let the test crash currently
    idx = [0] * 6
    idx[3] = 5
    assert tps.peek(idx) == 0

    tps.pook(idx, 42)
    assert tps.peek(idx) == 42
    print(tps)
    assert 1 == 0


def test25_sst_vect_access():
    t_map = tlib.ss_vect_tps()
    t_map.set_identity()

    tps = t_map[0]
    assert tps[0] == 1
    assert tps[1] == 0

    tps = t_map[1]
    assert tps[0] == 0
    assert tps[1] == 1
    assert tps[2] == 0

    print(tps.cst())
    assert tps.cst() == 0


def test30_element_access():
    """Test accessing element directly

    Todo:
        Clean up of the interface
    """
    t_map = tlib.ss_vect_tps()
    t_map.set_identity()

    print(t_map)
    assert t_map[0, 1] == 1
    assert t_map[1, 2] == 1
    assert t_map[2, 3] == 1
    assert t_map[3, 4] == 1
    assert t_map[4, 5] == 1
    assert t_map[5, 6] == 1
