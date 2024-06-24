from dataclasses import dataclass

import thor_scsi.lib as ts


@dataclass
class index_class:
    # Configuration space coordinates.
    X, Y, Z = [
        ts.spatial_index.X,
        ts.spatial_index.Y,
        ts.spatial_index.Z
    ]

    # Phase-space coordinates.
    [x, px, y, py, ct, delta] = [
        ts.phase_space_index_internal.x,
        ts.phase_space_index_internal.px,
        ts.phase_space_index_internal.y,
        ts.phase_space_index_internal.py,
        ts.phase_space_index_internal.ct,
        ts.phase_space_index_internal.delta,
    ]

    
__all__ = [index_class]
