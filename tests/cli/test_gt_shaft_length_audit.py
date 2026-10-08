"""The audit and crop share the actual GT disk-edge length definition."""

from ogo.cli.CheckFEModelBC import audit_femur_sideways
from ogo.fea.femur import DISTAL_FEMUR_NODE_SET, GREATER_TROCHANTER_NODE_SET


def test_shaft_length_uses_disk_edge_not_a_percentile_inside_the_disk():
    sets = {
        GREATER_TROCHANTER_NODE_SET: {
            "count": 100, "bounds": {"z_min": 40},
            "z_percentiles": {"p05": 44},
        },
        DISTAL_FEMUR_NODE_SET: {
            "count": 100, "z_percentiles": {"p50": 20},
            "planarity": {"rms_mm": 0},
        },
    }
    checks = audit_femur_sideways(sets, {}, 0.1)
    length = next(check for check in checks if check["name"] == "post-GT-support shaft length is positive")
    assert length["passed"]
    assert "gt_z_min_minus_distal_z_median=20 mm" in length["detail"]
