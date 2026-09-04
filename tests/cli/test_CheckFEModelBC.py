from ogo.cli.CheckFEModelBC import audit_spine_compression


def _flat_set(z):
    return {
        "count": 4,
        "bounds": {
            "x_min": 0.0,
            "x_max": 1.0,
            "y_min": 0.0,
            "y_max": 1.0,
            "z_min": float(z),
            "z_max": float(z),
        },
        "centroid": {"x": 0.5, "y": 0.5, "z": float(z)},
        "planarity": {"rms_mm": 0.0, "normal": [0.0, 0.0, 1.0]},
    }


def test_spine_bc_audit_requires_all_inferior_fixture_axes():
    node_sets = {
        "body_top": _flat_set(10.0),
        "body_bottom": _flat_set(0.0),
    }
    constraints = {
        "top_displacement": {"axes": ["z"], "values": [-0.1]},
        "bottom_fixed_x": {"axes": ["x"], "values": [0.0]},
        "bottom_fixed_y": {"axes": ["y"], "values": [0.0]},
        "bottom_fixed_z": {"axes": ["z"], "values": [0.0]},
    }

    checks = audit_spine_compression(node_sets, constraints, tolerance=1.0e-5)
    by_name = {check["name"]: check for check in checks}

    assert by_name["bottom_fixed_x fixes x"]["passed"] is True
    assert by_name["bottom_fixed_y fixes y"]["passed"] is True
    assert by_name["bottom_fixed_z fixes z"]["passed"] is True
