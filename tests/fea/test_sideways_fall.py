import ast
import inspect
import textwrap

import pytest


def test_sideways_fall_icp_path_does_not_pre_rotate_input_images():
    pytest.importorskip("vtk")
    pytest.importorskip("vtkbone")

    from ogo.fea import femur

    source = textwrap.dedent(inspect.getsource(femur.sidewaysFallFe))
    tree = ast.parse(source)
    pre_rotation_calls = [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Attribute) and node.attr == "preRotateImage"
    ]

    assert pre_rotation_calls == []


def test_sideways_fall_icp_path_uses_voxel_surface_points():
    pytest.importorskip("vtk")
    pytest.importorskip("vtkbone")

    from ogo.fea import femur

    source = textwrap.dedent(inspect.getsource(femur.sidewaysFallFe))
    tree = ast.parse(source)
    called_names = {
        node.func.id
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }
    marching_cube_calls = [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Attribute) and node.attr == "marchingCubes"
    ]

    assert "surface_points_from_vtk_mask" in called_names
    assert "scale_reference_point_cloud_to_sample" in called_names
    assert marching_cube_calls == []


def test_sideways_fall_icp_path_uses_deterministic_point_cloud_transform():
    pytest.importorskip("vtk")
    pytest.importorskip("vtkbone")

    from ogo.fea import femur

    source = textwrap.dedent(inspect.getsource(femur.sidewaysFallFe))
    tree = ast.parse(source)
    called_names = {
        node.func.id
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }
    helper_calls = [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Attribute) and node.attr == "iterativeClosestPoint"
    ]

    assert "estimate_rigid_icp" in called_names
    assert helper_calls == []


def test_sideways_fall_distal_support_is_bbox_relative_recipe_plane():
    pytest.importorskip("vtk")
    pytest.importorskip("vtkbone")

    from ogo.fea import femur

    bounds = (
        -76.45140061071066,
        -31.451400610710657,
        -102.34001525211478,
        -36.34001525211478,
        -0.39488812795613626,
        83.60511187204386,
    )
    plane = femur.bbox_relative_oriented_contact_plane(
        bounds,
        center_fraction=femur.DISTAL_SHAFT_FIXTURE_CENTER_FRACTION,
        size_fraction=femur.DISTAL_SHAFT_FIXTURE_SIZE_FRACTION,
        normal=femur.DISTAL_SHAFT_FIXTURE_NORMAL,
        size_bounds_axes=("x", "y"),
        shape="anatomy",
    )

    assert femur.DISTAL_SHAFT_FIXTURE_CENTER_FRACTION == (0.5, 0.5, -0.1)
    assert femur.DISTAL_SHAFT_FIXTURE_SIZE_FRACTION == (1.0, 1.0)
    assert plane["center"] == pytest.approx(
        (-53.95140061071066, -69.34001525211478, -8.794888127956137)
    )
    assert plane["normal"] == pytest.approx(
        femur.DISTAL_SHAFT_FIXTURE_NORMAL,
        abs=1.0e-12,
    )
    assert plane["size"] == pytest.approx((63.750, 45.0), abs=0.2)


def test_sideways_fall_fits_contact_canvas_before_generating_pmma():
    pytest.importorskip("vtk")
    pytest.importorskip("vtkbone")

    from ogo.fea import femur

    source = textwrap.dedent(inspect.getsource(femur.sidewaysFallFe))
    tree = ast.parse(source)
    called_names = {
        node.func.id
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }

    assert "projected_material_disk_required_bounds" in called_names
    assert "fit_vtk_images_to_physical_bounds" in called_names
