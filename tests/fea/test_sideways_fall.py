import ast
import inspect
import textwrap

import pytest


def test_sideways_fall_icp_path_does_not_pre_rotate_input_images():
    pytest.importorskip("vtk")
    pytest.importorskip("vtkbone")

    from ogo.fea import sideways_fall

    source = textwrap.dedent(inspect.getsource(sideways_fall.sidewaysFallFe))
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

    from ogo.fea import sideways_fall

    source = textwrap.dedent(inspect.getsource(sideways_fall.sidewaysFallFe))
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
