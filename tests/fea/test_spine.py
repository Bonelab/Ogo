from tests.fea import helpers
from ogo.fea import spine
import numpy as np
import pytest
from types import SimpleNamespace


def _mock_registration_inputs(monkeypatch):
    points = np.array([[0, 0, 0], [40, 0, 0], [0, 35, 0], [0, 0, 30]], dtype=float)
    monkeypatch.setattr(spine, "surface_points_from_vtk_mask", lambda *a, **k: points)
    monkeypatch.setattr(spine.ogo, "readPolyData", lambda *a: None)
    monkeypatch.setattr(spine, "polydata_points", lambda *a: points)
    monkeypatch.setattr(spine, "sample_points", lambda *a, **k: points)
    return SimpleNamespace(GetOutput=lambda: None)


def test_preparation_and_registration_preserve_spine_labels(tmp_path, monkeypatch):
    vtk = pytest.importorskip("vtk")
    from vtk.util.numpy_support import numpy_to_vtk
    from ogo.fea.boundary import vtk_image_to_numpy

    def image(array):
        result = vtk.vtkImageData()
        result.SetDimensions(array.shape)
        result.SetSpacing(1, 1, 1)
        result.GetPointData().SetScalars(numpy_to_vtk(array.ravel(order="F"), deep=True))
        return result

    labels = np.zeros((30, 40, 30), dtype=np.uint8)
    labels[10:20, 10:20, 10:20] = 2
    labels[13:17, 20:30, 13:17] = 3
    density = image(np.full(labels.shape, 200, dtype=np.float32))
    mask = image(labels)
    monkeypatch.setattr(spine, "read", lambda path: SimpleNamespace(
        GetOutput=lambda: density if path == "density" else mask))
    prepared = spine._prepare_spine_inputs(
        "density", "mask", 2, 3, None, None, 1.0, 0.0, 0.5, "body")
    assert prepared[6]["body_cleanup_removed_voxels"] == 0
    reference = tmp_path / "reference.vtk"
    reference.touch()
    identity = vtk.vtkTransform()
    monkeypatch.setattr(spine, "get_icp_with_scaling", lambda *a, **k: identity)
    registered = spine._register_spine_inputs(
        *prepared[:6], str(reference), 1.0, None, "0.8,0.8,0.75",
        "1.2,1.2,1.3", "numpy", 8000, 50, 4, 2.0, True, "body")
    assert np.count_nonzero(vtk_image_to_numpy(registered[0])) == np.count_nonzero(labels == 2)
    assert np.count_nonzero(vtk_image_to_numpy(registered[1])) == np.count_nonzero(labels == 3)
    assert registered[3] is None
    np.testing.assert_array_equal(registered[4], np.eye(3))


def test_spine_icp_preserves_native_orientation_at_initialization(monkeypatch):
    body = _mock_registration_inputs(monkeypatch)
    calls = []

    def estimate(**kwargs):
        calls.append(kwargs)
        return dict(rotation=np.eye(3), translation=np.zeros(3), iterations=1, mean_distance=0)

    monkeypatch.setattr(spine, "estimate_rigid_icp", estimate)
    spine.get_icp_with_scaling(body, "reference.vtk")
    assert len(calls) == 1
    assert calls[0]["start_by_matching_centroids_only"] is True


def test_spine_icp_rejects_invalid_process_without_axis_swapped_retries(monkeypatch):
    body = _mock_registration_inputs(monkeypatch)
    calls = []

    def estimate(**kwargs):
        calls.append(kwargs)
        return dict(rotation=np.diag([-1., -1., 1.]), translation=np.zeros(3), iterations=1, mean_distance=0)

    monkeypatch.setattr(spine, "estimate_rigid_icp", estimate)
    centroids = iter([np.zeros(3), np.array([0., 30., 0.])])
    monkeypatch.setattr(spine, "_mask_centroid_physical", lambda *a: next(centroids))
    monkeypatch.setattr(spine, "_icp_process_orientation_metrics", lambda *a: dict(
        status="failed", axial_offset_mm=30, transverse_offset_mm=5))
    with pytest.raises(ValueError, match="orientation"):
        spine.get_icp_with_scaling(body, "reference.vtk", process=body, scale_factors=1)
    assert len(calls) == 1


@pytest.mark.parametrize("rotation", [
    np.diag([1., -1., -1.]),
    np.diag([-1., -1., 1.]),
    np.array([[1., 0., 0.], [0., 0., -1.], [0., 1., 0.]]),
])
def test_spine_registration_rejects_anatomical_axis_flips(rotation):
    with pytest.raises(ValueError, match="orientation"):
        spine.check_spine_registration_orientation(rotation, [0, 0, 0], [0, 30, 0])


def test_spine_registration_accepts_small_anatomical_rotation():
    angle = np.deg2rad(32)
    rotation = np.array([[1., 0., 0.], [0., np.cos(angle), -np.sin(angle)],
                         [0., np.sin(angle), np.cos(angle)]])
    spine.check_spine_registration_orientation(rotation, [0, 0, 0], [0, 30, 0])


def test_default_spine_cap_geometry_matches_maintained_settings():
    assert spine.DEFAULT_SPINE_PMMA_THICKNESS_MM == 10
    assert spine.DEFAULT_SPINE_PMMA_INTRUSION_MM == 6
    assert spine.DEFAULT_SPINE_REGISTRATION_BACKEND == "numpy"
    assert spine.DEFAULT_SPINE_REGISTRATION_LANDMARKS == 8000
    assert spine.DEFAULT_SPINE_REGISTRATION_ITERATIONS == 50
    assert spine.DEFAULT_SPINE_ICP_TARGET == "body"


def test_spine_reference_path_depends_on_icp_target():
    assert spine.default_spine_reference_path("body").name == "L4_BODY_SPINE_COMPRESSION_REF.vtk"
    assert spine.default_spine_reference_path("vertebra").name == "L4_FULL_VERTEBRA_SPINE_COMPRESSION_REF.vtk"


def test_benchmark_presets_match_spinefe_notebook_settings():
    linear = helpers.benchmark_linear_params()
    nonlinear = helpers.benchmark_nonlinear_params()

    assert linear["fe_displacement"] == -0.2
    assert linear["target_displacement_percent"] == 0.68
    assert linear["elastic_E_func"] == "kopperdahl_trab_E"
    assert linear["yield_comp_func"] is None
    assert linear["pmma_yield_compression"] is None
    assert nonlinear["fe_displacement"] == -2.0
    assert nonlinear["target_displacement_percent"] == 4.0
    assert nonlinear["yield_comp_func"] == "kopperdahl_trab_yc"
    assert nonlinear["pmma_yield_compression"] == 70.0


def test_process_vertebra_rejects_unknown_settings_before_reading_inputs():
    with pytest.raises(TypeError, match="registration_landmakrs"):
        spine.process_vertebra("missing_mask", "missing_image", "model.n88model", 2, 3,
                               "reference.vtk", registration_landmakrs=8000)
