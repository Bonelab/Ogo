"""The sidecar and saved BC nodes use the same generated-disk origin."""

from types import SimpleNamespace

import numpy as np
import pytest
import vtk
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray

from tests.fea.test_gt_disk_crop import images, vtk_image


def test_short_scan_measurement_is_available_without_cropping():
    from ogo.fea.shaft_geometry import measure_available_shaft

    material, disk = images(distal_index=25)
    result = measure_available_shaft(material, disk, requested_length_mm=20,
                                     distal_boundary={"required_flat_face_z_mm": 31.5})
    assert result["available_shaft_length_mm"] == 15
    assert not result["eligible_for_requested_length"]
    assert result["status"] == "too_short"
    assert result["measured_shaft_length_mm"] is None


def model_with_nodes(distal_z):
    xyz = np.array([[0, 0, 40.], [1, 0, 40.], [0, 0, distal_z], [1, 0, distal_z]])
    points = vtk.vtkPoints()
    points.SetData(numpy_to_vtk(xyz))
    sets = {"Greater_Trochanter_PMMA_Nodes": numpy_to_vtkIdTypeArray(np.array([0, 1])),
            "Distal_Femur_Nodes": numpy_to_vtkIdTypeArray(np.array([2, 3]))}
    return SimpleNamespace(GetPoints=lambda: points, GetNodeSet=lambda name: sets[name])


def test_final_measurement_uses_gt_minimum_and_distal_median():
    from ogo.fea.shaft_geometry import verify_model_shaft

    meta = {"requested_shaft_length_mm": 20., "retained_length_mm": 20.,
            "distal_length_origin_z": 40., "cut_z_mm": 20.,
            "safe_distal_face_z_mm": 19.}
    verified = verify_model_shaft(model_with_nodes(20), meta)
    assert verified["measured_shaft_length_mm"] == 20
    assert verified["status"] == "verified"
    with pytest.raises(ValueError, match="shaft"):
        verify_model_shaft(model_with_nodes(19), meta)


def test_oblique_end_reduces_usable_coverage_not_the_requested_length():
    from ogo.fea.shaft_geometry import measure_available_shaft

    material, disk = images()
    result = measure_available_shaft(
        material, disk, requested_length_mm=20,
        distal_boundary={"required_flat_face_z_mm": 28.2},
    )
    assert result["raw_available_shaft_length_mm"] == 40
    assert result["safe_distal_face_z_mm"] == pytest.approx(28.5)
    assert result["available_shaft_length_mm"] == 18
    assert not result["eligible_for_requested_length"]
    assert result["requested_shaft_length_mm"] == 20
    assert result["status"] == "too_short"


def test_source_boundary_is_required_for_coverage():
    from ogo.fea.shaft_geometry import measure_available_shaft

    material, disk = images()
    with pytest.raises(ValueError, match="distal.*boundary"):
        measure_available_shaft(material, disk, requested_length_mm=20)


def test_capture_uses_full_native_voxel_faces_and_nonzero_extent():
    from ogo.fea.shaft_geometry import capture_distal_scan_face

    data = np.zeros((5, 6, 8), dtype=np.uint8)
    data[1:3, 2:4, 2:] = 1
    data[4, 5, 0] = 2  # Another label must not define this femur's end.
    image = vtk_image(data, origin=(10, 20, 30), spacing=(.5, .75, 2))
    image.SetExtent(5, 9, 10, 15, 20, 27)
    face = capture_distal_scan_face(image, labels={1})
    points = face["points_xyz"]
    assert points[:, 0].min() == pytest.approx(12.75)
    assert points[:, 0].max() == pytest.approx(13.75)
    assert points[:, 1].min() == pytest.approx(28.625)
    assert points[:, 1].max() == pytest.approx(30.125)
    assert np.all(points[:, 2] == 73)
    assert face["native_distal_slice_index"] == 22
    assert face["native_face_voxels"] == 4


def test_boundary_uses_inverse_icp_and_resampling_clearance():
    from ogo.fea.shaft_geometry import aligned_distal_scan_boundary

    angle = np.deg2rad(30)
    matrix = np.eye(4)
    matrix[:3, :3] = [[1, 0, 0], [0, np.cos(angle), -np.sin(angle)],
                        [0, np.sin(angle), np.cos(angle)]]
    matrix[:3, 3] = [10, 20, 30]
    points = np.array([[0, 0, 5], [0, 8, 5], [4, 8, 5], [4, 0, 5]], float)
    boundary = aligned_distal_scan_boundary(
        {"points_xyz": points, "native_face_voxels": 1}, matrix,
        resampling_spacing=(1, 1, 1), output_spacing=(1, 1, 1),
    )
    inverse = np.linalg.inv(matrix)
    transformed = points @ inverse[:3, :3].T + inverse[:3, 3]
    clearance = .5 * np.abs(inverse[2, :3]).sum() + .5
    assert boundary["aligned_face_z_min_mm"] == pytest.approx(transformed[:, 2].min())
    assert boundary["aligned_face_z_max_mm"] == pytest.approx(transformed[:, 2].max())
    assert boundary["resampling_clearance_mm"] == pytest.approx(clearance)
    assert boundary["required_flat_face_z_mm"] == pytest.approx(transformed[:, 2].max() + clearance)


def test_planar_bc_nodes_alone_do_not_verify_complete_section():
    from ogo.fea.shaft_geometry import verify_model_shaft

    meta = {"retained_length_mm": 20., "distal_length_origin_z": 40.,
            "cut_z_mm": 20., "safe_distal_face_z_mm": 21.}
    with pytest.raises(ValueError, match="original.*boundary"):
        verify_model_shaft(model_with_nodes(20), meta)
