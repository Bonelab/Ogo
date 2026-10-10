"""Shaft length is measured from the generated GT disk, not its fixture box."""

import numpy as np
import pytest

from ogo.fea.femur import crop_vtk_images_to_greater_trochanter_length
from ogo.util.vtk_image import vtk_image_to_numpy


def vtk_image(data, *, origin=(0, 0, 0), spacing=(1, 1, 1)):
    import vtk
    from vtk.util.numpy_support import numpy_to_vtk

    image = vtk.vtkImageData()
    image.SetDimensions(data.shape)
    image.SetOrigin(*origin)
    image.SetSpacing(*spacing)
    image.GetPointData().SetScalars(numpy_to_vtk(data.ravel(order="F"), deep=True))
    return image


def images(distal_index=0, spacing=1.0):
    bone = np.zeros((12, 12, 70), dtype=np.uint16)
    bone[4:8, 4:8, distal_index:65] = 1
    disk = np.zeros_like(bone)
    disk[3:9, 0:4, 40:46] = 5000
    geometry = {"origin": (2.0, 3.0, 7.0), "spacing": (spacing,) * 3}
    return (
        vtk_image(np.maximum(bone, disk), **geometry),
        vtk_image(disk, **geometry),
    )


def flat_boundary(image):
    """Known complete flat faces for synthetic crop-only tests."""
    data = vtk_image_to_numpy(image)
    index = int(np.argwhere(data != 0)[:, 2].min())
    z = image.GetOrigin()[2] + (image.GetExtent()[4] + index - .5) * image.GetSpacing()[2]
    return {"required_flat_face_z_mm": z}


@pytest.mark.parametrize("spacing", [0.5, 1.0, 2.0])
def test_length_uses_generated_disk_and_voxel_faces(spacing):
    material, disk = images(spacing=spacing)
    cropped, face, meta = crop_vtk_images_to_greater_trochanter_length(
        [material, disk], material, gt_support_vtk=disk, retained_length_mm=20,
        distal_boundary=flat_boundary(material),
    )
    expected_gt_edge = 7.0 + (40 - 0.5) * spacing
    assert meta["greater_trochanter_support_distal_edge_z"] == pytest.approx(expected_gt_edge)
    assert meta["cut_z_mm"] == pytest.approx(expected_gt_edge - 20)
    face_data = vtk_image_to_numpy(face)
    assert np.count_nonzero(face_data) == 16
    assert np.all(np.argwhere(face_data)[:, 2] == 0)
    assert cropped[0].GetOrigin()[2] - spacing / 2 == pytest.approx(expected_gt_edge - 20)
    assert np.count_nonzero(vtk_image_to_numpy(cropped[1])) == 6 * 4 * 6


def test_exactly_twenty_mm_coverage_still_has_a_distal_face():
    material, disk = images(distal_index=20)
    _, face, meta = crop_vtk_images_to_greater_trochanter_length(
        [material, disk], material, gt_support_vtk=disk, retained_length_mm=20,
        distal_boundary=flat_boundary(material),
    )
    assert meta["available_below_gt_support_distal_edge_mm"] == pytest.approx(20)
    assert np.count_nonzero(vtk_image_to_numpy(face)) == 16


def test_insufficient_coverage_is_not_solved_as_a_shorter_model():
    material, disk = images(distal_index=21)
    with pytest.raises(ValueError, match="too short"):
        crop_vtk_images_to_greater_trochanter_length(
            [material, disk], material, gt_support_vtk=disk, retained_length_mm=20,
            distal_boundary=flat_boundary(material),
        )


def test_empty_gt_support_cannot_define_length():
    material, disk = images()
    disk = vtk_image(np.zeros((12, 12, 70), dtype=np.uint16), origin=(2, 3, 7))
    with pytest.raises(ValueError, match="empty GT"):
        crop_vtk_images_to_greater_trochanter_length(
            [material, disk], material, gt_support_vtk=disk, retained_length_mm=20,
            distal_boundary=flat_boundary(material),
        )


def test_contact_canvas_extent_preserves_absolute_model_coordinates():
    material, disk = images()
    for image in (material, disk):
        image.SetExtent(10, 21, 20, 31, 30, 99)
    cropped, _, meta = crop_vtk_images_to_greater_trochanter_length(
        [material, disk], material, gt_support_vtk=disk, retained_length_mm=20,
        distal_boundary=flat_boundary(material),
    )
    assert meta["cut_z_mm"] == pytest.approx(7 + 30 + 40 - 0.5 - 20)
    assert cropped[0].GetOrigin()[2] - 0.5 == pytest.approx(meta["cut_z_mm"])


def test_partial_distal_face_is_rejected_even_with_twenty_mm_raw_span():
    material, disk = images()
    with pytest.raises(ValueError, match="too short"):
        crop_vtk_images_to_greater_trochanter_length(
            [material, disk], material, gt_support_vtk=disk, retained_length_mm=20,
            distal_boundary={"required_flat_face_z_mm": 30.5},
        )


def test_eligible_crop_preserves_disks_and_clears_original_end():
    material, disk = images()
    _, _, meta = crop_vtk_images_to_greater_trochanter_length(
        [material, disk], material, gt_support_vtk=disk, retained_length_mm=20,
        distal_boundary={"required_flat_face_z_mm": 19.2},
    )
    assert meta["retained_length_mm"] == 20
    assert meta["available_shaft_length_mm"] == 27
    assert meta["cut_z_mm"] >= meta["safe_distal_face_z_mm"]
