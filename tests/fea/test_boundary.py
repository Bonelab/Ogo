import pytest

np = pytest.importorskip("numpy")

from ogo.fea.boundary import (  # noqa: E402
    bbox_relative_contact_bounds,
    bbox_relative_contact_direction,
    bbox_relative_contact_plane,
    bounds_with_reference_extent,
    foreground_voxel_center_bounds_from_mask,
    generate_bone_cap_mask,
    generate_projected_material_disk_mask,
    pad_vtk_images_to_physical_bounds,
    projected_material_disk_required_bounds,
    smooth_label_mask_vtk,
)


def test_bbox_relative_contact_plane_uses_authored_axis_convention():
    bounds = (10, 50, -20, 80, 100, 200)

    low_y_plane = bbox_relative_contact_plane(
        bounds,
        center_fraction=(0.5, 0.02, 0.25),
        size_fraction=(1.1, 0.5),
        projection_axis="y",
        shape="rectangle",
    )
    high_y_plane = bbox_relative_contact_plane(
        bounds,
        center_fraction=(0.5, 1.02, 0.25),
        size_fraction=(1.1, 0.5),
        projection_axis="y",
        shape="rectangle",
    )

    assert bbox_relative_contact_direction((0.5, 0.02, 0.25), projection_axis="y") == "down"
    assert bbox_relative_contact_direction((0.5, 1.02, 0.25), projection_axis="y") == "up"
    assert low_y_plane["center"] == pytest.approx((30, -18, 125))
    assert low_y_plane["normal"] == pytest.approx((0, 1, 0))
    assert low_y_plane["u_axis"] == pytest.approx((0, 0, -1))
    assert low_y_plane["v_axis"] == pytest.approx((-1, 0, 0))
    assert low_y_plane["size"] == pytest.approx((110, 20))
    assert low_y_plane["shape"] == "rectangle"
    assert high_y_plane["normal"] == pytest.approx((0, -1, 0))
    assert high_y_plane["u_axis"] == pytest.approx((0, 0, 1))
    assert high_y_plane["v_axis"] == pytest.approx((-1, 0, 0))


def test_bbox_relative_contact_bounds_can_enforce_square_footprint():
    bounds = (10, 50, -20, 80, 100, 200)

    contact_bounds = bbox_relative_contact_bounds(
        bounds,
        center_fraction=(0.5, 1.1, 0.25),
        size_fraction=(1.1, 0.5),
        projection_axis="y",
        shape="square",
    )

    assert contact_bounds == pytest.approx((8, 52, -20, 80, 103, 147))


def test_bounds_with_reference_extent_preserves_current_center():
    current_bounds = (10, 20, -5, 35, 100, 133)
    reference_bounds = (0, 49, 0, 43, 0, 34)

    bounds = bounds_with_reference_extent(current_bounds, reference_bounds)

    assert bounds == pytest.approx((-9.5, 39.5, -6.5, 36.5, 99.5, 133.5))


def test_foreground_voxel_center_bounds_from_mask_uses_xyz_array_order():
    mask = np.zeros((8, 9, 10), dtype=np.uint8)
    mask[2:6, 3:8, 1:7] = 1

    bounds = foreground_voxel_center_bounds_from_mask(
        mask,
        origin=(10, 20, 30),
        spacing=(2, 3, 4),
    )

    assert bounds == pytest.approx((14, 20, 29, 41, 34, 54))


def test_pad_vtk_images_to_physical_bounds_preserves_world_coordinates():
    vtk = pytest.importorskip("vtk")
    from vtk.util.numpy_support import numpy_to_vtk

    from ogo.util.vtk_image import vtk_image_to_numpy

    values = np.zeros((3, 3, 3), dtype=np.uint8)
    values[0, 0, 0] = 7
    image = vtk.vtkImageData()
    image.SetDimensions(values.shape)
    image.SetOrigin(10.0, 20.0, 30.0)
    image.SetSpacing(2.0, 1.0, 5.0)
    image.GetPointData().SetScalars(numpy_to_vtk(values.ravel(order="F"), deep=True))

    padded, padding = pad_vtk_images_to_physical_bounds(
        [image],
        desired_bounds=(6.0, 14.0, 20.0, 24.0, 30.0, 40.0),
        constants=[3],
    )

    out = padded[0]
    data = vtk_image_to_numpy(out)

    assert out.GetDimensions() == (5, 5, 3)
    assert out.GetOrigin() == pytest.approx((6.0, 20.0, 30.0))
    assert padding["lower"] == (2, 0, 0)
    assert padding["upper"] == (0, 2, 0)
    assert data[2, 0, 0] == 7
    assert data[0, 4, 0] == 3


def test_projected_material_disk_required_bounds_includes_lateral_plane_extent():
    vtk = pytest.importorskip("vtk")
    from vtk.util.numpy_support import numpy_to_vtk

    mask = np.zeros((5, 5, 5), dtype=np.uint8)
    mask[2, 2, 2] = 1
    image = vtk.vtkImageData()
    image.SetDimensions(mask.shape)
    image.SetOrigin(0.0, 0.0, 0.0)
    image.SetSpacing(1.0, 1.0, 1.0)
    image.GetPointData().SetScalars(numpy_to_vtk(mask.ravel(order="F"), deep=True))

    bounds = projected_material_disk_required_bounds(
        image,
        center=(2.0, 2.0, 6.0),
        normal=(0.0, 0.0, -1.0),
        u_axis=(1.0, 0.0, 0.0),
        v_axis=(0.0, 1.0, 0.0),
        size=(8.0, 6.0),
        thickness=4.0,
        intrusion=2.0,
    )

    assert bounds == pytest.approx((-3.0, 7.0, -2.0, 6.0, -1.0, 5.0))


def test_smooth_label_mask_vtk_preserves_likely_label_regions():
    vtk = pytest.importorskip("vtk")

    labels = np.zeros((9, 9, 9), dtype=np.uint8)
    labels[2:7, 2:7, 2:7] = 1
    labels[4:7, 4:7, 4:7] = 2

    from vtk.util.numpy_support import numpy_to_vtk

    from ogo.util.vtk_image import vtk_image_to_numpy

    image = vtk.vtkImageData()
    image.SetDimensions(labels.shape)
    image.SetSpacing(1.0, 1.0, 1.0)
    image.SetOrigin(0.0, 0.0, 0.0)
    image.GetPointData().SetScalars(
        numpy_to_vtk(labels.ravel(order="F"), deep=True)
    )

    smoothed = smooth_label_mask_vtk(image, sigma_mm=1.0, threshold=0.5)
    values = vtk_image_to_numpy(smoothed)

    assert set(np.unique(values)) <= {0, 1, 2}
    assert values[3, 3, 3] == 1
    assert values[5, 5, 5] == 2


def test_fit_cap_follows_each_contact_ray():
    mask = np.zeros((12, 8, 8), dtype=bool)
    mask[4:7, 2:6, 2:6] = True
    mask[5, 1, 3] = True

    cap = generate_bone_cap_mask(mask, axis="x", direction="down", thickness=3, shape="fit")

    assert cap.any()
    assert not np.any(cap & mask)
    for y in range(mask.shape[1]):
        for z in range(mask.shape[2]):
            body_x = np.flatnonzero(mask[:, y, z])
            cap_x = np.flatnonzero(cap[:, y, z])
            if cap_x.size == 0:
                continue
            assert body_x.size > 0
            assert cap_x.max() < body_x.min()


def test_box_cap_uses_projected_contact_extent():
    mask = np.zeros((12, 8, 8), dtype=bool)
    mask[4:7, 2:6, 2:6] = True

    cap = generate_bone_cap_mask(mask, axis="x", direction="up", thickness=2, shape="box")

    assert cap.sum() == 2 * 4 * 4
    assert np.all(cap[7:9, 2:6, 2:6])
    assert not np.any(cap & mask)


@pytest.mark.parametrize("shape", ["box", "round"])
def test_shaped_caps_touch_irregular_bone_surface(shape):
    mask = np.zeros((14, 9, 9), dtype=bool)
    mask[3:6, 2:7, 2:7] = True
    mask[6, 4, 4] = True
    mask[7, 5, 4] = True

    cap = generate_bone_cap_mask(mask, axis="x", direction="up", thickness=3, shape=shape)

    assert cap.any()
    assert not np.any(cap & mask)
    for y in range(mask.shape[1]):
        for z in range(mask.shape[2]):
            cap_x = np.flatnonzero(cap[:, y, z])
            if cap_x.size == 0:
                continue
            body_x = np.flatnonzero(mask[:, y, z])
            assert body_x.size > 0
            assert cap_x.min() == body_x.max() + 1


def test_round_cap_stays_inside_projected_box():
    mask = np.zeros((12, 8, 8), dtype=bool)
    mask[4:7, 2:6, 2:6] = True

    round_cap = generate_bone_cap_mask(mask, axis="x", direction="up", thickness=2, shape="round")
    box_cap = generate_bone_cap_mask(mask, axis="x", direction="up", thickness=2, shape="box")

    assert round_cap.any()
    assert round_cap.sum() < box_cap.sum()
    assert np.all(round_cap <= box_cap)
    assert not np.any(round_cap & mask)


def test_bone_cap_intrusion_limits_global_surface_depth():
    mask = np.zeros((14, 9, 9), dtype=bool)
    mask[4:6, 2:7, 2:7] = True
    mask[7, 4:6, 4:6] = True

    cap = generate_bone_cap_mask(
        mask,
        axis="x",
        direction="up",
        thickness=2,
        shape="fit",
        intrusion=1,
    )

    assert cap.any()
    assert not np.any(cap & mask)
    cap_x = np.where(cap)[0]
    assert cap_x.min() >= 7
    assert cap_x.max() <= 9
    assert cap[:, 2, 2].sum() == 0
    assert cap[:, 4, 4].sum() > 0


def test_bone_cap_zero_intrusion_stays_outside_global_surface():
    mask = np.zeros((12, 8, 8), dtype=bool)
    mask[4:7, 2:6, 2:6] = True

    cap = generate_bone_cap_mask(
        mask,
        axis="x",
        direction="up",
        thickness=2,
        shape="box",
        intrusion=0,
    )

    assert np.all(cap[7:9, 2:6, 2:6])
    assert not np.any(cap[:7])
    assert not np.any(cap & mask)


def test_projected_bone_cap_closes_two_voxel_footprint_gaps():
    mask = np.zeros((12, 12, 12), dtype=bool)
    mask[5:7, 2:5, 2:10] = True
    mask[5:7, 7:10, 2:10] = True

    cap = generate_bone_cap_mask(
        mask,
        axis="x",
        direction="up",
        thickness=4,
        shape="fit",
        intrusion=2,
    )

    assert cap.any()
    assert not np.any(cap & mask)
    assert np.all(cap[7:9, 5:7, 2:10])
    cap_x = np.where(cap)[0]
    assert cap_x.max() - cap_x.min() + 1 <= 4


def test_projected_bone_cap_intrusion_keeps_requested_total_thickness():
    mask = np.zeros((24, 12, 12), dtype=bool)
    mask[8:10, 2:10, 2:10] = True
    mask[14:16, 5:7, 5:7] = True

    cap = generate_bone_cap_mask(
        mask,
        axis="x",
        direction="up",
        thickness=6,
        shape="fit",
        intrusion=5,
    )

    assert cap.any()
    assert not np.any(cap & mask)
    cap_x = np.where(cap)[0]
    assert cap_x.max() - cap_x.min() + 1 <= 6


def test_anatomy_bone_cap_uses_axis_aligned_surface_depth():
    mask = np.zeros((16, 9, 9), dtype=bool)
    mask[7:9, 2:6, 2:6] = True
    mask[5:7, 6, 3] = True
    mask[2:4, 7, 3] = True

    cap = generate_bone_cap_mask(
        mask,
        axis="x",
        direction="up",
        thickness=5,
        shape="anatomy",
        intrusion=3,
    )

    assert not np.any(cap & mask)
    assert np.flatnonzero(cap[:, 3, 3]).tolist() == [9, 10]
    assert np.flatnonzero(cap[:, 6, 3]).tolist() == [7, 8, 9, 10]
    assert np.flatnonzero(cap[:, 7, 3]).tolist() == []
    for y in range(mask.shape[1]):
        for z in range(mask.shape[2]):
            assert np.flatnonzero(cap[:, y, z]).size <= 5


def test_anatomy_bone_cap_uses_axis_aligned_surface_depth_inferior():
    mask = np.zeros((18, 9, 9), dtype=bool)
    mask[7:9, 2:6, 2:6] = True
    mask[9:11, 6, 3] = True
    mask[12:14, 7, 3] = True

    cap = generate_bone_cap_mask(
        mask,
        axis="x",
        direction="down",
        thickness=5,
        shape="anatomy",
        intrusion=3,
    )

    assert not np.any(cap & mask)
    assert np.flatnonzero(cap[:, 3, 3]).tolist() == [5, 6]
    assert np.flatnonzero(cap[:, 6, 3]).tolist() == [5, 6, 7, 8]
    assert np.flatnonzero(cap[:, 7, 3]).tolist() == []
    for y in range(mask.shape[1]):
        for z in range(mask.shape[2]):
            assert np.flatnonzero(cap[:, y, z]).size <= 5


def test_vtk_bone_cap_uses_requested_output_value():
    vtk = pytest.importorskip("vtk")
    from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

    from ogo.fea.boundary import generate_bone_cap_vtk

    labels = np.zeros((12, 8, 8), dtype=np.uint16)
    labels[4:7, 2:6, 2:6] = 12
    vtk_arr = numpy_to_vtk(np.swapaxes(labels, 0, 2).ravel(order="F"), deep=True)
    image = vtk.vtkImageData()
    image.SetDimensions(list(np.swapaxes(labels, 0, 2).shape))
    image.SetSpacing((1, 1, 1))
    image.SetOrigin((0, 0, 0))
    image.GetPointData().SetScalars(vtk_arr)

    cap = generate_bone_cap_vtk(
        image, label_value=12, axis="x", direction="up", thickness=2, shape="box", output_value=5000
    )
    cap_np = vtk_to_numpy(cap.GetPointData().GetScalars()).reshape(cap.GetDimensions(), order="F")

    assert set(np.unique(cap_np)) == {0, 5000}


def test_projected_material_disk_uses_fixed_thickness_and_anatomy_intrusion():
    active = np.zeros((14, 14, 14), dtype=bool)
    active[4:10, 6:9, 4:10] = True
    active[7, 6:8, 7] = False
    active[7, 8:10, 7] = True
    active[8, 6:10, 7] = False
    active[8, 10:12, 7] = True

    disk = generate_projected_material_disk_mask(
        active,
        spacing=(1, 1, 1),
        origin=(0, 0, 0),
        center=(7, 1, 7),
        normal=(0, 1, 0),
        u_axis=(1, 0, 0),
        v_axis=(0, 0, 1),
        size=(8, 8),
        shape="square",
        thickness=6,
        intrusion=3,
        anatomy_constrained=True,
    )

    assert np.flatnonzero(disk[5, :, 5]).tolist() == [3, 4, 5]
    assert np.flatnonzero(disk[7, :, 7]).tolist() == [3, 4, 5, 6, 7]
    assert not np.any(disk[8, :, 7])
    assert not np.any(disk & active)


def test_projected_material_disk_vtk_preserves_geometry_and_output_value():
    vtk = pytest.importorskip("vtk")
    from ogo.fea.boundary import generate_projected_material_disk_vtk
    from ogo.util.vtk_image import numpy_to_vtk_image, vtk_image_to_numpy

    material = np.zeros((10, 10, 10), dtype=np.uint16)
    material[3:7, 5:7, 3:7] = 12

    template = vtk.vtkImageData()
    template.SetDimensions(material.shape)
    template.SetOrigin(-2, 4, 8)
    template.SetSpacing(1, 1, 1)
    image = numpy_to_vtk_image(material, template, vtk_array_type=vtk.VTK_UNSIGNED_SHORT)

    disk = generate_projected_material_disk_vtk(
        image,
        center=(2, 12, 12),
        normal=(0, -1, 0),
        u_axis=(1, 0, 0),
        v_axis=(0, 0, 1),
        size=(8, 8),
        shape="square",
        thickness=4,
        intrusion=2,
        anatomy_constrained=True,
        output_value=5000,
    )

    disk_data = vtk_image_to_numpy(disk)
    assert disk.GetExtent() == image.GetExtent()
    assert disk.GetSpacing() == image.GetSpacing()
    assert disk.GetOrigin() == image.GetOrigin()
    assert set(np.unique(disk_data)) == {0, 5000}
    assert not np.any((disk_data != 0) & (material != 0))


def test_projected_material_disk_vtk_keeps_pmma_out_of_protected_mask():
    vtk = pytest.importorskip("vtk")
    from ogo.fea.boundary import generate_projected_material_disk_vtk
    from ogo.util.vtk_image import numpy_to_vtk_image, vtk_image_to_numpy

    material = np.zeros((10, 10, 10), dtype=np.uint16)
    material[3:7, 5:7, 3:7] = 12

    femur_mask = np.zeros_like(material)
    femur_mask[3:7, 4:7, 3:7] = 1

    template = vtk.vtkImageData()
    template.SetDimensions(material.shape)
    template.SetOrigin(-2, 4, 8)
    template.SetSpacing(1, 1, 1)
    material_image = numpy_to_vtk_image(material, template, vtk_array_type=vtk.VTK_UNSIGNED_SHORT)
    mask_image = numpy_to_vtk_image(femur_mask, template, vtk_array_type=vtk.VTK_UNSIGNED_CHAR)

    disk = generate_projected_material_disk_vtk(
        material_image,
        surface_vtk_image=mask_image,
        exclusion_vtk_image=mask_image,
        center=(2, 12, 12),
        normal=(0, -1, 0),
        u_axis=(1, 0, 0),
        v_axis=(0, 0, 1),
        size=(8, 8),
        shape="square",
        thickness=4,
        intrusion=2,
        anatomy_constrained=True,
        output_value=5000,
    )

    disk_data = vtk_image_to_numpy(disk)
    assert np.any(disk_data != 0)
    assert not np.any((disk_data != 0) & (femur_mask != 0))
