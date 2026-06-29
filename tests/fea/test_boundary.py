import pytest

np = pytest.importorskip("numpy")

from ogo.fea.boundary import (
    generate_bone_cap_mask,
    generate_fixture_cap_mask,
)


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
    mask[2:4, 6, 3] = True

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
    assert np.flatnonzero(cap[:, 6, 3]).tolist() == [4, 5, 6, 7, 8, 9, 10]


def test_anatomy_bone_cap_uses_axis_aligned_surface_depth_inferior():
    mask = np.zeros((18, 9, 9), dtype=bool)
    mask[7:9, 2:6, 2:6] = True
    mask[12:14, 6, 3] = True

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
    assert np.flatnonzero(cap[:, 6, 3]).tolist() == [5, 6, 7, 8, 9, 10, 11]


def test_fixture_box_cap_keeps_full_footprint_over_irregular_surface():
    mask = np.zeros((14, 9, 9), dtype=bool)
    mask[3:6, 2:7, 2:7] = True
    mask[6, 4, 4] = True
    mask[7, 5, 4] = True

    cap = generate_fixture_cap_mask(mask, axis="x", direction="up", thickness=3, shape="box")

    assert np.all(cap[8:11, 2:7, 2:7])
    assert cap.sum() == 3 * 5 * 5


def test_fixture_cap_intrudes_across_contact_plane():
    mask = np.zeros((14, 9, 9), dtype=bool)
    mask[3:6, 2:7, 2:7] = True

    cap = generate_fixture_cap_mask(
        mask,
        axis="x",
        direction="up",
        thickness=3,
        shape="box",
        intrusion=2,
    )

    assert np.all(cap[4:7, 2:7, 2:7])
    assert not np.any(cap[:4])
    assert not np.any(cap[7:])


def test_fixture_cap_intrusion_keeps_total_thickness():
    mask = np.zeros((18, 9, 9), dtype=bool)
    mask[8:10, 2:7, 2:7] = True

    half_intrusion = generate_fixture_cap_mask(
        mask,
        axis="x",
        direction="up",
        thickness=6,
        shape="box",
        intrusion=3,
    )
    full_intrusion = generate_fixture_cap_mask(
        mask,
        axis="x",
        direction="up",
        thickness=6,
        shape="box",
        intrusion=6,
    )

    half_x = np.flatnonzero(half_intrusion[:, 3, 3])
    full_x = np.flatnonzero(full_intrusion[:, 3, 3])
    assert half_x.tolist() == [7, 8, 9, 10, 11, 12]
    assert full_x.tolist() == [4, 5, 6, 7, 8, 9]
    assert half_x.size == 6
    assert full_x.size == 6


def test_fixture_contact_crop_removes_unsupported_columns():
    mask = np.zeros((14, 9, 9), dtype=bool)
    mask[3:6, 2:7, 4] = True
    mask[3:6, 4, 2:7] = True
    mask[0, 3, 3] = True

    full = generate_fixture_cap_mask(mask, axis="x", direction="up", thickness=2, shape="box")
    cropped = generate_fixture_cap_mask(
        mask,
        axis="x",
        direction="up",
        thickness=2,
        shape="box",
        intrusion=2,
        crop_to_contact=True,
    )

    assert full[:, 3, 3].any()
    assert not cropped[:, 3, 3].any()
    assert cropped[:, 4, 3].any()


def test_fixture_round_cap_is_not_contact_cropped():
    mask = np.zeros((14, 11, 11), dtype=bool)
    mask[3:6, 2:9, 5] = True
    mask[3:6, 5, 2:9] = True

    fixture = generate_fixture_cap_mask(mask, axis="x", direction="up", thickness=2, shape="round")
    contact = generate_bone_cap_mask(mask, axis="x", direction="up", thickness=2, shape="round")

    assert fixture[:, 3, 3].any()
    assert not contact[:, 3, 3].any()
    assert not np.any(contact & mask)


def test_vtk_roi_cap_uses_requested_output_value():
    vtk = pytest.importorskip("vtk")
    from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

    from ogo.fea.boundary import generate_bone_cap_vtk, label_foreground_in_bounds_vtk

    labels = np.zeros((12, 8, 8), dtype=np.uint16)
    labels[4:7, 2:6, 2:6] = 12
    vtk_arr = numpy_to_vtk(np.swapaxes(labels, 0, 2).ravel(order="F"), deep=True)
    image = vtk.vtkImageData()
    image.SetDimensions(list(np.swapaxes(labels, 0, 2).shape))
    image.SetSpacing((1, 1, 1))
    image.SetOrigin((0, 0, 0))
    image.GetPointData().SetScalars(vtk_arr)

    roi = label_foreground_in_bounds_vtk(image, (0, 11, 2, 5, 2, 5), label_value=1)
    cap = generate_bone_cap_vtk(
        roi, label_value=1, axis="x", direction="up", thickness=2, shape="box", output_value=5000
    )
    cap_np = vtk_to_numpy(cap.GetPointData().GetScalars()).reshape(cap.GetDimensions(), order="F")

    assert set(np.unique(cap_np)) == {0, 5000}


def test_vtk_conversion_preserves_template_extent():
    vtk = pytest.importorskip("vtk")
    from vtk.util.numpy_support import numpy_to_vtk

    from ogo.fea.boundary import generate_fixture_cap_vtk

    labels = np.zeros((6, 5, 4), dtype=np.uint16)
    labels[2:4, 1:4, 1:3] = 1
    vtk_arr = numpy_to_vtk(np.swapaxes(labels, 0, 2).ravel(order="F"), deep=True)
    image = vtk.vtkImageData()
    image.SetExtent(10, 13, 20, 24, 30, 35)
    image.SetSpacing((0.7, 0.8, 0.9))
    image.SetOrigin((1.0, 2.0, 3.0))
    image.GetPointData().SetScalars(vtk_arr)

    cap = generate_fixture_cap_vtk(
        image,
        label_value=1,
        axis="x",
        direction="up",
        thickness=1,
        shape="box",
        output_value=5000,
    )

    assert cap.GetExtent() == image.GetExtent()
    assert cap.GetSpacing() == image.GetSpacing()
    assert cap.GetOrigin() == image.GetOrigin()
