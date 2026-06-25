import pytest

from ogo.fea import femur


def test_side_suffix_and_rotation_match_compact_outputs():
    assert femur.side_suffix(1) == "LF"
    assert femur.side_suffix(2) == "RF"
    assert femur.side_rotation(1) == 90
    assert femur.side_rotation(2) == -90


def test_sideways_fall_output_name_uses_compact_side_suffix():
    assert femur.sideways_fall_output_name("density.n88model", 1) == "density_LF.n88model"
    assert femur.sideways_fall_output_name("density.n88model", 2) == "density_RF.n88model"


def test_invalid_femur_side_is_rejected():
    with pytest.raises(ValueError):
        femur.side_suffix(3)


def test_tilted_side_support_vector_has_unit_length():
    np = pytest.importorskip("numpy")

    vector = femur.tilted_side_support_vector(-20)

    assert np.linalg.norm(vector) == pytest.approx(1.0)


def test_distal_cut_z_from_reference_uses_superior_bound():
    assert femur.distal_cut_z_from_reference((0, 1, 0, 1, -20, 120), 100) == pytest.approx(20)


def test_swap_xz_footprint_preserves_center_and_swaps_lengths():
    swapped = femur.swap_xz_footprint((10, 50, 1, 3, 100, 120))

    assert swapped == pytest.approx((20, 40, 1, 3, 90, 130))


def test_expand_z_footprint_preserves_center():
    expanded = femur.expand_z_footprint((20, 40, 1, 3, 90, 130), 20)

    assert expanded == pytest.approx((20, 40, 1, 3, 80, 140))


def test_expand_xz_footprint_preserves_center():
    expanded = femur.expand_xz_footprint(
        (20, 40, 1, 3, 90, 130),
        x_extension_mm=10,
        z_extension_mm=20,
    )

    assert expanded == pytest.approx((15, 45, 1, 3, 80, 140))


def test_flat_femur_shaft_crop_rejects_short_distal_coverage():
    np = pytest.importorskip("numpy")
    vtk = pytest.importorskip("vtk")
    from vtk.util.numpy_support import numpy_to_vtk

    data = np.zeros((4, 4, 4), dtype=np.uint8)
    data[:, :, 2:] = 1
    image = vtk.vtkImageData()
    image.SetDimensions(4, 4, 4)
    image.SetOrigin(0, 0, 0)
    image.SetSpacing(1, 1, 1)
    image.GetPointData().SetScalars(numpy_to_vtk(data.ravel(order="F"), deep=True))

    with pytest.raises(ValueError, match="does not include enough distal shaft"):
        femur.standardize_femur_shaft_length(
            image,
            image,
            reference_bounds=(0, 0, 0, 0, 0, 10),
            retained_length_mm=9,
            cut_mode="fixed_length",
        )


def test_flat_femur_shaft_crop_zeroes_voxels_below_cut():
    np = pytest.importorskip("numpy")
    vtk = pytest.importorskip("vtk")
    from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

    data = np.ones((3, 3, 5), dtype=np.uint8)
    image = vtk.vtkImageData()
    image.SetDimensions(3, 3, 5)
    image.SetOrigin(0, 0, 0)
    image.SetSpacing(1, 1, 1)
    image.GetPointData().SetScalars(numpy_to_vtk(data.ravel(order="F"), deep=True))

    _, cropped, meta = femur.standardize_femur_shaft_length(
        image,
        image,
        reference_bounds=(0, 0, 0, 0, 0, 4),
        retained_length_mm=2,
        cut_mode="fixed_length",
    )
    out = vtk_to_numpy(cropped.GetPointData().GetScalars()).reshape((3, 3, 5), order="F")

    assert meta["cut_z"] == pytest.approx(2)
    assert not any("femoral-head-side region" in warning for warning in meta["warnings"])
    assert not np.any(out[:, :, :2])
    assert np.all(out[:, :, 2:])


def test_femur_crop_warns_when_cut_enters_head_side_region():
    np = pytest.importorskip("numpy")
    vtk = pytest.importorskip("vtk")
    from vtk.util.numpy_support import numpy_to_vtk

    data = np.zeros((40, 40, 100), dtype=np.uint8)
    data[:, 0:8, 20:70] = 1
    data[:, 8:25, 10:95] = 1
    image = vtk.vtkImageData()
    image.SetDimensions(data.shape)
    image.SetOrigin(0, 0, 0)
    image.SetSpacing(1, 1, 1)
    image.GetPointData().SetScalars(numpy_to_vtk(data.ravel(order="F"), deep=True))

    _, _, meta = femur.standardize_femur_shaft_length(
        image,
        image,
        reference_bounds=(0, 0, 0, 0, 0, 90),
        retained_length_mm=30,
        cut_mode="fixed_length",
    )

    assert any("femoral-head-side region" in warning for warning in meta["warnings"])
    assert any("unusually short" in warning for warning in meta["warnings"])


def test_femur_crop_does_not_warn_short_length_above_threshold():
    np = pytest.importorskip("numpy")
    vtk = pytest.importorskip("vtk")
    from vtk.util.numpy_support import numpy_to_vtk

    data = np.zeros((40, 40, 100), dtype=np.uint8)
    data[:, 10:30, 20:95] = 1
    image = vtk.vtkImageData()
    image.SetDimensions(data.shape)
    image.SetOrigin(0, 0, 0)
    image.SetSpacing(1, 1, 1)
    image.GetPointData().SetScalars(numpy_to_vtk(data.ravel(order="F"), deep=True))

    _, _, meta = femur.standardize_femur_shaft_length(
        image,
        image,
        reference_bounds=(0, 0, 0, 0, 0, 90),
        retained_length_mm=56,
        cut_mode="fixed_length",
    )

    assert not any("unusually short" in warning for warning in meta["warnings"])


def test_femur_input_padding_adds_only_missing_foreground_margin():
    np = pytest.importorskip("numpy")
    vtk = pytest.importorskip("vtk")
    from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

    data = np.zeros((5, 6, 7), dtype=np.uint8)
    data[0:3, 2:5, 0:4] = 1
    image = vtk.vtkImageData()
    image.SetDimensions(data.shape)
    image.SetOrigin(0, 0, 0)
    image.SetSpacing(1, 2, 1)
    image.GetPointData().SetScalars(numpy_to_vtk(data.ravel(order="F"), deep=True))

    padded, meta = femur.pad_vtk_images_to_foreground_margin([image], image, margin_mm=2)
    out = vtk_to_numpy(padded[0].GetPointData().GetScalars()).reshape(
        padded[0].GetDimensions(),
        order="F",
    )
    coords = np.array(np.where(out != 0))

    assert meta["lower"] == (2, 0, 2)
    assert meta["upper"] == (0, 0, 0)
    assert tuple(coords.min(axis=1)) == (2, 2, 2)
    assert out.shape[0] - 1 - coords.max(axis=1)[0] >= 2
    assert out.shape[1] - 1 - coords.max(axis=1)[1] >= 1
    assert out.shape[2] - 1 - coords.max(axis=1)[2] >= 2


def test_lesser_trochanter_cut_uses_distal_area_peak():
    np = pytest.importorskip("numpy")
    vtk = pytest.importorskip("vtk")
    from vtk.util.numpy_support import numpy_to_vtk

    data = np.zeros((60, 80, 90), dtype=np.uint8)
    x = np.arange(data.shape[0])[:, None]
    y = np.arange(data.shape[1])[None, :]
    for z in range(10, 86):
        radius = 10
        y_center = 35
        if 58 <= z <= 64:
            radius = 19
        if 73 <= z <= 79:
            y_center = 50
            radius = 12
        section = ((x - 30) ** 2 + (y - y_center) ** 2) <= radius**2
        data[:, :, z] = section.astype(np.uint8)

    image = vtk.vtkImageData()
    image.SetDimensions(data.shape)
    image.SetOrigin(0, 0, 0)
    image.SetSpacing(1, 1, 1)
    image.GetPointData().SetScalars(numpy_to_vtk(data.ravel(order="F"), deep=True))

    meta = femur.detect_lesser_trochanter_cut_z(image)

    assert meta["greater_trochanter_z"] == pytest.approx(76)
    assert meta["lesser_trochanter_z"] == pytest.approx(61)
    assert meta["cut_z"] == pytest.approx(61)


def test_lesser_trochanter_cut_uses_peak_plateau_center_and_percent_offset():
    np = pytest.importorskip("numpy")
    vtk = pytest.importorskip("vtk")
    from vtk.util.numpy_support import numpy_to_vtk

    data = np.zeros((80, 90, 100), dtype=np.uint8)
    x = np.arange(data.shape[0])[:, None]
    y = np.arange(data.shape[1])[None, :]
    for z in range(10, 96):
        radius = 10
        y_center = 35
        if 54 <= z <= 62:
            radius = 20
        if 78 <= z <= 82:
            y_center = 58
            radius = 12
        section = ((x - 40) ** 2 + (y - y_center) ** 2) <= radius**2
        data[:, :, z] = section.astype(np.uint8)

    image = vtk.vtkImageData()
    image.SetDimensions(data.shape)
    image.SetOrigin(0, 0, 0)
    image.SetSpacing(1, 1, 1)
    image.GetPointData().SetScalars(numpy_to_vtk(data.ravel(order="F"), deep=True))

    meta = femur.detect_lesser_trochanter_cut_z(
        image,
        distal_offset_percent=50.0,
        max_distal_to_greater_mm=45.0,
    )

    assert meta["greater_trochanter_z"] == pytest.approx(80)
    assert meta["lesser_trochanter_z"] == pytest.approx(58)
    assert meta["distal_offset_mm"] == pytest.approx(11)
    assert meta["cut_z"] == pytest.approx(47)


def test_cortical_compartment_mask_uses_default_labels():
    np = pytest.importorskip("numpy")
    vtk = pytest.importorskip("vtk")
    from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

    data = np.zeros((3, 3, 3), dtype=np.uint8)
    data[0, :, :] = 1
    data[1:, :, :] = 2
    image = vtk.vtkImageData()
    image.SetDimensions(data.shape)
    image.GetPointData().SetScalars(numpy_to_vtk(data.ravel(order="F"), deep=True))

    cortical = femur.cortical_compartment_mask(image)
    out = vtk_to_numpy(cortical.GetPointData().GetScalars()).reshape(data.shape, order="F")

    assert np.all(out[0, :, :] == 1)
    assert not np.any(out[1:, :, :])


def test_cortical_compartment_mask_requires_trab_and_cort_labels():
    np = pytest.importorskip("numpy")
    vtk = pytest.importorskip("vtk")
    from vtk.util.numpy_support import numpy_to_vtk

    data = np.ones((3, 3, 3), dtype=np.uint8)
    image = vtk.vtkImageData()
    image.SetDimensions(data.shape)
    image.GetPointData().SetScalars(numpy_to_vtk(data.ravel(order="F"), deep=True))

    with pytest.raises(ValueError, match="missing required label"):
        femur.cortical_compartment_mask(image)
