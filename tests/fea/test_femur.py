import pytest

from ogo.fea import femur
from ogo.fea.image_io import write_vtk_image_with_sitk_geometry


@pytest.mark.parametrize("options", [
    ["--femur_cut_mode", "post_icp_oblique_ratio"],
    ["--femur_bbox_ratio", "1", "1.2", "none"],
    ["--femur_lesser_trochanter_distal_offset", "20"],
])
def test_builder_rejects_obsolete_crop_methods(monkeypatch, options):
    monkeypatch.setattr("sys.argv", ["hip", "density.nii.gz", "mask.nii.gz"] + options)
    with pytest.raises(SystemExit) as error:
        femur.main()
    assert error.value.code == 2


def test_builder_has_one_gt_disk_relative_recipe(monkeypatch):
    captured = []
    monkeypatch.setattr("sys.argv", ["hip", "density.nii.gz", "mask.nii.gz",
                                    "--femur_greater_trochanter_distal_length", "15"])
    monkeypatch.setattr(femur, "sidewaysFallFe", captured.append)
    femur.main()
    args = captured[0]
    assert args.femur_cut_mode == "greater_trochanter_length"
    assert args.femur_greater_trochanter_distal_length == 15
    assert args.femur_shaft_length == 120
    assert not hasattr(args, "femur_bbox_ratio")


def test_builder_defaults_match_current_short_shaft_protocol(monkeypatch):
    captured = []
    monkeypatch.setattr("sys.argv", ["hip", "density.nii.gz", "mask.nii.gz"])
    monkeypatch.setattr(femur, "sidewaysFallFe", captured.append)
    femur.main()
    assert captured[0].femur_greater_trochanter_distal_length == 10
    assert femur.DEFAULT_FEMUR_PISTOIA_CRITICAL_VOLUME_PERCENT == 11.2


def _vtk_image_from_array(data, *, origin=(0, 0, 0), spacing=(1, 1, 1)):
    vtk = pytest.importorskip("vtk")
    from vtk.util.numpy_support import numpy_to_vtk

    image = vtk.vtkImageData()
    image.SetDimensions(data.shape)
    image.SetOrigin(*origin)
    image.SetSpacing(*spacing)
    image.GetPointData().SetScalars(numpy_to_vtk(data.ravel(order="F"), deep=True))
    return image


def _polydata_from_points(points):
    vtk = pytest.importorskip("vtk")

    vtk_points = vtk.vtkPoints()
    vertices = vtk.vtkCellArray()
    for point in points:
        point_id = vtk_points.InsertNextPoint(*point)
        vertices.InsertNextCell(1)
        vertices.InsertCellPoint(point_id)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(vtk_points)
    polydata.SetVerts(vertices)
    return polydata


def test_preparation_and_fixed_registration_preserve_full_femur(tmp_path, monkeypatch):
    np = pytest.importorskip("numpy")
    from ogo.fea.boundary import vtk_image_to_numpy

    captured = []
    monkeypatch.setattr("sys.argv", ["hip", "density.nii.gz", "mask.nii.gz"])
    monkeypatch.setattr(femur, "sidewaysFallFe", captured.append)
    femur.main()
    args = captured[0]
    density = _vtk_image_from_array(np.full((12, 12, 150), 200, dtype=np.float32))
    labels = np.zeros((12, 12, 150), dtype=np.uint8)
    labels[3:9, 3:9, 5:145] = 1
    mask = _vtk_image_from_array(labels)
    monkeypatch.setattr(femur.ogo, "readNii", lambda path: density if path == args.calibrated_image else mask)
    transform = tmp_path / "fixed_icp.json"
    femur.write_icp_transform(
        transform, matrix=np.eye(4),
        icp_transform={"rotation": np.eye(3).tolist(), "iterations": 0, "mean_distance": 0},
        reference_scale={}, femur_side=1, rough_crop={},
    )
    args.femur_icp_transform_in = str(transform)

    prepared = femur._prepare_femur_inputs(args)
    full_mask = prepared[1]
    registration_mask = prepared[4]
    assert np.count_nonzero(vtk_image_to_numpy(full_mask)) == np.count_nonzero(labels)
    assert np.count_nonzero(vtk_image_to_numpy(registration_mask)) < np.count_nonzero(labels)
    registered = femur._register_femur_inputs(args, str(tmp_path / "model.n88model"), *prepared)
    assert np.count_nonzero(vtk_image_to_numpy(registered[1])) == np.count_nonzero(labels)
    assert registered[4] is False


def test_side_suffix_matches_compact_outputs():
    assert femur.side_suffix(1) == "LF"
    assert femur.side_suffix(2) == "RF"


def test_sideways_fall_output_name_uses_compact_side_suffix():
    assert femur.sideways_fall_output_name("density.n88model", 1) == "density_LF.n88model"
    assert femur.sideways_fall_output_name("density.n88model", 2) == "density_RF.n88model"


def test_pistoia_mask_output_path_uses_model_stem(tmp_path):
    assert (
        femur.pistoia_mask_output_path(tmp_path / "case_LF.n88model")
        == tmp_path / "case_LF_pistoia_mask.nii.gz"
    )


def test_write_vtk_image_with_sitk_geometry_preserves_origin_and_spacing(tmp_path):
    np = pytest.importorskip("numpy")
    sitk = pytest.importorskip("SimpleITK")

    image = _vtk_image_from_array(
        np.ones((3, 4, 5), dtype=np.uint8),
        origin=(-10.5, 20.25, 3.75),
        spacing=(0.8, 0.9, 1.1),
    )
    out = write_vtk_image_with_sitk_geometry(image, tmp_path / "mask.nii.gz")
    reread = sitk.ReadImage(str(out))

    assert reread.GetSize() == (3, 4, 5)
    assert reread.GetOrigin() == pytest.approx((-10.5, 20.25, 3.75))
    assert reread.GetSpacing() == pytest.approx((0.8, 0.9, 1.1))


def test_icp_transform_sidecar_round_trips_matrix(tmp_path):
    np = pytest.importorskip("numpy")

    matrix = np.eye(4)
    matrix[:3, 3] = [1.5, -2.0, 3.25]
    path = tmp_path / "sub-001_left_icp.json"

    femur.write_icp_transform(
        path,
        matrix=matrix,
        icp_transform={"iterations": 12, "mean_distance": 0.25},
        reference_scale={"source": "test"},
        femur_side=1,
        rough_crop={"retained_length_mm": 120.0},
    )

    data, loaded = femur.load_icp_transform(path)

    assert data["type"] == "ogo.femur.icp_transform"
    assert data["femur_side"] == 1
    assert data["icp"] == {"iterations": 12, "mean_distance": 0.25}
    assert femur.matrix4x4_to_numpy(loaded) == pytest.approx(matrix)


def test_invalid_femur_side_is_rejected():
    with pytest.raises(ValueError):
        femur.side_suffix(3)


def test_hip_fixture_defaults_use_fixed_thickness_and_intrusion():
    assert femur.DEFAULT_PMMA_THICKNESS_MM == pytest.approx(10.0)
    assert femur.DEFAULT_PMMA_INTRUSION_MM == pytest.approx(6.0)


def test_bbox_relative_fixture_bounds_scale_lateral_axes_to_model_bbox():
    bounds = (10, 50, -20, 80, 100, 200)

    fixture_bounds = femur.bbox_relative_fixture_bounds(
        bounds,
        center_fraction=(0.5, 1.1, 0.25),
        size_fraction=(1.1, 0.5),
        projection_axis="y",
    )

    assert fixture_bounds == pytest.approx((8, 52, -20, 80, 100, 150))


def test_bbox_relative_fixture_bounds_can_enforce_square_footprint():
    bounds = (10, 50, -20, 80, 100, 200)

    fixture_bounds = femur.bbox_relative_fixture_bounds(
        bounds,
        center_fraction=(0.5, 1.1, 0.25),
        size_fraction=(1.1, 0.5),
        projection_axis="y",
        shape="square",
    )

    assert fixture_bounds == pytest.approx((8, 52, -20, 80, 103, 147))


def test_foreground_voxel_center_bounds_match_workflow_bbox_convention():
    np = pytest.importorskip("numpy")

    mask = np.zeros((8, 9, 10), dtype=np.uint8)
    mask[2:6, 3:8, 1:7] = 1

    bounds = femur.foreground_voxel_center_bounds_from_mask(
        mask,
        origin=(10, 20, 30),
        spacing=(2, 3, 4),
    )

    assert bounds == pytest.approx((14, 20, 29, 41, 34, 54))


def test_reference_grid_from_output_to_input_matrix_uses_transformed_surface_bounds():
    np = pytest.importorskip("numpy")

    points = np.asarray(
        [
            [10.2, -4.0, 2.5],
            [14.7, 6.1, 8.2],
            [12.0, 0.5, 5.0],
        ],
        dtype=float,
    )
    matrix = np.eye(4)
    matrix[:3, 3] = [2.0, -3.0, 1.0]

    origin, size = femur.reference_grid_from_output_to_input_matrix(
        points,
        matrix,
        spacing=(1.0, 1.0, 1.0),
        margin_voxels=2,
    )

    # The matrix maps output coordinates to input coordinates, so the grid is
    # built from the inverse-transformed input points.
    expected = points - np.asarray([2.0, -3.0, 1.0])
    expected_lo = expected.min(axis=0) - 2.0
    expected_hi = expected.max(axis=0) + 2.0
    assert origin == pytest.approx(tuple(expected_lo))
    assert size == tuple((np.ceil(expected_hi - expected_lo).astype(int) + 1).tolist())


def test_bbox_relative_fixture_direction_uses_plane_side_fraction():
    assert femur.bbox_relative_fixture_direction((0.5, 0.02, 0.5), projection_axis="y") == "down"
    assert femur.bbox_relative_fixture_direction((0.5, 1.02, 0.5), projection_axis="y") == "up"


def test_bbox_relative_fixture_plane_uses_authored_y_plane_axes():
    bounds = (10, 50, -20, 80, 100, 200)

    low_y_plane = femur.bbox_relative_fixture_plane(
        bounds,
        center_fraction=(0.5, 0.02, 0.25),
        size_fraction=(1.1, 0.5),
        projection_axis="y",
        shape="rectangle",
    )
    high_y_plane = femur.bbox_relative_fixture_plane(
        bounds,
        center_fraction=(0.5, 1.02, 0.25),
        size_fraction=(1.1, 0.5),
        projection_axis="y",
        shape="rectangle",
    )

    assert low_y_plane["center"] == pytest.approx((30, -18, 125))
    assert low_y_plane["normal"] == pytest.approx((0, 1, 0))
    assert low_y_plane["u_axis"] == pytest.approx((0, 0, -1))
    assert low_y_plane["v_axis"] == pytest.approx((-1, 0, 0))
    assert low_y_plane["size"] == pytest.approx((110, 20))
    assert low_y_plane["shape"] == "rectangle"
    assert high_y_plane["normal"] == pytest.approx((0, -1, 0))
    assert high_y_plane["u_axis"] == pytest.approx((0, 0, 1))
    assert high_y_plane["v_axis"] == pytest.approx((-1, 0, 0))


def test_proximal_sideways_fall_fixture_plane_does_not_reach_distal_shaft():
    bounds = (-30, 30, -40, 80, -95, 120)

    plane = femur.proximal_sideways_fall_fixture_plane(
        bounds,
        center_fraction=femur.FEMORAL_HEAD_FIXTURE_CENTER_FRACTION,
    )

    assert plane["footprint"] == "proximal_only"
    assert plane["center"][2] == pytest.approx(80.0)
    assert plane["size"] == pytest.approx((80.0, 70.0))
    assert plane["center"][2] - plane["size"][0] / 2.0 == pytest.approx(40.0)
    assert plane["center"][2] + plane["size"][0] / 2.0 == pytest.approx(120.0)


def test_scale_reference_to_sample_principal_lengths_clips_and_scales_origin():
    reference = _polydata_from_points(
        [
            (-1.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, -2.0, 0.0),
            (0.0, 2.0, 0.0),
            (0.0, 0.0, -3.0),
            (0.0, 0.0, 3.0),
        ]
    )
    sample = _polydata_from_points(
        [
            (-0.5, 0.0, 0.0),
            (0.5, 0.0, 0.0),
            (0.0, -4.0, 0.0),
            (0.0, 4.0, 0.0),
            (0.0, 0.0, -6.0),
            (0.0, 0.0, 6.0),
        ]
    )

    scaled, metadata = femur.scale_reference_to_sample_principal_lengths(
        reference,
        sample,
        min_scale=(0.8, 0.8, 0.75),
        max_scale=(1.2, 1.2, 1.3),
    )

    assert metadata["scale_factors"] == pytest.approx([0.8, 1.2, 1.3])
    assert scaled.GetBounds() == pytest.approx((-0.8, 0.8, -2.4, 2.4, -3.9, 3.9))


def test_scale_reference_point_cloud_to_sample_preserves_reference_center():
    np = pytest.importorskip("numpy")

    reference = _polydata_from_points(
        [
            (9.0, 20.0, 30.0),
            (11.0, 20.0, 30.0),
            (10.0, 18.0, 30.0),
            (10.0, 22.0, 30.0),
            (10.0, 20.0, 27.0),
            (10.0, 20.0, 33.0),
        ]
    )
    sample_points = np.asarray(
        [
            (-1.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, -4.0, 0.0),
            (0.0, 4.0, 0.0),
            (0.0, 0.0, -6.0),
            (0.0, 0.0, 6.0),
        ]
    )

    scaled, metadata = femur.scale_reference_point_cloud_to_sample(
        reference,
        sample_points,
        min_scale=(0.8, 0.8, 0.75),
        max_scale=(1.2, 1.2, 1.3),
    )

    points = np.asarray(
        [scaled.GetPoint(point_id) for point_id in range(scaled.GetNumberOfPoints())]
    )
    assert points.mean(axis=0) == pytest.approx((10.0, 20.0, 30.0))
    assert metadata["reference_center"] == pytest.approx([10.0, 20.0, 30.0])
    assert metadata["scale_factors"] == pytest.approx([1.0, 1.2, 1.3])


def test_surface_points_from_vtk_mask_uses_voxel_surface_and_stride_sampling():
    np = pytest.importorskip("numpy")

    mask = np.ones((3, 3, 3), dtype=np.uint8)
    image = _vtk_image_from_array(mask, origin=(10, 20, 30), spacing=(2, 3, 4))

    full = femur.surface_points_from_vtk_mask(image, max_points=None)
    sampled = femur.surface_points_from_vtk_mask(
        image,
        max_points=5,
        sample_mode="stride",
    )

    assert full.shape == (26, 3)
    assert not np.any(np.all(np.isclose(full, (12, 23, 34)), axis=1))
    assert sampled == pytest.approx(full[[0, 6, 12, 18, 24]])


def test_scale_reference_point_cloud_to_sample_matches_voxel_surface_lengths():
    np = pytest.importorskip("numpy")

    reference = _polydata_from_points(
        [
            (-1.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, -2.0, 0.0),
            (0.0, 2.0, 0.0),
            (0.0, 0.0, -3.0),
            (0.0, 0.0, 3.0),
        ]
    )
    sample_points = np.asarray(
        [
            (-0.5, 0.0, 0.0),
            (0.5, 0.0, 0.0),
            (0.0, -4.0, 0.0),
            (0.0, 4.0, 0.0),
            (0.0, 0.0, -6.0),
            (0.0, 0.0, 6.0),
        ]
    )

    scaled, metadata = femur.scale_reference_point_cloud_to_sample(
        reference,
        sample_points,
        min_scale=(0.8, 0.8, 0.75),
        max_scale=(1.2, 1.2, 1.3),
    )

    assert metadata["source"] == "voxel_surface_point_cloud"
    assert metadata["scale_factors"] == pytest.approx([0.8, 1.2, 1.3])
    assert scaled.GetBounds() == pytest.approx((-0.8, 0.8, -2.4, 2.4, -3.9, 3.9))


def test_fixed_proximal_length_crop_keeps_proximal_side():
    np = pytest.importorskip("numpy")
    from ogo.util.vtk_image import vtk_image_to_numpy

    density = np.zeros((20, 20, 180), dtype=np.float32)
    mask = np.zeros_like(density, dtype=np.uint8)
    density[5:15, 5:15, 10:170] = 700.0
    mask[5:15, 5:15, 10:170] = 2

    cropped, crop_face, meta = femur.crop_vtk_images_to_fixed_proximal_length(
        [_vtk_image_from_array(density)],
        _vtk_image_from_array(mask),
        retained_length_mm=100.0,
        labels={2},
    )

    cropped_mask = vtk_image_to_numpy(cropped[0]) != 0
    coords = np.argwhere(cropped_mask)
    assert meta["method"] == "fixed_proximal_length"
    assert meta["retained_length_mm"] == pytest.approx(100.0)
    assert meta["status"] == "cropped"
    assert coords[:, 2].min() == 0
    assert coords[:, 2].max() == 99
    assert cropped[0].GetOrigin()[2] == pytest.approx(70.0)
    assert vtk_image_to_numpy(crop_face).sum() > 0


def test_straight_crop_face_support_uses_central_fraction_on_flat_face():
    np = pytest.importorskip("numpy")
    from ogo.util.vtk_image import vtk_image_to_numpy

    active = np.zeros((12, 12, 8), dtype=np.uint8)
    active[1:11, 1:11, 3:7] = 1
    crop_face = np.zeros_like(active)
    crop_face[1:11, 1:11, 3] = 1

    surface = femur.straight_crop_face_support_surface_vtk(
        _vtk_image_from_array(crop_face),
        _vtk_image_from_array(active),
        support_fraction=0.9,
    )
    surface_data = vtk_image_to_numpy(surface) != 0

    assert surface_data.sum() == 9 * 9
    assert np.all(surface_data[1:10, 1:10, 3])
    assert not np.any(surface_data[:, :, :3])
    assert not np.any(surface_data[:, :, 4:])


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
    assert padded[0].GetOrigin() == pytest.approx((-2, 0, -2))
    assert padded[0].GetExtent() == (0, 6, 0, 5, 0, 8)
    assert tuple(coords.min(axis=1)) == (2, 2, 2)
    assert out.shape[0] - 1 - coords.max(axis=1)[0] >= 2
    assert out.shape[1] - 1 - coords.max(axis=1)[1] >= 1
    assert out.shape[2] - 1 - coords.max(axis=1)[2] >= 2


def test_greater_trochanter_length_crop_uses_gt_support_distal_edge_origin():
    np = pytest.importorskip("numpy")
    vtk = pytest.importorskip("vtk")
    from ogo.util.vtk_image import vtk_image_to_numpy
    from vtk.util.numpy_support import numpy_to_vtk

    data = np.zeros((60, 80, 130), dtype=np.uint8)
    x = np.arange(data.shape[0])[:, None]
    y = np.arange(data.shape[1])[None, :]
    for z in range(10, 126):
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

    cropped, crop_face, meta = femur.crop_vtk_images_to_greater_trochanter_length(
        [image],
        image,
        gt_support_vtk=_vtk_image_from_array(data * ((np.arange(130) >= 74) & (np.arange(130) <= 79))[None, None, :]),
        retained_length_mm=10.0,
        distal_boundary={"required_flat_face_z_mm": 9.5},
        labels={1},
    )

    cropped_mask = vtk_image_to_numpy(cropped[0]) != 0
    assert meta["method"] == "greater_trochanter_length"
    assert meta["greater_trochanter_support_distal_edge_z"] == pytest.approx(73.5)
    assert meta["distal_length_origin_z"] == pytest.approx(73.5)
    assert meta["cut_z_mm"] == pytest.approx(63.5)
    assert meta["available_below_gt_support_distal_edge_mm"] == pytest.approx(64)
    assert meta["retained_length_mm"] == pytest.approx(10)
    assert cropped[0].GetOrigin()[2] == pytest.approx(64)
    assert np.argwhere(cropped_mask)[:, 2].min() == 0
    assert vtk_image_to_numpy(crop_face).sum() > 0


def test_greater_trochanter_length_crop_rejects_short_femur():
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

    with pytest.raises(ValueError, match="too short for the requested greater-trochanter"):
        femur.crop_vtk_images_to_greater_trochanter_length(
            [image],
            image,
            gt_support_vtk=_vtk_image_from_array(data * ((np.arange(90) >= 74) & (np.arange(90) <= 79))[None, None, :]),
            retained_length_mm=80.0,
            distal_boundary={"required_flat_face_z_mm": 9.5},
            labels={1},
        )


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
