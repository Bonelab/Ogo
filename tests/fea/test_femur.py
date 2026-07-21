import pytest

from ogo.fea import femur


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


def test_side_suffix_matches_compact_outputs():
    assert femur.side_suffix(1) == "LF"
    assert femur.side_suffix(2) == "RF"


def test_sideways_fall_output_name_uses_compact_side_suffix():
    assert femur.sideways_fall_output_name("density.n88model", 1) == "density_LF.n88model"
    assert femur.sideways_fall_output_name("density.n88model", 2) == "density_RF.n88model"


def test_invalid_femur_side_is_rejected():
    with pytest.raises(ValueError):
        femur.side_suffix(3)


def test_hip_fixture_defaults_use_fixed_thickness_and_intrusion():
    assert femur.DEFAULT_PMMA_THICKNESS_MM == pytest.approx(10.0)
    assert femur.DEFAULT_PMMA_INTRUSION_MM == pytest.approx(6.0)


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


def test_bbox_relative_oriented_contact_plane_keeps_measured_normal():
    np = pytest.importorskip("numpy")
    bounds = (10, 50, -20, 80, 100, 200)

    plane = femur.bbox_relative_oriented_contact_plane(
        bounds,
        center_fraction=(0.25, 0.5, 0.1),
        size_fraction=(0.4, 0.5),
        normal=(0.0, -0.5, 0.8660254038),
        shape="anatomy",
    )

    assert plane["center"] == pytest.approx((20.0, 30.0, 110.0))
    assert plane["normal"] == pytest.approx((0.0, -0.5, 0.8660254038))
    assert plane["outward_normal"] == pytest.approx((0.0, 0.5, -0.8660254038))
    assert plane["v_axis"] == pytest.approx((-1.0, 0.0, 0.0))
    assert np.cross(plane["u_axis"], plane["v_axis"]) == pytest.approx(plane["normal"])
    assert plane["size"] == pytest.approx((54.641016152, 20.0))
    assert plane["shape"] == "anatomy"


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


def test_bbox_ratio_crop_matches_recipe_order_and_tracks_cut_face():
    np = pytest.importorskip("numpy")
    from ogo.util.vtk_image import vtk_image_to_numpy

    density = np.zeros((24, 20, 24), dtype=np.float32)
    mask = np.zeros_like(density, dtype=np.uint8)
    density[2:22, 4:14, 2:18] = 700.0
    mask[2:22, 4:14, 2:18] = 2

    cropped, crop_face, meta = femur.crop_vtk_images_to_bbox_ratio(
        [_vtk_image_from_array(density)],
        _vtk_image_from_array(mask),
        bbox_ratio=femur.DEFAULT_FEMUR_BBOX_RATIO,
        bbox_crop_from=femur.DEFAULT_FEMUR_BBOX_CROP_FROM,
        labels={2},
    )

    cropped_density = vtk_image_to_numpy(cropped[0])
    crop_face_data = vtk_image_to_numpy(crop_face)

    assert cropped_density.shape == (20, 10, 13)
    assert cropped[0].GetOrigin() == pytest.approx((2, 4, 5))
    assert meta["ratio_xyz"] == (None, 1.0, 1.3)
    assert meta["crop_from_xyz"] == (None, None, "min")
    assert meta["crop_slices_xyz"] == ((2, 22), (4, 14), (5, 18))
    assert crop_face_data.sum() == 20 * 10
    assert np.all(crop_face_data[:, :, 0] == 1)
    assert not np.any(crop_face_data[:, :, 1:])


def test_bbox_ratio_crop_uses_shortest_unit_axis_as_reference():
    np = pytest.importorskip("numpy")
    from ogo.util.vtk_image import vtk_image_to_numpy

    density = np.zeros((24, 20, 24), dtype=np.float32)
    mask = np.zeros_like(density, dtype=np.uint8)
    density[2:22, 4:14, 2:18] = 700.0
    mask[2:22, 4:14, 2:18] = 2

    cropped, _crop_face, meta = femur.crop_vtk_images_to_bbox_ratio(
        [_vtk_image_from_array(density)],
        _vtk_image_from_array(mask),
        bbox_ratio=(1.0, 1.0, 1.0),
        labels={2},
    )

    assert vtk_image_to_numpy(cropped[0]).shape == (10, 10, 10)
    assert cropped[0].GetOrigin() == pytest.approx((7, 4, 5))
    assert meta["reference_axis"] == "y"
    assert meta["reference_length_mm"] == pytest.approx(10.0)


def test_bbox_ratio_crop_recomputes_reference_after_distal_shaft_crop():
    np = pytest.importorskip("numpy")
    from ogo.util.vtk_image import vtk_image_to_numpy

    density = np.zeros((12, 60, 90), dtype=np.float32)
    mask = np.zeros_like(density, dtype=np.uint8)
    density[2:10, 5:55, 0:10] = 700.0
    mask[2:10, 5:55, 0:10] = 2
    density[2:10, 20:40, 10:80] = 700.0
    mask[2:10, 20:40, 10:80] = 2

    cropped, _crop_face, meta = femur.crop_vtk_images_to_bbox_ratio(
        [_vtk_image_from_array(density)],
        _vtk_image_from_array(mask),
        bbox_ratio=femur.DEFAULT_FEMUR_BBOX_RATIO,
        bbox_crop_from=femur.DEFAULT_FEMUR_BBOX_CROP_FROM,
        labels={2},
    )

    cropped_mask = vtk_image_to_numpy(cropped[0]) != 0
    coords = np.argwhere(cropped_mask)
    final_size = coords.max(axis=0) - coords.min(axis=0) + 1

    assert final_size[1] == 20
    assert final_size[2] <= 26
    assert meta["crop_iterations"] > 1
    assert meta["reference_length_mm"] == pytest.approx(20.0)
    assert meta["final_bbox_ratio_xyz"][2] <= 1.3


def test_proximal_box_ratio_crop_uses_proximal_transverse_width():
    np = pytest.importorskip("numpy")
    from ogo.util.vtk_image import vtk_image_to_numpy

    density = np.zeros((110, 120, 120), dtype=np.float32)
    mask = np.zeros_like(density, dtype=np.uint8)
    density[30:82, 45:73, 70:110] = 700.0
    mask[30:82, 45:73, 70:110] = 2
    density[46:66, 45:105, 10:70] = 700.0
    mask[46:66, 45:105, 10:70] = 2

    cropped, crop_face, meta = femur.crop_vtk_images_to_proximal_box_ratio(
        [_vtk_image_from_array(density)],
        _vtk_image_from_array(mask),
        ratio=1.2,
        labels={2},
    )

    cropped_mask = vtk_image_to_numpy(cropped[0]) != 0
    coords = np.argwhere(cropped_mask)
    final_size = coords.max(axis=0) - coords.min(axis=0) + 1

    assert meta["reference_width_mm"] == pytest.approx(52.0)
    assert meta["reference_axis"] == "x"
    assert meta["target_length_mm"] == pytest.approx(62.4)
    assert final_size[2] == 62
    assert final_size[1] == 60
    assert vtk_image_to_numpy(crop_face).sum() > 0


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


def test_post_icp_flat_ratio_crop_uses_y_width_for_z_length():
    np = pytest.importorskip("numpy")
    from ogo.util.vtk_image import vtk_image_to_numpy

    density = np.zeros((32, 80, 160), dtype=np.float32)
    mask = np.zeros_like(density, dtype=np.uint8)
    density[5:25, 10:70, 20:150] = 700.0
    mask[5:25, 10:70, 20:150] = 2

    cropped, crop_face, meta = femur.crop_vtk_images_to_flat_post_icp_ratio(
        [_vtk_image_from_array(density)],
        _vtk_image_from_array(mask),
        ratio=1.2,
        labels={2},
    )

    cropped_mask = vtk_image_to_numpy(cropped[0]) != 0
    coords = np.argwhere(cropped_mask)
    final_size = coords.max(axis=0) - coords.min(axis=0) + 1

    assert meta["method"] == "post_icp_flat_ratio"
    assert meta["reference_width_mm"] == pytest.approx(60.0)
    assert meta["target_length_mm"] == pytest.approx(72.0)
    assert meta["status"] == "cropped"
    assert final_size[1] == 60
    assert final_size[2] == 72
    assert final_size[2] / final_size[1] == pytest.approx(1.2)
    assert vtk_image_to_numpy(crop_face).sum() > 0


def test_post_icp_oblique_ratio_crop_uses_transformed_crop_face_angle():
    np = pytest.importorskip("numpy")
    from ogo.util.vtk_image import vtk_image_to_numpy

    density = np.zeros((32, 80, 160), dtype=np.float32)
    mask = np.zeros_like(density, dtype=np.uint8)
    rough_face = np.zeros_like(mask)
    for x in range(5, 25):
        for y in range(10, 70):
            z = 20 + (x - 5) // 5
            density[x, y, z:150] = 700.0
            mask[x, y, z:150] = 2
            rough_face[x, y, z] = 1

    cropped, crop_face, meta = femur.crop_vtk_images_to_oblique_post_icp_ratio(
        [_vtk_image_from_array(density)],
        _vtk_image_from_array(mask),
        _vtk_image_from_array(rough_face),
        ratio=1.2,
        labels={2},
    )

    cropped_mask = vtk_image_to_numpy(cropped[0]) != 0
    crop_face_mask = vtk_image_to_numpy(crop_face) != 0
    crop_face_coords = np.argwhere(crop_face_mask)

    assert meta["method"] == "post_icp_oblique_ratio"
    assert meta["reference_width_mm"] == pytest.approx(60.0)
    assert meta["target_length_mm"] == pytest.approx(72.0)
    assert meta["status"] == "cropped"
    assert meta["plane_normal"][2] > 0.8
    assert crop_face_coords[:, 2].max() > crop_face_coords[:, 2].min()
    assert cropped_mask.sum() > 0
    assert crop_face_mask.sum() > 0


def test_post_icp_oblique_ratio_short_case_keeps_transformed_crop_face():
    np = pytest.importorskip("numpy")
    from ogo.util.vtk_image import vtk_image_to_numpy

    density = np.zeros((20, 40, 50), dtype=np.float32)
    mask = np.zeros_like(density, dtype=np.uint8)
    rough_face = np.zeros_like(mask)
    density[4:16, 5:35, 12:42] = 700.0
    mask[4:16, 5:35, 12:42] = 2
    rough_face[4:16, 5:35, 12] = 1

    cropped, crop_face, meta = femur.crop_vtk_images_to_oblique_post_icp_ratio(
        [_vtk_image_from_array(density)],
        _vtk_image_from_array(mask),
        _vtk_image_from_array(rough_face),
        ratio=1.2,
        labels={2},
    )

    cropped_mask = vtk_image_to_numpy(cropped[0]) != 0
    crop_face_mask = vtk_image_to_numpy(crop_face) != 0

    assert meta["method"] == "post_icp_oblique_ratio"
    assert meta["status"] == "short"
    assert cropped_mask.sum() == np.count_nonzero(mask)
    assert crop_face_mask.sum() == np.count_nonzero(rough_face)
    assert np.unique(np.argwhere(crop_face_mask)[:, 2]).tolist() == [0]


def test_post_icp_oblique_ratio_allows_resampled_face_offset_from_mask():
    np = pytest.importorskip("numpy")
    from ogo.util.vtk_image import vtk_image_to_numpy

    density = np.zeros((32, 80, 160), dtype=np.float32)
    mask = np.zeros_like(density, dtype=np.uint8)
    rough_face = np.zeros_like(mask)
    density[5:25, 10:70, 20:150] = 700.0
    mask[5:25, 10:70, 20:150] = 2
    rough_face[5:25, 10:70, 19] = 1

    cropped, crop_face, meta = femur.crop_vtk_images_to_oblique_post_icp_ratio(
        [_vtk_image_from_array(density)],
        _vtk_image_from_array(mask),
        _vtk_image_from_array(rough_face),
        ratio=1.2,
        labels={2},
    )

    assert meta["method"] == "post_icp_oblique_ratio"
    assert meta["status"] == "cropped"
    assert vtk_image_to_numpy(cropped[0]).sum() > 0
    assert vtk_image_to_numpy(crop_face).sum() > 0


def test_thicken_z_crop_face_into_mask_returns_mask_constrained_slab():
    np = pytest.importorskip("numpy")
    from ogo.util.vtk_image import vtk_image_to_numpy

    mask = np.zeros((8, 8, 8), dtype=np.uint8)
    mask[2:6, 2:6, 1:7] = 1
    face = np.zeros_like(mask)
    face[2:6, 2:6, 1] = 1

    slab = femur.thicken_z_crop_face_into_mask(
        _vtk_image_from_array(face),
        _vtk_image_from_array(mask),
        thickness_voxels=3,
    )
    slab_data = vtk_image_to_numpy(slab) != 0

    assert slab_data.sum() == 4 * 4 * 3
    assert np.all(slab_data[2:6, 2:6, 1:4])
    assert not np.any(slab_data[:, :, 0])
    assert not np.any(slab_data[:, :, 4:])


def test_thicken_z_crop_face_falls_back_to_distal_mask_slice():
    np = pytest.importorskip("numpy")
    from ogo.util.vtk_image import vtk_image_to_numpy

    mask = np.zeros((8, 8, 8), dtype=np.uint8)
    mask[2:6, 2:6, 2:7] = 1
    empty_face = np.zeros_like(mask)

    slab = femur.thicken_z_crop_face_into_mask(
        _vtk_image_from_array(empty_face),
        _vtk_image_from_array(mask),
        thickness_voxels=2,
    )
    slab_data = vtk_image_to_numpy(slab) != 0

    assert slab_data.sum() == 4 * 4 * 2
    assert np.all(slab_data[2:6, 2:6, 2:4])
    assert not np.any(slab_data[:, :, :2])
    assert not np.any(slab_data[:, :, 4:])


def test_crop_face_support_vector_points_away_from_retained_bone():
    np = pytest.importorskip("numpy")

    mask = np.zeros((10, 10, 10), dtype=np.uint8)
    crop_face = np.zeros_like(mask)
    for x in range(2, 8):
        for y in range(2, 8):
            z = 2 + (x - 2) // 3
            crop_face[x, y, z] = 1
            mask[x, y, z:9] = 1

    vector = femur.crop_face_support_vector(
        _vtk_image_from_array(crop_face),
        _vtk_image_from_array(mask),
    )

    assert np.linalg.norm(vector) == pytest.approx(1.0)
    assert vector[2] < -0.6


def test_crop_face_contact_plane_uses_transformed_crop_face_angle():
    np = pytest.importorskip("numpy")

    mask = np.zeros((12, 12, 12), dtype=np.uint8)
    crop_face = np.zeros_like(mask)
    for x in range(3, 9):
        for y in range(3, 9):
            z = 2 + (x - 3) // 3
            crop_face[x, y, z] = 1
            mask[x, y, z:10] = 1

    plane = femur.crop_face_contact_plane(
        _vtk_image_from_array(crop_face),
        _vtk_image_from_array(mask),
    )

    assert np.linalg.norm(plane["normal"]) == pytest.approx(1.0)
    assert plane["normal"][2] > 0.6
    assert plane["normal"][0] < -0.2
    assert plane["outward_normal"][2] < -0.6
    assert plane["size"][0] > 0
    assert plane["size"][1] > 0


def test_projected_crop_face_surface_mask_keeps_first_contact_layer():
    np = pytest.importorskip("numpy")
    from ogo.util.vtk_image import vtk_image_to_numpy

    active = np.zeros((10, 10, 10), dtype=np.uint8)
    active[3:7, 3:7, 4:8] = 1
    plane = {
        "center": (4.5, 4.5, 3.0),
        "normal": (0.0, 0.0, 1.0),
        "u_axis": (1.0, 0.0, 0.0),
        "v_axis": (0.0, 1.0, 0.0),
        "size": (4.0, 4.0),
        "shape": "anatomy",
    }

    surface = femur.projected_crop_face_surface_vtk(
        _vtk_image_from_array(active),
        plane,
        intrusion=6.0,
    )
    surface_data = vtk_image_to_numpy(surface) != 0

    assert surface_data.sum() == 16
    assert np.all(surface_data[3:7, 3:7, 4])
    assert not np.any(surface_data[:, :, 5:])


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
    assert padded[0].GetOrigin() == pytest.approx((-2, 0, -2))
    assert padded[0].GetExtent() == (0, 6, 0, 5, 0, 8)
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
