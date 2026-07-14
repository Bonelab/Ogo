import pytest


def _vtk_image_from_array(data, *, origin=(0, 0, 0), spacing=(1, 1, 1)):
    vtk = pytest.importorskip("vtk")
    from vtk.util.numpy_support import numpy_to_vtk

    image = vtk.vtkImageData()
    image.SetDimensions(data.shape)
    image.SetOrigin(*origin)
    image.SetSpacing(*spacing)
    image.GetPointData().SetScalars(numpy_to_vtk(data.ravel(order="F"), deep=True))
    return image


def test_surface_points_from_vtk_mask_uses_voxel_surface_centers():
    np = pytest.importorskip("numpy")

    from ogo.fea.alignment import surface_points_from_vtk_mask

    mask = np.ones((3, 3, 3), dtype=np.uint8)
    image = _vtk_image_from_array(mask, origin=(10, 20, 30), spacing=(2, 3, 4))

    full = surface_points_from_vtk_mask(image, max_points=None)
    sampled = surface_points_from_vtk_mask(
        image,
        max_points=5,
        sample_mode="stride",
    )

    assert full.shape == (26, 3)
    assert not np.any(np.all(np.isclose(full, (12, 23, 34)), axis=1))
    assert sampled == pytest.approx(full[[0, 6, 12, 18, 24]])


def test_estimate_rigid_icp_recovers_simple_translation():
    np = pytest.importorskip("numpy")

    from ogo.fea.alignment import estimate_rigid_icp

    moving = np.asarray(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    fixed = moving + np.asarray([3.0, -2.0, 5.0])

    transform = estimate_rigid_icp(
        moving_points=moving,
        fixed_points=fixed,
        iterations=5,
        tolerance=1.0e-9,
    )

    assert transform["rotation"] == pytest.approx(np.eye(3))
    assert transform["translation"] == pytest.approx((3.0, -2.0, 5.0))


def test_invert_point_transform_reverses_row_vector_transform():
    np = pytest.importorskip("numpy")

    from ogo.fea.alignment import invert_point_transform

    rotation = np.asarray(
        [
            [0.0, -1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    translation = np.asarray([3.0, -2.0, 5.0])
    points = np.asarray([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])

    inverse_rotation, inverse_translation = invert_point_transform(rotation, translation)
    transformed = points @ rotation.T + translation
    restored = transformed @ inverse_rotation.T + inverse_translation

    assert restored == pytest.approx(points)


def test_output_grid_for_point_transform_uses_transformed_voxel_centers():
    np = pytest.importorskip("numpy")

    from ogo.fea.alignment import output_grid_for_point_transform

    points = np.asarray(
        [
            [10.0, 20.0, 30.0],
            [14.0, 20.0, 30.0],
            [10.0, 26.0, 30.0],
            [10.0, 20.0, 38.0],
        ]
    )
    rotation = np.eye(3)
    translation = np.asarray([1.0, -2.0, 3.0])

    origin, size = output_grid_for_point_transform(
        points,
        rotation=rotation,
        translation=translation,
        spacing=(2.0, 3.0, 4.0),
        margin_voxels=1,
    )

    assert origin == pytest.approx((9.0, 15.0, 29.0))
    assert size == (5, 5, 5)


def test_resample_vtk_image_with_point_transform_uses_explicit_output_grid():
    np = pytest.importorskip("numpy")

    from ogo.fea.alignment import resample_vtk_image_with_point_transform
    from ogo.util.vtk_image import vtk_image_to_numpy

    data = np.zeros((3, 3, 3), dtype=np.uint8)
    data[1, 1, 1] = 7
    image = _vtk_image_from_array(data, origin=(0, 0, 0), spacing=(1, 1, 1))

    shifted = resample_vtk_image_with_point_transform(
        image,
        rotation=np.eye(3),
        translation=np.asarray([2.0, 0.0, 0.0]),
        output_origin=(0.0, 0.0, 0.0),
        output_size=(5, 3, 3),
        output_spacing=(1.0, 1.0, 1.0),
        interpolation="nearest",
    )

    shifted_array = vtk_image_to_numpy(shifted)
    assert shifted.GetOrigin() == pytest.approx((0.0, 0.0, 0.0))
    assert shifted.GetSpacing() == pytest.approx((1.0, 1.0, 1.0))
    assert shifted.GetDimensions() == (5, 3, 3)
    assert shifted_array[3, 1, 1] == 7
    assert int(np.count_nonzero(shifted_array)) == 1
