"""Shared point-cloud helpers for reference-frame FE alignment."""


def polydata_points(polydata):
    """Return VTK polydata point coordinates as an ``(n, 3)`` NumPy array."""
    import numpy as np
    from vtk.util.numpy_support import vtk_to_numpy

    points = polydata.GetPoints()
    if points is None or points.GetNumberOfPoints() == 0:
        raise ValueError("Polydata contains no points.")
    coordinates = np.asarray(vtk_to_numpy(points.GetData()), dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError("Polydata points must have shape (n, 3).")
    return coordinates[np.all(np.isfinite(coordinates), axis=1)]


def polydata_from_points(points):
    """Create a vertex-only VTK polydata from an ``(n, 3)`` point cloud."""
    import numpy as np
    import vtk

    coordinates = np.asarray(points, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3 or coordinates.shape[0] < 3:
        raise ValueError("points must have shape (n, 3) with at least three rows.")

    vtk_points = vtk.vtkPoints()
    vertices = vtk.vtkCellArray()
    for point in coordinates:
        point_id = vtk_points.InsertNextPoint(float(point[0]), float(point[1]), float(point[2]))
        vertices.InsertNextCell(1)
        vertices.InsertCellPoint(point_id)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(vtk_points)
    polydata.SetVerts(vertices)
    return polydata


def sample_points(points, *, max_points=None, mode="linspace", offset=0):
    """Sample a point cloud using deterministic stride or linspace modes."""
    import numpy as np

    coordinates = np.asarray(points, dtype=float)
    if max_points is None or int(max_points) <= 0 or coordinates.shape[0] <= int(max_points):
        return coordinates

    token = str(mode).strip().lower()
    if token == "stride":
        step = max(1, int(np.ceil(coordinates.shape[0] / int(max_points))))
        start = int(offset) % step
        indices = np.arange(start, coordinates.shape[0], step, dtype=int)
        if indices.size < int(max_points):
            fallback = np.linspace(0, coordinates.shape[0] - 1, int(max_points), dtype=int)
            indices = np.unique(np.concatenate([indices, fallback]))
        indices = indices[: int(max_points)]
    else:
        indices = np.linspace(0, coordinates.shape[0] - 1, int(max_points), dtype=int)
    return coordinates[indices]


def _binary_erosion_6(mask):
    import numpy as np

    padded = np.pad(mask, 1, mode="constant", constant_values=False)
    center = padded[1:-1, 1:-1, 1:-1]
    eroded = center.copy()
    eroded &= padded[:-2, 1:-1, 1:-1]
    eroded &= padded[2:, 1:-1, 1:-1]
    eroded &= padded[1:-1, :-2, 1:-1]
    eroded &= padded[1:-1, 2:, 1:-1]
    eroded &= padded[1:-1, 1:-1, :-2]
    eroded &= padded[1:-1, 1:-1, 2:]
    return eroded


def surface_points_from_vtk_mask(
    vtk_mask,
    *,
    max_points=8000,
    sample_mode="stride",
    sample_offset=0,
):
    """Return physical x/y/z points from 6-connected voxel surface centers."""
    import numpy as np

    from ogo.util.vtk_image import vtk_image_to_numpy

    mask = np.asarray(vtk_image_to_numpy(vtk_mask), dtype=bool)
    if not np.any(mask):
        raise ValueError("mask contains no foreground voxels.")

    surface = mask & ~_binary_erosion_6(mask)
    surface_zyx = np.transpose(surface if np.any(surface) else mask, (2, 1, 0))
    coords_xyz = np.argwhere(surface_zyx)[:, [2, 1, 0]].astype(float)
    points = np.asarray(vtk_mask.GetOrigin(), dtype=float) + coords_xyz * np.asarray(
        vtk_mask.GetSpacing(),
        dtype=float,
    )
    return sample_points(
        points,
        max_points=max_points,
        mode=sample_mode,
        offset=sample_offset,
    )


def point_cloud_axis_lengths(points):
    """Return approximate principal-axis diameters for an ``(n, 3)`` point cloud."""
    import numpy as np

    coordinates = np.asarray(points, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3 or coordinates.shape[0] < 3:
        raise ValueError("points must have shape (n, 3) with at least three rows.")
    centered = coordinates - coordinates.mean(axis=0)
    covariance = np.cov(centered.T)
    eigvals = np.linalg.eigvalsh(covariance)
    return np.sqrt(np.maximum(eigvals, 0.0)) * 2.0


def estimate_rigid_icp(
    *,
    moving_points,
    fixed_points,
    iterations=50,
    tolerance=1.0e-4,
    start_by_matching_centroids_only=False,
    convergence="delta",
    distance_mode="mean",
    initial_transform=None,
):
    """Estimate a rigid transform from ``moving_points`` to ``fixed_points``."""
    import numpy as np

    moving = _points_array(moving_points, "moving_points")
    fixed = _points_array(fixed_points, "fixed_points")
    if initial_transform is not None:
        rotation = np.asarray(initial_transform["rotation"], dtype=float)
        translation = np.asarray(initial_transform["translation"], dtype=float)
        if rotation.shape != (3, 3):
            raise ValueError("initial_transform rotation must have shape (3, 3).")
        if translation.shape != (3,):
            raise ValueError("initial_transform translation must contain three values.")
    elif start_by_matching_centroids_only:
        rotation = np.eye(3)
        translation = fixed.mean(axis=0) - moving.mean(axis=0)
    else:
        rotation, translation = _best_initial_transform(moving, fixed)

    previous_error = np.inf
    used_iterations = 0
    convergence_token = str(convergence).strip().lower()
    distance_token = str(distance_mode).strip().lower()
    for used_iterations in range(1, max(1, int(iterations)) + 1):
        transformed = moving @ rotation.T + translation
        matched = fixed[_nearest_indices(transformed, fixed)]
        step_rotation, step_translation = _kabsch(transformed, matched)
        rotation = step_rotation @ rotation
        translation = step_rotation @ translation + step_translation
        distances = np.linalg.norm(transformed - matched, axis=1)
        if distance_token == "mean":
            error = float(np.mean(distances))
        elif distance_token == "rms":
            error = float(np.sqrt(np.mean(distances * distances)))
        else:
            raise ValueError("distance_mode must be 'mean' or 'rms'.")

        if convergence_token == "delta":
            converged = abs(previous_error - error) <= float(tolerance)
        elif convergence_token in {"absolute", "abs"}:
            converged = error <= float(tolerance)
        else:
            raise ValueError("convergence must be 'delta' or 'absolute'.")
        if converged:
            previous_error = error
            break
        previous_error = error

    return {
        "rotation": rotation,
        "translation": translation,
        "iterations": used_iterations,
        "mean_distance": previous_error,
    }


def invert_point_transform(rotation, translation):
    """Invert ``points @ rotation.T + translation`` row-vector transforms."""
    import numpy as np

    rotation_arr = np.asarray(rotation, dtype=float)
    translation_arr = np.asarray(translation, dtype=float)
    if rotation_arr.shape != (3, 3):
        raise ValueError("rotation must have shape (3, 3).")
    if translation_arr.shape != (3,):
        raise ValueError("translation must contain three values.")
    inverse_rotation = rotation_arr.T
    inverse_translation = -translation_arr @ rotation_arr
    return inverse_rotation, inverse_translation


def output_grid_for_point_transform(
    points,
    *,
    rotation,
    translation,
    spacing,
    margin_voxels=0,
):
    """Return an explicit output image grid for transformed physical points."""
    import numpy as np

    coordinates = _points_array(points, "points")
    rotation_arr = np.asarray(rotation, dtype=float)
    translation_arr = np.asarray(translation, dtype=float)
    spacing_arr = np.asarray(spacing, dtype=float)
    if rotation_arr.shape != (3, 3):
        raise ValueError("rotation must have shape (3, 3).")
    if translation_arr.shape != (3,):
        raise ValueError("translation must contain three values.")
    if spacing_arr.shape != (3,) or np.any(spacing_arr <= 0.0):
        raise ValueError("spacing must contain three positive values.")

    transformed = coordinates @ rotation_arr.T + translation_arr
    margin = int(max(0, margin_voxels))
    lower = transformed.min(axis=0) - margin * spacing_arr
    upper = transformed.max(axis=0) + margin * spacing_arr
    size = np.maximum(1, np.ceil((upper - lower) / spacing_arr).astype(int) + 1)
    return tuple(float(value) for value in lower), tuple(int(value) for value in size)


def resample_vtk_image_with_point_transform(
    vtk_image,
    *,
    rotation,
    translation,
    output_origin,
    output_size,
    output_spacing,
    interpolation="nearest",
):
    """Resample a VTK image onto an explicit grid after a point transform."""
    import numpy as np
    import SimpleITK as sitk
    import vtk
    from vtk.util.numpy_support import numpy_to_vtk

    from ogo.util.vtk_image import vtk_image_to_numpy

    rotation_arr = np.asarray(rotation, dtype=float)
    translation_arr = np.asarray(translation, dtype=float)
    if rotation_arr.shape != (3, 3):
        raise ValueError("rotation must have shape (3, 3).")
    if translation_arr.shape != (3,):
        raise ValueError("translation must contain three values.")

    array_xyz = vtk_image_to_numpy(vtk_image)
    image = sitk.GetImageFromArray(np.transpose(array_xyz, (2, 1, 0)))
    image.SetSpacing(tuple(float(value) for value in vtk_image.GetSpacing()))
    image.SetOrigin(tuple(float(value) for value in vtk_image.GetOrigin()))

    transform = sitk.AffineTransform(3)
    transform.SetMatrix(tuple(float(value) for value in rotation_arr.reshape(-1)))
    transform.SetTranslation(tuple(float(value) for value in translation_arr))

    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(image)
    resampler.SetSize([int(value) for value in output_size])
    resampler.SetOutputSpacing(tuple(float(value) for value in output_spacing))
    resampler.SetOutputOrigin(tuple(float(value) for value in output_origin))
    resampler.SetOutputDirection((1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
    resampler.SetTransform(transform.GetInverse())
    mode = str(interpolation).strip().lower()
    if mode == "nearest":
        resampler.SetInterpolator(sitk.sitkNearestNeighbor)
    elif mode in {"linear", "trilinear"}:
        resampler.SetInterpolator(sitk.sitkLinear)
    elif mode in {"cubic", "bspline", "spline"}:
        resampler.SetInterpolator(sitk.sitkBSpline)
    else:
        raise ValueError("interpolation must be nearest, linear, or bspline.")
    resampler.SetDefaultPixelValue(0)
    resampled_zyx = sitk.GetArrayFromImage(resampler.Execute(image))
    resampled_xyz = np.transpose(resampled_zyx, (2, 1, 0))

    output = vtk.vtkImageData()
    output.SetDimensions(tuple(int(value) for value in output_size))
    output.SetSpacing(tuple(float(value) for value in output_spacing))
    output.SetOrigin(tuple(float(value) for value in output_origin))
    output.GetPointData().SetScalars(
        numpy_to_vtk(resampled_xyz.ravel(order="F"), deep=True)
    )
    return output


def _points_array(points, name):
    import numpy as np

    array = np.asarray(points, dtype=float)
    if array.ndim != 2 or array.shape[1] != 3:
        raise ValueError(f"{name} must have shape (n, 3).")
    if array.shape[0] < 3:
        raise ValueError(f"{name} must contain at least three points.")
    return array[np.all(np.isfinite(array), axis=1)]


def _best_initial_transform(moving, fixed):
    import numpy as np

    moving_center = moving.mean(axis=0)
    fixed_center = fixed.mean(axis=0)
    candidates = [
        (np.eye(3), fixed_center - moving_center),
        _initial_pca_transform(moving, fixed),
    ]
    return min(candidates, key=lambda item: _nearest_mean_distance(moving, fixed, *item))


def _initial_pca_transform(moving, fixed):
    import numpy as np

    moving_center = moving.mean(axis=0)
    fixed_center = fixed.mean(axis=0)
    moving_axes = _principal_axes(moving - moving_center)
    fixed_axes = _principal_axes(fixed - fixed_center)
    rotation = fixed_axes @ moving_axes.T
    if np.linalg.det(rotation) < 0:
        fixed_axes[:, -1] *= -1
        rotation = fixed_axes @ moving_axes.T
    translation = fixed_center - rotation @ moving_center
    return rotation, translation


def _principal_axes(points):
    import numpy as np

    _, _, vh = np.linalg.svd(points, full_matrices=False)
    axes = vh.T
    if np.linalg.det(axes) < 0:
        axes[:, -1] *= -1
    return axes


def _nearest_indices(query, target):
    try:
        from scipy.spatial import cKDTree
    except ImportError:
        return _nearest_indices_numpy(query, target)
    return cKDTree(target).query(query, workers=-1)[1]


def _nearest_indices_numpy(query, target):
    import numpy as np

    out = np.empty((query.shape[0],), dtype=int)
    chunk = 512
    for start in range(0, query.shape[0], chunk):
        stop = min(start + chunk, query.shape[0])
        diff = query[start:stop, None, :] - target[None, :, :]
        out[start:stop] = np.argmin(np.sum(diff * diff, axis=2), axis=1)
    return out


def _nearest_mean_distance(moving, fixed, rotation, translation):
    import numpy as np

    transformed = moving @ rotation.T + translation
    matched = fixed[_nearest_indices(transformed, fixed)]
    return float(np.mean(np.linalg.norm(transformed - matched, axis=1)))


def _kabsch(moving, fixed):
    import numpy as np

    moving_center = moving.mean(axis=0)
    fixed_center = fixed.mean(axis=0)
    h = (moving - moving_center).T @ (fixed - fixed_center)
    u, _, vt = np.linalg.svd(h)
    rotation = vt.T @ u.T
    if np.linalg.det(rotation) < 0:
        vt[-1, :] *= -1
        rotation = vt.T @ u.T
    translation = fixed_center - rotation @ moving_center
    return rotation, translation
