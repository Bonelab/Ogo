"""Femur-specific defaults and helpers for sideways-fall FE generation.

Default femur workflow:
1. Read the calibrated density image and whole-femur mask. The public wrapper
   chooses left/right from ``--side`` and the lower-level script uses
   ``femur_side=1`` for left and ``2`` for right.
2. Run ICP to the bundled side-specific femur reference from the cropped and
   padded input geometry, then build an explicit reference-frame output grid
   from the transformed femur surface. Density is resampled with B-spline
   interpolation; bone, crop-face, and compartment labels use nearest-neighbor
   interpolation on the same grid.
3. Smooth the transformed femur mask with one binary close/open pass only when
   at least one input spacing dimension is coarser than 2 mm. If a compartment
   mask is supplied, the derived cortical binary mask follows the same rule.
4. Standardize the distal shaft by cutting on a flat model-grid z plane. The
   default cut mode detects the lesser-trochanter cross-section peak and keeps
   the femur to 50 mm distal to that landmark. If the scan does not include the
   required distal field of view, model generation fails instead of silently
   using a shorter shaft. ``fixed_length`` mode is available for debugging.
5. Generate two geometric PMMA fixtures from bbox-relative contact planes:
   a femoral-head loading fixture on the high-y side and a greater-trochanter
   contact fixture on the low-y side. Defaults are 10 mm PMMA thickness and
   6 mm intrusion through that fixed thickness. The square fixture footprints
   scale with the generated model bbox, and the fixture masks themselves do
   not overwrite bone voxels.
6. Apply sideways-fall boundary conditions: prescribed displacement at the
   femoral-head PMMA cap toward the greater trochanter, loading-direction
   constraint at the greater-trochanter PMMA cap, and distal shaft constraints
   to remove rigid-body motion.
7. Build materials with the same shared bone/PMMA material-table helper used by
   the spine workflow. If no compartment mask is supplied, femur bone is one
   trabecular-style region. If supplied, cortical=1 and trabecular=2 by default.
"""

from pathlib import Path


LEFT_FEMUR = 1
RIGHT_FEMUR = 2

DEFAULT_FEMUR_ISO_RESOLUTION_MM = 1.0
DEFAULT_FEMUR_FE_DISPLACEMENT = -4.0
DEFAULT_FEMUR_TARGET_DISPLACEMENT_PERCENT = 4.0
DEFAULT_FEMUR_MASK_SMOOTHING_SPACING_THRESHOLD_MM = 2.0
FEMORAL_HEAD_FIXTURE_WIDTH_EXTENSION_MM = 10.0
FEMORAL_HEAD_FIXTURE_LONG_AXIS_EXTENSION_MM = 80.0
SIDEWAYS_FALL_FIXTURE_SIZE_FRACTION = (1.1, 1.1)
FEMORAL_HEAD_FIXTURE_CENTER_FRACTION = (0.5, 1.0278767152, 0.9650933279)
GREATER_TROCHANTER_FIXTURE_CENTER_FRACTION = (0.5, 0.0151983839, 0.6551527108)
DISTAL_SHAFT_FIXTURE_CENTER_FRACTION = (0.1958514895, 0.5708725889, 0.0601674635)
DISTAL_SHAFT_FIXTURE_SIZE_FRACTION = (0.3640088821, 0.6351022953)
DEFAULT_FEMUR_BBOX_RATIO = (1.0, 1.2, None)
DEFAULT_FEMUR_BBOX_CROP_FROM = (None, "min", None)
DEFAULT_FEMUR_REFERENCE_MIN_SCALE = (0.8, 0.8, 0.75)
DEFAULT_FEMUR_REFERENCE_MAX_SCALE = (1.2, 1.2, 1.3)
DEFAULT_PMMA_THICKNESS_MM = 10.0
DEFAULT_PMMA_INTRUSION_MM = 6.0
DEFAULT_FEMUR_INPUT_MARGIN_MM = DEFAULT_PMMA_THICKNESS_MM + DEFAULT_PMMA_INTRUSION_MM
DEFAULT_FEMUR_SHAFT_LENGTH_MM = 100.0
DEFAULT_FEMUR_CUT_MODE = "bbox_ratio"
DEFAULT_LESSER_TROCHANTER_DISTAL_OFFSET_MM = 50.0
DEFAULT_CORTICAL_LABEL = 1
DEFAULT_TRABECULAR_LABEL = 2

FEMORAL_HEAD_NODE_SET = "Femoral_Head_PMMA_Nodes"
GREATER_TROCHANTER_NODE_SET = "Greater_Trochanter_PMMA_Nodes"
DISTAL_FEMUR_NODE_SET = "Distal_Femur_Nodes"
SIDEWAYS_FALL_NODE_SETS = [
    FEMORAL_HEAD_NODE_SET,
    GREATER_TROCHANTER_NODE_SET,
    DISTAL_FEMUR_NODE_SET,
]


def side_suffix(femur_side):
    """Return the compact output stem for a femur side."""
    if femur_side == LEFT_FEMUR:
        return "LF"
    if femur_side == RIGHT_FEMUR:
        return "RF"
    raise ValueError("femur_side must be 1 for left or 2 for right.")


def sideways_fall_output_name(output_file, femur_side):
    """Return the compact side-specific sideways-fall output path."""
    output_path = Path(output_file)
    return str(output_path.with_name(f"{output_path.stem}_{side_suffix(femur_side)}.n88model"))


def _axis_index(axis):
    axis_map = {"x": 0, "y": 1, "z": 2}
    try:
        return axis_map[str(axis).lower()]
    except KeyError as exc:
        raise ValueError("projection_axis must be 'x', 'y', or 'z'.") from exc


def _axis_bounds(model_bounds, axis):
    index = 2 * int(axis)
    return float(model_bounds[index]), float(model_bounds[index + 1])


def foreground_voxel_center_bounds_from_mask(mask, *, origin, spacing):
    """Return x/y/z physical bounds of nonzero voxel centers.

    Bbox-relative workflow fixtures are authored against the occupied voxel
    centers, matching the image-space convention used by interactive workflow
    replay. FE mesh/node bounds are wider by the element faces and should not
    be used to resolve those recipe planes.
    """
    import numpy as np

    active = np.asarray(mask) != 0
    coords = np.argwhere(active)
    if coords.size == 0:
        raise ValueError("Cannot compute foreground bounds from an empty mask.")
    origin = np.asarray(origin, dtype=np.float64)
    spacing = np.asarray(spacing, dtype=np.float64)
    if origin.shape != (3,) or spacing.shape != (3,):
        raise ValueError("origin and spacing must contain three values.")
    lo = origin + coords.min(axis=0).astype(np.float64) * spacing
    hi = origin + coords.max(axis=0).astype(np.float64) * spacing
    return (
        float(lo[0]),
        float(hi[0]),
        float(lo[1]),
        float(hi[1]),
        float(lo[2]),
        float(hi[2]),
    )


def foreground_voxel_center_bounds(vtk_image):
    """Return physical foreground bounds for a VTK image using voxel centers."""
    from ogo.util.vtk_image import vtk_image_to_numpy

    return foreground_voxel_center_bounds_from_mask(
        vtk_image_to_numpy(vtk_image),
        origin=vtk_image.GetOrigin(),
        spacing=vtk_image.GetSpacing(),
    )


def matrix4x4_to_numpy(matrix):
    """Return a 4 x 4 NumPy matrix from a VTK matrix or array-like value."""
    import numpy as np

    if hasattr(matrix, "GetElement"):
        return np.asarray(
            [[float(matrix.GetElement(row, col)) for col in range(4)] for row in range(4)],
            dtype=np.float64,
        )
    values = np.asarray(matrix, dtype=np.float64)
    if values.shape != (4, 4):
        raise ValueError("transform matrix must have shape 4 x 4.")
    return values


def reference_grid_from_output_to_input_matrix(
    points_xyz,
    output_to_input_matrix,
    *,
    spacing,
    margin_voxels=4,
):
    """Return an explicit output grid for a transformed surface point cloud.

    Registration matrices used by image resampling map output coordinates back
    to input coordinates. To build the output lattice, first map input surface
    points through the inverse matrix, then pad the transformed bounds by a
    fixed number of voxels. This keeps the final model grid stable and
    independent of toolkit-specific automatic crop rules.
    """
    import numpy as np

    points = np.asarray(points_xyz, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 3 or points.shape[0] == 0:
        raise ValueError("points_xyz must have shape (n, 3) with at least one point.")
    matrix = matrix4x4_to_numpy(output_to_input_matrix)
    inverse = np.linalg.inv(matrix)
    homogeneous = np.column_stack([points, np.ones(points.shape[0], dtype=np.float64)])
    transformed = (homogeneous @ inverse.T)[:, :3]
    spacing_arr = np.asarray(spacing, dtype=np.float64)
    if spacing_arr.shape != (3,) or np.any(spacing_arr <= 0.0):
        raise ValueError("spacing must contain three positive values.")
    margin = max(0, int(margin_voxels))
    lower = transformed.min(axis=0) - margin * spacing_arr
    upper = transformed.max(axis=0) + margin * spacing_arr
    size = np.maximum(1, np.ceil((upper - lower) / spacing_arr).astype(int) + 1)
    return tuple(float(value) for value in lower), tuple(int(value) for value in size)


def reference_grid_from_vtk_mask(
    vtk_mask,
    output_to_input_matrix,
    *,
    margin_voxels=4,
):
    """Return an explicit reference-frame output grid for a transformed mask."""
    points = surface_points_from_vtk_mask(vtk_mask, max_points=None)
    return reference_grid_from_output_to_input_matrix(
        points,
        output_to_input_matrix,
        spacing=vtk_mask.GetSpacing(),
        margin_voxels=margin_voxels,
    )


def transform_resample_vtk_image_to_reference_grid(
    vtk_image,
    output_to_input_matrix,
    *,
    output_origin,
    output_size,
    output_spacing,
    interpolation="nearest",
):
    """Resample a VTK image onto an explicit reference-frame grid."""
    import numpy as np
    import SimpleITK as sitk
    import vtk

    from ogo.util.vtk_image import vtk_image_to_numpy

    matrix = matrix4x4_to_numpy(output_to_input_matrix)
    array_zyx = vtk_image_to_numpy(vtk_image, processing_order=True)
    image = sitk.GetImageFromArray(array_zyx)
    image.SetSpacing(tuple(float(value) for value in vtk_image.GetSpacing()))
    image.SetOrigin(tuple(float(value) for value in vtk_image.GetOrigin()))
    image.SetDirection((1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))

    transform = sitk.AffineTransform(3)
    transform.SetMatrix(tuple(float(value) for value in matrix[:3, :3].reshape(-1)))
    transform.SetTranslation(tuple(float(value) for value in matrix[:3, 3]))

    resampler = sitk.ResampleImageFilter()
    resampler.SetSize([int(value) for value in output_size])
    resampler.SetOutputSpacing(tuple(float(value) for value in output_spacing))
    resampler.SetOutputOrigin(tuple(float(value) for value in output_origin))
    resampler.SetOutputDirection((1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
    resampler.SetTransform(transform)
    resampler.SetDefaultPixelValue(0)
    resampler.SetInterpolator(
        sitk.sitkNearestNeighbor
        if interpolation == "nearest"
        else sitk.sitkBSpline
        if interpolation == "bspline"
        else sitk.sitkLinear
    )
    out_zyx = sitk.GetArrayFromImage(resampler.Execute(image))
    out_xyz = np.transpose(out_zyx, (2, 1, 0))
    vtk_type = vtk.VTK_UNSIGNED_CHAR if interpolation == "nearest" else vtk_image.GetScalarType()
    out = _vtk_image_from_array(
        out_xyz,
        vtk_image,
        origin=tuple(float(value) for value in output_origin),
        vtk_array_type=vtk_type,
    )
    out.SetSpacing(tuple(float(value) for value in output_spacing))
    return out


def bbox_relative_fixture_bounds(
    model_bounds,
    *,
    center_fraction,
    size_fraction,
    projection_axis="y",
    shape="rectangle",
):
    """Return physical contact ROI bounds from model-bbox-relative plane values.

    ``center_fraction`` has x/y/z fractions relative to the model bbox. The two
    ``size_fraction`` values scale the lateral axes of the projection plane in
    coordinate order. When ``shape="square"``, the smaller scaled lateral
    length is used on both axes. The projection axis itself spans the whole
    model bbox so the cap builder can find the nearest supported anatomy
    column.
    """
    from ogo.fea.boundary import bbox_relative_contact_bounds

    return bbox_relative_contact_bounds(
        model_bounds,
        center_fraction=center_fraction,
        size_fraction=size_fraction,
        projection_axis=projection_axis,
        shape=shape,
    )


def bbox_relative_fixture_direction(center_fraction, *, projection_axis="y"):
    """Return the cap side implied by a bbox-relative plane fraction."""
    from ogo.fea.boundary import bbox_relative_contact_direction

    return bbox_relative_contact_direction(
        center_fraction,
        projection_axis=projection_axis,
    )


def _axis_extent_from_vector(model_bounds, vector):
    import numpy as np

    spans = np.asarray(
        [
            float(model_bounds[1]) - float(model_bounds[0]),
            float(model_bounds[3]) - float(model_bounds[2]),
            float(model_bounds[5]) - float(model_bounds[4]),
        ],
        dtype=float,
    )
    return float(np.sum(np.abs(np.asarray(vector, dtype=float)) * spans))


def _fixture_plane_axes(projection_axis, normal_sign):
    """Return stable in-plane axes for bbox-relative fixture planes."""
    axis = _axis_index(projection_axis)
    sign = 1.0 if float(normal_sign) >= 0.0 else -1.0
    if axis == 0:
        return (0.0, 0.0, sign), (0.0, -1.0, 0.0)
    if axis == 1:
        return (0.0, 0.0, -sign), (-1.0, 0.0, 0.0)
    return (1.0, 0.0, 0.0), (0.0, sign, 0.0)


def bbox_relative_fixture_plane(
    model_bounds,
    *,
    center_fraction,
    size_fraction,
    projection_axis="y",
    shape="square",
):
    """Return a physical contact-plane definition from bbox-relative values."""
    from ogo.fea.boundary import bbox_relative_contact_plane

    return bbox_relative_contact_plane(
        model_bounds,
        center_fraction=center_fraction,
        size_fraction=size_fraction,
        projection_axis=projection_axis,
        shape=shape,
    )


def _scale_triplet(values, name):
    import numpy as np

    parsed = np.asarray(values, dtype=float)
    if parsed.shape != (3,):
        raise ValueError(f"{name} must contain three x/y/z values.")
    return parsed


def principal_axis_lengths(polydata):
    """Return approximate principal-axis diameters for a VTK polydata surface."""
    import numpy as np
    from vtk.util.numpy_support import vtk_to_numpy

    points = polydata.GetPoints()
    if points is None or points.GetNumberOfPoints() == 0:
        raise ValueError("Cannot measure principal axes for empty polydata.")
    coordinates = np.asarray(vtk_to_numpy(points.GetData()), dtype=float)
    centered = coordinates - coordinates.mean(axis=0)
    covariance = np.cov(centered.T)
    eigvals = np.linalg.eigvalsh(covariance)
    return np.sqrt(np.maximum(eigvals, 0.0)) * 2.0


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
    """Sample a point cloud using the same small deterministic modes as workflows."""
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
    """Return physical x/y/z points from the 6-connected voxel mask surface."""
    import numpy as np
    from ogo.util.vtk_image import vtk_image_to_numpy

    mask = np.asarray(vtk_image_to_numpy(vtk_mask), dtype=bool)
    if not np.any(mask):
        raise ValueError("mask contains no foreground voxels.")

    surface = mask & ~_binary_erosion_6(mask)
    # Match the workflow-replay path: enumerate the voxel surface in z/y/x
    # NumPy order, then convert those indices back to physical x/y/z points.
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


def scale_reference_point_cloud_to_sample(
    reference_polydata,
    sample_surface_points,
    *,
    max_points=8000,
    reference_sample_mode="linspace",
    min_scale=DEFAULT_FEMUR_REFERENCE_MIN_SCALE,
    max_scale=DEFAULT_FEMUR_REFERENCE_MAX_SCALE,
):
    """Scale sampled reference points to the sampled voxel-surface point cloud."""
    import numpy as np

    reference_points = sample_points(
        polydata_points(reference_polydata),
        max_points=max_points,
        mode=reference_sample_mode,
    )
    sample_points_array = np.asarray(sample_surface_points, dtype=float)
    reference_lengths = point_cloud_axis_lengths(reference_points)
    sample_lengths = point_cloud_axis_lengths(sample_points_array)
    scale = sample_lengths / np.maximum(reference_lengths, 1.0e-6)
    min_values = _scale_triplet(min_scale, "min_scale")
    max_values = _scale_triplet(max_scale, "max_scale")
    scale = np.clip(scale, min_values, max_values)

    scaled_points = reference_points * scale
    return polydata_from_points(scaled_points), {
        "source": "voxel_surface_point_cloud",
        "reference_axis_lengths": reference_lengths.tolist(),
        "sample_axis_lengths": sample_lengths.tolist(),
        "scale_factors": scale.tolist(),
        "min_scale": min_values.tolist(),
        "max_scale": max_values.tolist(),
    }


def scale_reference_to_sample_principal_lengths(
    reference_polydata,
    sample_polydata,
    *,
    min_scale=DEFAULT_FEMUR_REFERENCE_MIN_SCALE,
    max_scale=DEFAULT_FEMUR_REFERENCE_MAX_SCALE,
):
    """Scale a femur reference so its principal lengths match the sample."""
    import numpy as np
    import vtk

    reference_lengths = principal_axis_lengths(reference_polydata)
    sample_lengths = principal_axis_lengths(sample_polydata)
    scale = sample_lengths / np.maximum(reference_lengths, 1.0e-6)
    min_values = _scale_triplet(min_scale, "min_scale")
    max_values = _scale_triplet(max_scale, "max_scale")
    scale = np.clip(scale, min_values, max_values)

    transform = vtk.vtkTransform()
    transform.Scale(float(scale[0]), float(scale[1]), float(scale[2]))

    transform_filter = vtk.vtkTransformPolyDataFilter()
    transform_filter.SetInputData(reference_polydata)
    transform_filter.SetTransform(transform)
    transform_filter.Update()

    return transform_filter.GetOutput(), {
        "reference_axis_lengths": reference_lengths.tolist(),
        "sample_axis_lengths": sample_lengths.tolist(),
        "scale_factors": scale.tolist(),
        "min_scale": min_values.tolist(),
        "max_scale": max_values.tolist(),
    }


def _parse_bbox_ratio_recipe(values):
    """Return bbox-ratio values in Ogo's x/y/z image order.

    Saved workflows store bbox ratios in recipe order:
    ``reference, constrained, free``. For the hip sideways-fall recipe the
    reference axis is y, the constrained crop axis is z, and x is left free.
    Ogo image arrays are x/y/z, so the mapped order is ``free, reference,
    constrained``.
    """
    if len(values) != 3:
        raise ValueError("bbox_ratio must contain reference/constrained/free values.")
    parsed = []
    for item in values:
        if item is None:
            parsed.append(None)
            continue
        token = str(item).strip().lower()
        if token in {"", "none", "null", "auto"}:
            parsed.append(None)
            continue
        value = float(item)
        if value <= 0:
            raise ValueError("bbox_ratio values must be positive or null.")
        parsed.append(value)
    reference, constrained, free = parsed
    return free, reference, constrained


def _crop_from_value(value):
    if value is None:
        return None
    token = str(value).strip().lower()
    if token in {"", "none", "null", "auto", "center", "centre"}:
        return None
    if token in {"min", "low", "lo", "start"}:
        return "min"
    if token in {"max", "high", "hi", "end"}:
        return "max"
    raise ValueError("bbox_crop_from values must be min, max, center, or null.")


def _parse_bbox_crop_from_recipe(values):
    """Return bbox crop-end values in Ogo's x/y/z image order.

    The recipe uses the same reference/constrained/free convention as the
    Slicer-authored workflow. For the hip sideways-fall recipe the mapped order
    is ``free, reference, constrained`` in Ogo's x/y/z image axes.
    """
    if values is None:
        return None, None, None
    if len(values) != 3:
        raise ValueError("bbox_crop_from must contain reference/constrained/free values.")
    reference, constrained, free = (_crop_from_value(item) for item in values)
    return free, reference, constrained


def _vtk_image_from_array(array, template_vtk_image, *, origin, vtk_array_type=None):
    """Create a zero-based VTK image from an x/y/z NumPy array."""
    import numpy as np
    import vtk
    from vtk.util.numpy_support import numpy_to_vtk

    data = np.ascontiguousarray(array)
    image = vtk.vtkImageData()
    image.SetDimensions(data.shape)
    image.SetOrigin(*origin)
    image.SetSpacing(*template_vtk_image.GetSpacing())
    if vtk_array_type is None:
        vtk_array_type = template_vtk_image.GetScalarType()
    scalars = numpy_to_vtk(data.ravel(order="F"), deep=True, array_type=vtk_array_type)
    image.GetPointData().SetScalars(scalars)
    return image


def resample_vtk_image_like_workflow(vtk_image, target_spacing_mm, *, interpolation="nearest"):
    """Resample a VTK image with the same output-size rule as workflow replay."""
    import numpy as np
    import SimpleITK as sitk
    import vtk
    from ogo.util.vtk_image import vtk_image_to_numpy

    target = float(target_spacing_mm)
    target_spacing = (target, target, target)
    array_zyx = vtk_image_to_numpy(vtk_image, processing_order=True)
    image = sitk.GetImageFromArray(array_zyx)
    image.SetSpacing(tuple(float(value) for value in vtk_image.GetSpacing()))
    image.SetOrigin(tuple(float(value) for value in vtk_image.GetOrigin()))
    original_size = np.asarray(image.GetSize(), dtype=np.int64)
    original_spacing = np.asarray(image.GetSpacing(), dtype=np.float64)
    new_spacing = np.asarray(target_spacing, dtype=np.float64)
    new_size = np.maximum(1, np.round(original_size * original_spacing / new_spacing)).astype(int)

    resampler = sitk.ResampleImageFilter()
    resampler.SetOutputSpacing(target_spacing)
    resampler.SetSize([int(value) for value in new_size])
    resampler.SetOutputOrigin(image.GetOrigin())
    resampler.SetOutputDirection(image.GetDirection())
    resampler.SetDefaultPixelValue(0)
    resampler.SetInterpolator(
        sitk.sitkNearestNeighbor
        if interpolation == "nearest"
        else sitk.sitkBSpline
        if interpolation == "bspline"
        else sitk.sitkLinear
    )
    out_zyx = sitk.GetArrayFromImage(resampler.Execute(image))
    out_xyz = np.transpose(out_zyx, (2, 1, 0))
    vtk_type = vtk.VTK_UNSIGNED_CHAR if interpolation == "nearest" else vtk_image.GetScalarType()
    out = _vtk_image_from_array(
        out_xyz,
        vtk_image,
        origin=vtk_image.GetOrigin(),
        vtk_array_type=vtk_type,
    )
    out.SetSpacing(*target_spacing)
    return out


def crop_vtk_images_to_bbox_ratio(
    vtk_images,
    vtk_mask,
    *,
    bbox_ratio=DEFAULT_FEMUR_BBOX_RATIO,
    bbox_crop_from=None,
    labels=None,
):
    """Crop VTK images to a mask-bbox ratio and return the retained crop face.

    Ratios are interpreted in Ogo's x/y/z VTK image order. The returned
    crop-face image has value 1 on the newly exposed one-sided crop surface.
    For the hip recipe this is the z-min face created by
    ``bbox_ratio=(1, 1.2, None)`` and ``bbox_crop_from=(None, "min", None)``.
    """
    import numpy as np
    import vtk

    from ogo.util.vtk_image import vtk_image_to_numpy

    images = list(vtk_images)
    mask_data = vtk_image_to_numpy(vtk_mask)
    if labels:
        active = np.isin(mask_data, sorted(int(label) for label in labels))
    else:
        active = mask_data != 0
    if not np.any(active):
        raise ValueError("Cannot apply bbox-ratio crop to an empty femur mask.")

    ratio_xyz = _parse_bbox_ratio_recipe(bbox_ratio)
    crop_from_xyz = _parse_bbox_crop_from_recipe(bbox_crop_from)
    numeric_axes = [axis for axis, value in enumerate(ratio_xyz) if value is not None]
    if not numeric_axes:
        empty_face = np.zeros(mask_data.shape, dtype=np.uint8)
        return images, _vtk_image_from_array(empty_face, vtk_mask, origin=vtk_mask.GetOrigin(), vtk_array_type=vtk.VTK_UNSIGNED_CHAR), {
            "enabled": False,
            "ratio_recipe": tuple(bbox_ratio),
            "ratio_xyz": ratio_xyz,
            "crop_from_xyz": crop_from_xyz,
        }

    reference_axes = [axis for axis in numeric_axes if np.isclose(float(ratio_xyz[axis]), 1.0)]
    if not reference_axes:
        raise ValueError("bbox_ratio must contain one preserved axis with value 1.")

    coords = np.argwhere(active)
    lo = coords.min(axis=0).astype(np.int64)
    hi = (coords.max(axis=0) + 1).astype(np.int64)
    size = hi - lo
    spacing = np.asarray(vtk_mask.GetSpacing(), dtype=np.float64)
    physical_size = size.astype(np.float64) * spacing
    reference_axis = min(reference_axes, key=lambda axis: float(physical_size[axis]))
    reference_length_mm = float(size[reference_axis]) * float(spacing[reference_axis])

    out_lo = lo.copy()
    out_hi = hi.copy()
    crop_surface = None
    for axis, axis_ratio in enumerate(ratio_xyz):
        if axis_ratio is None:
            continue
        target_mm = reference_length_mm * float(axis_ratio)
        target_voxels = max(1, int(round(target_mm / float(spacing[axis]))))
        target_voxels = min(int(size[axis]), target_voxels)
        mode = crop_from_xyz[axis]
        if mode == "min":
            start = int(hi[axis]) - target_voxels
            crop_surface = {"axis": axis, "side": "min", "local_index": 0, "normal_sign": -1.0}
        elif mode == "max":
            start = int(lo[axis])
            crop_surface = {
                "axis": axis,
                "side": "max",
                "local_index": target_voxels - 1,
                "normal_sign": 1.0,
            }
        else:
            center = 0.5 * (float(lo[axis]) + float(hi[axis]))
            start = int(round(center - 0.5 * float(target_voxels)))
        start = max(int(lo[axis]), min(start, int(hi[axis]) - target_voxels))
        out_lo[axis] = start
        out_hi[axis] = start + target_voxels

    slices = tuple(slice(int(out_lo[axis]), int(out_hi[axis])) for axis in range(3))
    origin = tuple(
        float(vtk_mask.GetOrigin()[axis]) + float(out_lo[axis]) * float(spacing[axis])
        for axis in range(3)
    )
    cropped_images = [
        _vtk_image_from_array(vtk_image_to_numpy(image)[slices], image, origin=origin)
        for image in images
    ]

    cropped_active = active[slices]
    crop_face = np.zeros(cropped_active.shape, dtype=np.uint8)
    if crop_surface is not None:
        axis = int(crop_surface["axis"])
        index = int(crop_surface["local_index"])
        face_slice = [slice(None), slice(None), slice(None)]
        face_slice[axis] = index
        crop_face[tuple(face_slice)] = cropped_active[tuple(face_slice)].astype(np.uint8)

    crop_face_image = _vtk_image_from_array(
        crop_face,
        vtk_mask,
        origin=origin,
        vtk_array_type=vtk.VTK_UNSIGNED_CHAR,
    )
    meta = {
        "enabled": True,
        "ratio_recipe": tuple(bbox_ratio),
        "crop_from_recipe": tuple(bbox_crop_from) if bbox_crop_from is not None else None,
        "ratio_xyz": ratio_xyz,
        "crop_from_xyz": crop_from_xyz,
        "reference_axis": "xyz"[reference_axis],
        "reference_length_mm": reference_length_mm,
        "input_bbox_xyz": tuple((int(lo[axis]), int(hi[axis])) for axis in range(3)),
        "crop_slices_xyz": tuple((int(out_lo[axis]), int(out_hi[axis])) for axis in range(3)),
        "output_shape_xyz": tuple(int(value) for value in cropped_active.shape),
        "output_origin": origin,
        "crop_surface": None
        if crop_surface is None
        else {
            "axis": "xyz"[int(crop_surface["axis"])],
            "side": crop_surface["side"],
            "normal_sign": crop_surface["normal_sign"],
        },
        "crop_face_voxels": int(crop_face.sum()),
    }
    return cropped_images, crop_face_image, meta


def crop_face_support_vector(crop_face_vtk, mask_vtk):
    """Estimate the outward normal of a transformed crop-face label image."""
    import numpy as np

    from ogo.util.vtk_image import vtk_image_to_numpy

    face = vtk_image_to_numpy(crop_face_vtk) != 0
    bone = vtk_image_to_numpy(mask_vtk) != 0
    if not np.any(face):
        raise ValueError("Cannot estimate crop-face support vector from an empty face mask.")
    if not np.any(bone):
        raise ValueError("Cannot orient crop-face support vector from an empty femur mask.")

    spacing = np.asarray(crop_face_vtk.GetSpacing(), dtype=np.float64)
    origin = np.asarray(crop_face_vtk.GetOrigin(), dtype=np.float64)
    face_coords = np.argwhere(face).astype(np.float64)
    bone_coords = np.argwhere(bone).astype(np.float64)
    face_points = origin + face_coords * spacing
    bone_points = origin + bone_coords * spacing
    if face_points.shape[0] < 3:
        raise ValueError("At least three crop-face voxels are required to estimate a surface normal.")

    centered = face_points - face_points.mean(axis=0)
    _u, _s, vh = np.linalg.svd(centered, full_matrices=False)
    normal = np.asarray(vh[-1], dtype=np.float64)
    norm = float(np.linalg.norm(normal))
    if norm <= 0.0:
        raise ValueError("Cannot estimate crop-face support vector from degenerate face points.")
    normal /= norm

    # The vector should point out of the retained bone, toward the cropped-away
    # side. Orient it by comparing the face centroid to the retained-bone
    # centroid.
    into_bone = bone_points.mean(axis=0) - face_points.mean(axis=0)
    if float(np.dot(normal, into_bone)) > 0.0:
        normal = -normal
    return tuple(float(value) for value in normal)


def _image_index_points(vtk_image, indices):
    """Return physical center coordinates for VTK image array indices."""
    import numpy as np

    extent = vtk_image.GetExtent()
    offset = np.asarray([extent[0], extent[2], extent[4]], dtype=np.float64)
    spacing = np.asarray(vtk_image.GetSpacing(), dtype=np.float64)
    origin = np.asarray(vtk_image.GetOrigin(), dtype=np.float64)
    return origin + (np.asarray(indices, dtype=np.float64) + offset) * spacing


def _unit_vector(vector, name):
    import numpy as np

    values = np.asarray(vector, dtype=np.float64)
    if values.shape != (3,) or not np.all(np.isfinite(values)):
        raise ValueError(f"{name} must be a finite 3-vector.")
    norm = float(np.linalg.norm(values))
    if norm <= 0.0:
        raise ValueError(f"{name} must be non-zero.")
    return values / norm


def bbox_relative_oriented_contact_plane(
    model_bounds,
    *,
    center_fraction,
    size_fraction,
    normal,
    width_axis=(-1.0, 0.0, 0.0),
    shape="anatomy",
):
    """Return a bbox-relative plane that keeps a measured contact normal.

    The center and footprint scale with the generated model bounding box. The
    normal comes from measured geometry, such as a transformed crop face, so the
    plane can follow an oblique distal shaft cut while still using a stable
    recipe-sized contact patch.
    """
    import numpy as np

    if len(model_bounds) != 6:
        raise ValueError("model_bounds must contain x/y/z min/max values.")
    if len(center_fraction) != 3:
        raise ValueError("center_fraction must contain x/y/z fractions.")
    if len(size_fraction) != 2:
        raise ValueError("size_fraction must contain the two in-plane scale factors.")

    normal_vec = _unit_vector(normal, "contact plane normal")
    width_vec = np.asarray(width_axis, dtype=np.float64)
    if width_vec.shape != (3,) or not np.all(np.isfinite(width_vec)):
        raise ValueError("width_axis must be a finite 3-vector.")

    width_vec = width_vec - normal_vec * float(np.dot(width_vec, normal_vec))
    if float(np.linalg.norm(width_vec)) <= 1.0e-8:
        fallback_axes = np.eye(3, dtype=np.float64)
        width_vec = max(
            (axis - normal_vec * float(np.dot(axis, normal_vec)) for axis in fallback_axes),
            key=lambda candidate: float(np.linalg.norm(candidate)),
        )
    v_axis = _unit_vector(width_vec, "contact plane width axis")
    u_axis = _unit_vector(np.cross(v_axis, normal_vec), "contact plane long axis")

    center = []
    for axis in range(3):
        lo, hi = _axis_bounds(model_bounds, axis)
        center.append(lo + float(center_fraction[axis]) * (hi - lo))

    size = (
        _axis_extent_from_vector(model_bounds, u_axis) * float(size_fraction[0]),
        _axis_extent_from_vector(model_bounds, v_axis) * float(size_fraction[1]),
    )
    if size[0] <= 0.0 or size[1] <= 0.0:
        raise ValueError("size_fraction values must produce positive plane dimensions.")

    return {
        "center": tuple(float(value) for value in center),
        "normal": tuple(float(value) for value in normal_vec),
        "outward_normal": tuple(float(value) for value in -normal_vec),
        "u_axis": tuple(float(value) for value in u_axis),
        "v_axis": tuple(float(value) for value in v_axis),
        "size": tuple(float(value) for value in size),
        "shape": str(shape).strip().lower(),
    }


def _inside_plane_shape(shape, u_values, v_values, half_u, half_v, *, tolerance):
    import numpy as np

    token = str(shape).strip().lower()
    if token == "square":
        half = min(float(half_u), float(half_v))
        return (np.abs(u_values) <= half + tolerance) & (np.abs(v_values) <= half + tolerance)
    if token in {"round", "circle", "circular", "oval"}:
        safe_u = max(float(half_u) + tolerance, 1.0e-9)
        safe_v = max(float(half_v) + tolerance, 1.0e-9)
        return ((u_values / safe_u) ** 2 + (v_values / safe_v) ** 2) <= 1.0
    return (np.abs(u_values) <= float(half_u) + tolerance) & (
        np.abs(v_values) <= float(half_v) + tolerance
    )


def _plane_bucket_key(u_value, v_value, *, spacing):
    import math

    resolution = max(min(float(value) for value in spacing), 1.0e-6)
    return (
        int(math.floor(float(u_value) / resolution + 0.5)),
        int(math.floor(float(v_value) / resolution + 0.5)),
    )


def crop_face_contact_plane(crop_face_vtk, mask_vtk):
    """Fit the transformed distal crop face as a contact plane.

    The returned ``normal`` points from the crop face into the retained femur.
    ``outward_normal`` points toward the cropped-away shaft. Both are derived
    from the transformed crop-face label, so the shaft support angle follows the
    same registration transform as the femur.
    """
    import numpy as np

    from ogo.util.vtk_image import vtk_image_to_numpy

    face = vtk_image_to_numpy(crop_face_vtk) != 0
    bone = vtk_image_to_numpy(mask_vtk) != 0
    if not np.any(face):
        raise ValueError("Cannot fit a contact plane from an empty crop-face mask.")
    if not np.any(bone):
        raise ValueError("Cannot orient a crop-face contact plane from an empty femur mask.")

    face_indices = np.argwhere(face)
    if face_indices.shape[0] < 3:
        raise ValueError("At least three crop-face voxels are required to fit a contact plane.")
    face_points = _image_index_points(crop_face_vtk, face_indices)
    bone_points = _image_index_points(mask_vtk, np.argwhere(bone))

    center = face_points.mean(axis=0)
    centered = face_points - center
    _u_svd, _s_svd, vh = np.linalg.svd(centered, full_matrices=False)
    u_axis = _unit_vector(vh[0], "crop-face u axis")
    v_axis = _unit_vector(vh[1], "crop-face v axis")
    normal = _unit_vector(vh[-1], "crop-face normal")

    into_bone = bone_points.mean(axis=0) - center
    if float(np.dot(normal, into_bone)) < 0.0:
        normal = -normal
    if float(np.dot(np.cross(u_axis, v_axis), normal)) < 0.0:
        v_axis = -v_axis

    u_values = centered @ u_axis
    v_values = centered @ v_axis
    min_spacing = min(float(value) for value in crop_face_vtk.GetSpacing())
    size = (
        max(float(np.ptp(u_values)) + min_spacing, min_spacing),
        max(float(np.ptp(v_values)) + min_spacing, min_spacing),
    )
    return {
        "center": tuple(float(value) for value in center),
        "normal": tuple(float(value) for value in normal),
        "outward_normal": tuple(float(value) for value in -normal),
        "u_axis": tuple(float(value) for value in u_axis),
        "v_axis": tuple(float(value) for value in v_axis),
        "size": tuple(float(value) for value in size),
        "shape": "anatomy",
    }


def projected_crop_face_surface_vtk(material_vtk, plane, *, intrusion=0.0, output_value=1):
    """Project a fitted crop-face plane onto the first active femur surface."""
    import numpy as np
    import vtk

    from ogo.util.vtk_image import numpy_to_vtk_image, vtk_image_to_numpy

    active = vtk_image_to_numpy(material_vtk) != 0
    out = np.zeros(active.shape, dtype=np.uint8)
    if not np.any(active):
        return numpy_to_vtk_image(out, material_vtk, vtk_array_type=vtk.VTK_UNSIGNED_CHAR)

    center = np.asarray(plane["center"], dtype=np.float64)
    normal = _unit_vector(plane["normal"], "plane normal")
    u_axis = _unit_vector(plane["u_axis"], "plane u axis")
    v_axis = _unit_vector(plane["v_axis"], "plane v axis")
    size = plane.get("size", (0.0, 0.0))
    half_u = max(float(size[0]) / 2.0, 0.5)
    half_v = max(float(size[1]) / 2.0, 0.5)
    spacing = tuple(float(value) for value in material_vtk.GetSpacing())
    tolerance = max(min(spacing) * 0.75, 1.0e-6)

    indices = np.argwhere(active)
    points = _image_index_points(material_vtk, indices)
    relative = points - center
    distances = relative @ normal
    u_values = relative @ u_axis
    v_values = relative @ v_axis
    inside = _inside_plane_shape(
        plane.get("shape", "anatomy"),
        u_values,
        v_values,
        half_u,
        half_v,
        tolerance=tolerance,
    )
    candidates = inside & (distances >= -tolerance)
    if not np.any(candidates):
        return numpy_to_vtk_image(out, material_vtk, vtk_array_type=vtk.VTK_UNSIGNED_CHAR)

    best = {}
    for voxel, distance, u_value, v_value in zip(
        indices[candidates],
        distances[candidates],
        u_values[candidates],
        v_values[candidates],
        strict=True,
    ):
        if float(distance) < -tolerance:
            continue
        key = _plane_bucket_key(u_value, v_value, spacing=spacing)
        current = best.get(key)
        if current is None or float(distance) < current[0]:
            best[key] = (float(distance), np.asarray(voxel, dtype=np.int64))
    if not best:
        return numpy_to_vtk_image(out, material_vtk, vtk_array_type=vtk.VTK_UNSIGNED_CHAR)

    first_distance = min(value[0] for value in best.values())
    max_distance = first_distance + max(float(intrusion), 0.0) + tolerance
    kept = [voxel for distance, voxel in best.values() if distance <= max_distance]
    if kept:
        kept_indices = np.stack(kept, axis=0)
        out[tuple(kept_indices.T)] = np.uint8(output_value)
    return numpy_to_vtk_image(out, material_vtk, vtk_array_type=vtk.VTK_UNSIGNED_CHAR)


def mirror_polydata_x(polydata):
    """Return a left/right mirrored copy of polydata around its x-bounds center."""
    import vtk

    bounds = polydata.GetBounds()
    center_x = (bounds[0] + bounds[1]) / 2.0

    transform = vtk.vtkTransform()
    transform.Translate(center_x, 0.0, 0.0)
    transform.Scale(-1.0, 1.0, 1.0)
    transform.Translate(-center_x, 0.0, 0.0)

    transform_filter = vtk.vtkTransformPolyDataFilter()
    transform_filter.SetInputData(polydata)
    transform_filter.SetTransform(transform)
    transform_filter.Update()

    reverse = vtk.vtkReverseSense()
    reverse.SetInputData(transform_filter.GetOutput())
    reverse.ReverseCellsOn()
    reverse.ReverseNormalsOn()
    reverse.Update()

    mirrored = vtk.vtkPolyData()
    mirrored.DeepCopy(reverse.GetOutput())
    return mirrored


def tilted_side_support_vector(angle_degrees=-20):
    """Return the distal support vector used for the sideways-fall fixture."""
    import numpy as np

    theta = np.radians(angle_degrees)
    rotation_x = np.array(
        [
            [1, 0, 0],
            [0, np.cos(theta), -np.sin(theta)],
            [0, np.sin(theta), np.cos(theta)],
        ]
    )
    return rotation_x @ np.array([0, 0, -1])


def cortical_compartment_mask(compartment_vtk, *, cortical_label=1, trabecular_label=2):
    """Return a binary VTK mask for cortical voxels from a trab/cort label image."""
    import numpy as np
    import vtk

    from ogo.util.vtk_image import numpy_to_vtk_image, vtk_image_to_numpy

    data = np.rint(vtk_image_to_numpy(compartment_vtk)).astype(np.int32)
    cortical_label = int(cortical_label)
    trabecular_label = int(trabecular_label)
    present = set(int(value) for value in np.unique(data) if int(value) != 0)
    missing = [label for label in (cortical_label, trabecular_label) if label not in present]
    if missing:
        raise ValueError(
            "Compartment mask is missing required label(s): {}. "
            "Expected cortical={} and trabecular={}.".format(
                ", ".join(str(label) for label in missing),
                cortical_label,
                trabecular_label,
            )
        )
    return numpy_to_vtk_image(
        (data == cortical_label).astype(np.uint8),
        compartment_vtk,
        vtk_array_type=vtk.VTK_UNSIGNED_CHAR,
    )


def distal_cut_z_from_reference(reference_bounds, retained_length_mm):
    """Return a flat distal cut plane in aligned reference coordinates."""
    retained_length_mm = float(retained_length_mm)
    if retained_length_mm <= 0:
        raise ValueError("retained_length_mm must be positive.")
    return float(reference_bounds[5]) - retained_length_mm


def femur_z_coverage(vtk_mask):
    """Return physical z coverage for nonzero voxels in a VTK femur mask."""
    import numpy as np
    from ogo.util.vtk_image import vtk_image_to_numpy

    dims = vtk_mask.GetDimensions()
    extent = vtk_mask.GetExtent()
    spacing = vtk_mask.GetSpacing()
    origin = vtk_mask.GetOrigin()
    data = vtk_image_to_numpy(vtk_mask)
    occupied = np.where(np.any(data != 0, axis=(0, 1)))[0]
    if occupied.size == 0:
        raise ValueError("Cannot standardize femur shaft length from an empty mask.")
    z = origin[2] + (occupied + extent[4]) * spacing[2]
    return float(z.min()), float(z.max())


def femur_z_profile(vtk_mask):
    """Return per-z cross-section profile arrays for a VTK femur mask."""
    import numpy as np
    from ogo.util.vtk_image import vtk_image_to_numpy

    dims = vtk_mask.GetDimensions()
    extent = vtk_mask.GetExtent()
    spacing = vtk_mask.GetSpacing()
    origin = vtk_mask.GetOrigin()
    data = vtk_image_to_numpy(vtk_mask) != 0
    z_coords = origin[2] + (np.arange(dims[2]) + extent[4]) * spacing[2]
    y_coords = origin[1] + (np.arange(dims[1]) + extent[2]) * spacing[1]

    rows = []
    for z_index, z_value in enumerate(z_coords):
        plane = data[:, :, z_index]
        if not plane.any():
            continue
        y_indices = np.where(plane)[1]
        rows.append(
            (
                float(z_value),
                float(plane.sum()),
                float(y_coords[y_indices].min()),
                float(y_coords[y_indices].max()),
            )
        )
    if not rows:
        raise ValueError("Cannot compute femur profile from an empty mask.")

    profile = np.asarray(rows, dtype=float)
    return {
        "z": profile[:, 0],
        "area": profile[:, 1],
        "y_min": profile[:, 2],
        "y_max": profile[:, 3],
        "width_y": profile[:, 3] - profile[:, 2],
    }


def femoral_head_side_z_bounds(vtk_mask, *, y_fraction=0.2):
    """Estimate the z extent of the femoral-head-side y band in an aligned femur mask."""
    import numpy as np

    from ogo.util.vtk_image import vtk_image_to_numpy

    data = vtk_image_to_numpy(vtk_mask) != 0
    coords = np.array(np.where(data))
    if coords.size == 0:
        raise ValueError("Cannot estimate femoral-head-side bounds from an empty mask.")

    dims = vtk_mask.GetDimensions()
    extent = vtk_mask.GetExtent()
    spacing = vtk_mask.GetSpacing()
    origin = vtk_mask.GetOrigin()
    y = origin[1] + (np.arange(dims[1]) + extent[2]) * spacing[1]
    z = origin[2] + (np.arange(dims[2]) + extent[4]) * spacing[2]

    y_min = float(y[coords[1].min()])
    y_max = float(y[coords[1].max()])
    y_limit = y_min + (y_max - y_min) * float(y_fraction)
    head_side = data & (y[None, :, None] >= y_min) & (y[None, :, None] <= y_limit)
    head_coords = np.array(np.where(head_side))
    if head_coords.size == 0:
        raise ValueError("Cannot estimate femoral-head-side bounds from an empty y band.")

    return {
        "y_min": y_min,
        "y_max": y_limit,
        "z_min": float(z[head_coords[2].min()]),
        "z_max": float(z[head_coords[2].max()]),
    }


def pad_vtk_images_to_foreground_margin(
    vtk_images,
    vtk_mask,
    *,
    margin_mm=DEFAULT_FEMUR_INPUT_MARGIN_MM,
    constants=None,
):
    """Pad images so femur foreground has a physical margin to every extent face."""
    import math
    import numpy as np
    import vtk

    from ogo.util.vtk_image import vtk_image_to_numpy

    images = list(vtk_images)
    constants = [0] * len(images) if constants is None else list(constants)
    if len(constants) != len(images):
        raise ValueError("constants must match vtk_images length.")

    margin_mm = max(0.0, float(margin_mm))
    if margin_mm == 0.0:
        return images, {"lower": (0, 0, 0), "upper": (0, 0, 0), "margin_mm": 0.0}

    mask = vtk_image_to_numpy(vtk_mask) != 0
    coords = np.array(np.where(mask))
    if coords.size == 0:
        raise ValueError("Cannot pad femur images from an empty mask.")

    dims = np.array(mask.shape, dtype=int)
    mins = coords.min(axis=1)
    maxs = coords.max(axis=1)
    spacing = vtk_mask.GetSpacing()
    margin_voxels = np.array(
        [max(0, int(math.ceil(margin_mm / max(float(value), 1.0e-6)))) for value in spacing],
        dtype=int,
    )
    lower = np.maximum(0, margin_voxels - mins)
    upper = np.maximum(0, margin_voxels - ((dims - 1) - maxs))
    if not np.any(lower) and not np.any(upper):
        return images, {
            "lower": tuple(int(v) for v in lower),
            "upper": tuple(int(v) for v in upper),
            "margin_mm": margin_mm,
            "margin_voxels": tuple(int(v) for v in margin_voxels),
        }

    extent = vtk_mask.GetExtent()
    output_extent = (
        int(extent[0] - lower[0]),
        int(extent[1] + upper[0]),
        int(extent[2] - lower[1]),
        int(extent[3] + upper[1]),
        int(extent[4] - lower[2]),
        int(extent[5] + upper[2]),
    )
    padded = []
    for image, constant in zip(images, constants):
        pad = vtk.vtkImageConstantPad()
        pad.SetInputData(image)
        pad.SetOutputWholeExtent(output_extent)
        pad.SetConstant(float(constant))
        pad.Update()
        out = vtk.vtkImageData()
        out.DeepCopy(pad.GetOutput())
        dims_out = out.GetDimensions()
        spacing_out = out.GetSpacing()
        origin = image.GetOrigin()
        out.SetOrigin(
            float(origin[0]) + float(output_extent[0]) * float(spacing_out[0]),
            float(origin[1]) + float(output_extent[2]) * float(spacing_out[1]),
            float(origin[2]) + float(output_extent[4]) * float(spacing_out[2]),
        )
        out.SetExtent(0, dims_out[0] - 1, 0, dims_out[1] - 1, 0, dims_out[2] - 1)
        padded.append(out)

    return padded, {
        "lower": tuple(int(v) for v in lower),
        "upper": tuple(int(v) for v in upper),
        "margin_mm": margin_mm,
        "margin_voxels": tuple(int(v) for v in margin_voxels),
    }


def _smooth_profile(values, window=7):
    """Smooth one profile with an odd moving-average window."""
    import numpy as np

    values = np.asarray(values, dtype=float)
    window = max(1, int(window))
    if window % 2 == 0:
        window += 1
    if values.size < window or window == 1:
        return values.copy()
    pad = window // 2
    padded = np.pad(values, (pad, pad), mode="edge")
    return np.convolve(padded, np.ones(window, dtype=float) / float(window), mode="valid")


def _peak_center_z(z, values, indices, *, relative_height=0.95):
    """Return the center of the high plateau around a dominant local peak."""
    import numpy as np

    peak_index = int(indices[np.argmax(values[indices])])
    threshold = float(values[peak_index]) * float(relative_height)
    left = peak_index
    right = peak_index
    valid = set(int(i) for i in indices)
    while left - 1 in valid and values[left - 1] >= threshold:
        left -= 1
    while right + 1 in valid and values[right + 1] >= threshold:
        right += 1
    plateau = np.arange(left, right + 1, dtype=int)
    weights = np.maximum(values[plateau], 0.0)
    if float(weights.sum()) <= 0.0:
        return peak_index, float(z[peak_index])
    return peak_index, float(np.average(z[plateau], weights=weights))


def detect_lesser_trochanter_cut_z(
    vtk_mask,
    *,
    distal_offset_mm=0.0,
    distal_offset_percent=None,
    min_distal_to_greater_mm=8.0,
    max_distal_to_greater_mm=45.0,
    min_distal_coverage_mm=0.0,
    distal_fov_tolerance_mm=2.0,
):
    """Detect a lesser-trochanter-based flat distal cut plane.

    The aligned sideways-fall reference frame has femur length along z and the
    lateral greater-trochanter side at high y. The greater trochanter is found
    from the smoothed high-y profile. The lesser-trochanter/transition level is
    then the dominant smoothed cross-sectional area peak distal to that point.
    """
    import numpy as np

    profile = femur_z_profile(vtk_mask)
    z = profile["z"]
    area = _smooth_profile(profile["area"], window=7)
    y_max = _smooth_profile(profile["y_max"], window=7)

    z_min = float(z.min())
    z_max = float(z.max())
    if z_max - z_min < max_distal_to_greater_mm:
        raise ValueError(
            "Femur scan does not include enough proximal-distal coverage to identify "
            "the lesser trochanter."
        )

    proximal_mask = z >= (z_min + 0.55 * (z_max - z_min))
    if not np.any(proximal_mask):
        raise ValueError("Cannot identify greater trochanter from femur profile.")
    proximal_indices = np.where(proximal_mask)[0]
    greater_index = int(proximal_indices[np.argmax(y_max[proximal_indices])])
    greater_z = float(z[greater_index])

    distal_mask = (
        (z <= greater_z - float(min_distal_to_greater_mm))
        & (z >= greater_z - float(max_distal_to_greater_mm))
    )
    distal_indices = np.where(distal_mask)[0]
    if distal_indices.size < 5:
        raise ValueError(
            "Femur scan does not include the distal profile needed to identify "
            "the lesser trochanter."
        )

    lesser_index, lesser_z = _peak_center_z(z, area, distal_indices)
    if distal_offset_percent is not None:
        distal_offset_mm = abs(greater_z - lesser_z) * float(distal_offset_percent) / 100.0
    cut_z = lesser_z - float(distal_offset_mm)

    required_z_min = cut_z - float(min_distal_coverage_mm)
    if z_min > required_z_min:
        missing = z_min - required_z_min
        if missing <= float(distal_fov_tolerance_mm):
            cut_z = z_min + float(min_distal_coverage_mm)
        else:
            available = max(0.0, cut_z - z_min)
            raise ValueError(
                "Femur scan does not extend distal to the lesser-trochanter cut plane "
                f"(required {float(min_distal_coverage_mm):.1f} mm, available {available:.1f} mm)."
            )

    shaft_area = float(np.median(area[z < lesser_z - 10.0])) if np.any(z < lesser_z - 10.0) else 0.0
    if shaft_area > 0 and float(area[lesser_index]) < 1.08 * shaft_area:
        raise ValueError(
            "Could not identify a clear lesser-trochanter cross-section peak; "
            "the femur field of view is likely incomplete or the alignment failed."
        )

    return {
        "cut_z": float(cut_z),
        "lesser_trochanter_z": lesser_z,
        "greater_trochanter_z": greater_z,
        "mask_z_min": z_min,
        "mask_z_max": z_max,
        "retained_length_mm": float(z_max - cut_z),
        "distal_offset_mm": float(distal_offset_mm),
        "distal_offset_percent": None if distal_offset_percent is None else float(distal_offset_percent),
        "lesser_area": float(area[lesser_index]),
        "shaft_area_median": shaft_area,
    }


def flat_crop_vtk_image_below_z(vtk_image, cut_z):
    """Zero voxels below a physical z cut plane while preserving image metadata."""
    import numpy as np
    from ogo.util.vtk_image import numpy_to_vtk_image, vtk_image_to_numpy

    dims = vtk_image.GetDimensions()
    extent = vtk_image.GetExtent()
    spacing = vtk_image.GetSpacing()
    origin = vtk_image.GetOrigin()
    data = vtk_image_to_numpy(vtk_image, copy=True)
    z = origin[2] + (np.arange(dims[2]) + extent[4]) * spacing[2]
    data[:, :, z < float(cut_z)] = 0

    return numpy_to_vtk_image(data, vtk_image)


def standardize_femur_shaft_length(
    image_vtk,
    mask_vtk,
    *,
    reference_bounds,
    retained_length_mm=None,
    cut_mode="lesser_trochanter",
    lesser_trochanter_distal_offset_mm=50.0,
    lesser_trochanter_distal_offset_percent=None,
):
    """Apply a flat distal crop to aligned femur image and mask.

    In the default mode, the cut plane is detected from the lesser-trochanter
    cross-section profile. A fixed retained length can be requested explicitly
    for controlled debugging runs.
    """
    cut_mode = str(cut_mode).lower()
    if cut_mode == "lesser_trochanter":
        meta = detect_lesser_trochanter_cut_z(
            mask_vtk,
            distal_offset_mm=lesser_trochanter_distal_offset_mm,
            distal_offset_percent=lesser_trochanter_distal_offset_percent,
        )
        cut_z = meta["cut_z"]
        retained_length_mm = meta["retained_length_mm"]
    elif cut_mode == "fixed_length":
        if retained_length_mm is None:
            raise ValueError("retained_length_mm is required for fixed_length femur cuts.")
        cut_z = distal_cut_z_from_reference(reference_bounds, retained_length_mm)
        meta = {
            "cut_z": cut_z,
            "retained_length_mm": float(retained_length_mm),
        }
    else:
        raise ValueError("cut_mode must be 'lesser_trochanter' or 'fixed_length'.")

    z_min, z_max = femur_z_coverage(mask_vtk)
    warnings = []
    try:
        head_side_bounds = femoral_head_side_z_bounds(mask_vtk)
    except ValueError:
        head_side_bounds = None
    if head_side_bounds is not None:
        head_side_z_min = float(head_side_bounds["z_min"])
        if float(cut_z) > head_side_z_min + 2.0:
            warnings.append(
                "distal cut plane appears to enter the femoral-head-side region "
                f"(cut_z={float(cut_z):.3f}, head_side_z_min={head_side_z_min:.3f})"
            )
        meta["femoral_head_side_z_min"] = head_side_z_min
        meta["femoral_head_side_z_max"] = float(head_side_bounds["z_max"])
    if retained_length_mm is not None and float(retained_length_mm) < 40.0:
        warnings.append(
            "retained proximal femur length is unusually short "
            f"({float(retained_length_mm):.3f} mm)"
        )
    if z_min > cut_z:
        available = max(0.0, float(reference_bounds[5]) - z_min)
        raise ValueError(
            "Femur scan does not include enough distal shaft for the standard crop "
            f"(required retained length {float(retained_length_mm):.1f} mm, "
            f"available approximately {available:.1f} mm)."
        )
    meta.update(
        {
            "cut_z": cut_z,
            "mask_z_min": z_min,
            "mask_z_max": z_max,
            "retained_length_mm": float(retained_length_mm),
            "cut_mode": cut_mode,
            "warnings": warnings,
        }
    )
    return (
        flat_crop_vtk_image_below_z(image_vtk, cut_z),
        flat_crop_vtk_image_below_z(mask_vtk, cut_z),
        meta,
    )


def swap_xz_footprint(bounds):
    """Swap x/z footprint dimensions of physical bounds around the same center."""
    x_center = (float(bounds[0]) + float(bounds[1])) / 2.0
    z_center = (float(bounds[4]) + float(bounds[5])) / 2.0
    x_length = float(bounds[1]) - float(bounds[0])
    z_length = float(bounds[5]) - float(bounds[4])
    return (
        x_center - z_length / 2.0,
        x_center + z_length / 2.0,
        float(bounds[2]),
        float(bounds[3]),
        z_center - x_length / 2.0,
        z_center + x_length / 2.0,
    )


def expand_z_footprint(bounds, extension_mm):
    """Expand physical bounds along z around the same center."""
    z_center = (float(bounds[4]) + float(bounds[5])) / 2.0
    z_length = float(bounds[5]) - float(bounds[4]) + float(extension_mm)
    if z_length <= 0:
        raise ValueError("Expanded z footprint length must be positive.")
    return (
        float(bounds[0]),
        float(bounds[1]),
        float(bounds[2]),
        float(bounds[3]),
        z_center - z_length / 2.0,
        z_center + z_length / 2.0,
    )


def expand_xz_footprint(bounds, *, x_extension_mm=0.0, z_extension_mm=0.0):
    """Expand physical bounds along x and z around the same center."""
    x_center = (float(bounds[0]) + float(bounds[1])) / 2.0
    z_center = (float(bounds[4]) + float(bounds[5])) / 2.0
    x_length = float(bounds[1]) - float(bounds[0]) + float(x_extension_mm)
    z_length = float(bounds[5]) - float(bounds[4]) + float(z_extension_mm)
    if x_length <= 0 or z_length <= 0:
        raise ValueError("Expanded x/z footprint lengths must be positive.")
    return (
        x_center - x_length / 2.0,
        x_center + x_length / 2.0,
        float(bounds[2]),
        float(bounds[3]),
        z_center - z_length / 2.0,
        z_center + z_length / 2.0,
    )
