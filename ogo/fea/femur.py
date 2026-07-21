"""Femur-specific defaults and helpers for sideways-fall FE generation.

Shared point-cloud and boundary-contact primitives live in ``alignment`` and
``boundary``. They are re-exported here so older femur workflow imports keep
working while this module stays focused on femur defaults and femur-only
geometry.

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
4. Stabilize registration with a fixed 100 mm proximal rough crop after
   isotropic resampling, then standardize the distal shaft after ICP with a flat
   aligned-frame crop that keeps z length at 1.2 times the aligned y width.
   The post-ICP crop face becomes the distal support surface. ``bbox_ratio``,
   ``proximal_box_ratio``, ``lesser_trochanter``, and ``fixed_length`` modes are
   available for debugging or historical comparisons.
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

from ogo.fea.alignment import (
    estimate_rigid_icp,
    point_cloud_axis_lengths,
    polydata_from_points,
    polydata_points,
    sample_points,
    surface_points_from_vtk_mask,
)
from ogo.fea.boundary import (
    _axis_bounds,
    _axis_extent_from_vector,
    bbox_relative_contact_bounds as bbox_relative_fixture_bounds,
    bbox_relative_contact_direction as bbox_relative_fixture_direction,
    bbox_relative_contact_plane as bbox_relative_fixture_plane,
    foreground_voxel_center_bounds,
    foreground_voxel_center_bounds_from_mask,
)


LEFT_FEMUR = 1
RIGHT_FEMUR = 2

DEFAULT_FEMUR_ISO_RESOLUTION_MM = 1.0
DEFAULT_FEMUR_FE_DISPLACEMENT = -4.0
DEFAULT_FEMUR_TARGET_DISPLACEMENT_PERCENT = 4.0
DEFAULT_FEMUR_MASK_SMOOTHING_SPACING_THRESHOLD_MM = 2.0
FEMORAL_HEAD_FIXTURE_WIDTH_EXTENSION_MM = 10.0
FEMORAL_HEAD_FIXTURE_LONG_AXIS_EXTENSION_MM = 80.0
SIDEWAYS_FALL_FIXTURE_SHAPE = "anatomy"
SIDEWAYS_FALL_FIXTURE_SIZE_FRACTION = (1.0, 1.0)
FEMORAL_HEAD_FIXTURE_CENTER_FRACTION = (0.5, 1.1, 0.5)
GREATER_TROCHANTER_FIXTURE_CENTER_FRACTION = (0.5, -0.1, 0.5)
DISTAL_SHAFT_FIXTURE_CENTER_FRACTION = (0.5, 0.5, -0.1)
DISTAL_SHAFT_FIXTURE_SIZE_FRACTION = (1.0, 1.0)
DISTAL_SHAFT_FIXTURE_NORMAL = (0.0, -0.2588190451025207, 0.9659258262890683)
DEFAULT_FEMUR_BBOX_RATIO = (1.0, 1.3, None)
DEFAULT_FEMUR_BBOX_CROP_FROM = (None, "max", None)
DEFAULT_FEMUR_EXPERIMENTAL_RATIO = 1.2
DEFAULT_FEMUR_PROXIMAL_REFERENCE_DISTANCE_MM = 40.0
DEFAULT_FEMUR_PROXIMAL_REFERENCE_WIDTH = "max_xy"
DEFAULT_FEMUR_REFERENCE_MIN_SCALE = (0.8, 0.8, 0.75)
DEFAULT_FEMUR_REFERENCE_MAX_SCALE = (1.2, 1.2, 1.3)
DEFAULT_PMMA_THICKNESS_MM = 10.0
DEFAULT_PMMA_INTRUSION_MM = 6.0
DEFAULT_FEMUR_INPUT_MARGIN_MM = DEFAULT_PMMA_THICKNESS_MM + DEFAULT_PMMA_INTRUSION_MM
DEFAULT_FEMUR_SHAFT_LENGTH_MM = 120.0
DEFAULT_FEMUR_ROUGH_PRE_ICP_LENGTH_MM = 120.0
DEFAULT_FEMUR_CUT_MODE = "post_icp_flat_ratio"
PRE_ICP_CROP_MODES = {"bbox_ratio", "proximal_box_ratio"}
ROUGH_PRE_ICP_CROP_MODES = {"post_icp_flat_ratio"}
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


def target_displacement_percent():
    """Return the maintained hip sideways-fall reporting endpoint."""
    return DEFAULT_FEMUR_TARGET_DISPLACEMENT_PERCENT


def solve_report_profile():
    """Return femur-specific FAIM reporting settings for the shared solver path."""
    return {
        "report_profile": "femur",
        "analysis_var": "fy_ns1",
        "pistoia_vars": ["pis_fy_fail", "pis_stiffy"],
        "failure_axis": "y",
        "default_applied_displacement": DEFAULT_FEMUR_FE_DISPLACEMENT,
        "target_displacement_percent": target_displacement_percent(),
    }


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


def point_transform_to_vtk_matrix(rotation, translation):
    """Return a VTK matrix for ``p_out = p_in @ rotation.T + translation``."""
    import numpy as np
    import vtk

    rotation_arr = np.asarray(rotation, dtype=np.float64)
    translation_arr = np.asarray(translation, dtype=np.float64)
    if rotation_arr.shape != (3, 3):
        raise ValueError("rotation must have shape 3 x 3.")
    if translation_arr.shape != (3,):
        raise ValueError("translation must contain three values.")

    matrix = vtk.vtkMatrix4x4()
    matrix.Identity()
    for row in range(3):
        for col in range(3):
            matrix.SetElement(row, col, float(rotation_arr[row, col]))
        matrix.SetElement(row, 3, float(translation_arr[row]))
    return matrix


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


def _scale_triplet(values, name):
    import numpy as np

    parsed = np.asarray(values, dtype=float)
    if parsed.shape != (3,):
        raise ValueError(f"{name} must contain three x/y/z values.")
    return parsed


def principal_axis_lengths(polydata):
    """Return approximate principal-axis diameters for a VTK polydata surface."""
    try:
        return point_cloud_axis_lengths(polydata_points(polydata))
    except ValueError as exc:
        if "contains no points" in str(exc):
            raise ValueError("Cannot measure principal axes for empty polydata.") from exc
        raise


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

    reference_center = reference_points.mean(axis=0)
    scaled_points = (reference_points - reference_center) * scale + reference_center
    return polydata_from_points(scaled_points), {
        "source": "voxel_surface_point_cloud",
        "reference_center": reference_center.tolist(),
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


def _opposite_crop_end(value):
    if value == "min":
        return "max"
    if value == "max":
        return "min"
    return value


def _parse_bbox_crop_from_recipe(values):
    """Return bbox crop-end values in Ogo's x/y/z image order.

    The recipe uses the same reference/constrained/free convention as the
    Slicer-authored workflow. For the hip sideways-fall recipe the mapped order
    is ``free, reference, constrained`` in Ogo's x/y/z image axes. The workflow
    crop ends are stored in patient-space axis directions, while VTK image
    storage runs opposite on the free and constrained axes used by this recipe;
    those two ends are therefore flipped during the axis-order conversion.
    """
    if values is None:
        return None, None, None
    if len(values) != 3:
        raise ValueError("bbox_crop_from must contain reference/constrained/free values.")
    reference, constrained, free = (_crop_from_value(item) for item in values)
    return _opposite_crop_end(free), reference, _opposite_crop_end(constrained)


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
    """Resample a VTK image using the FE reference-grid output-size rule."""
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
    For the hip recipe this is the exposed shaft face created by
    ``bbox_ratio=(1, 1.3, None)`` and ``bbox_crop_from=(None, "max", None)``.
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
    active_for_crop = _largest_connected_component_mask(active)

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

    spacing = np.asarray(vtk_mask.GetSpacing(), dtype=np.float64)
    out_lo, out_hi, crop_surface, crop_meta = _stable_bbox_ratio_crop_bounds(
        active=active_for_crop,
        spacing=spacing,
        ratio_xyz=ratio_xyz,
        crop_from_xyz=crop_from_xyz,
        reference_axes=reference_axes,
    )

    slices = tuple(slice(int(out_lo[axis]), int(out_hi[axis])) for axis in range(3))
    origin = tuple(
        float(vtk_mask.GetOrigin()[axis]) + float(out_lo[axis]) * float(spacing[axis])
        for axis in range(3)
    )
    cropped_images = [
        _vtk_image_from_array(vtk_image_to_numpy(image)[slices], image, origin=origin)
        for image in images
    ]

    cropped_active = active_for_crop[slices]
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
        "reference_axis": crop_meta["reference_axis"],
        "reference_length_mm": crop_meta["reference_length_mm"],
        "input_bbox_xyz": crop_meta["input_bbox_xyz"],
        "crop_slices_xyz": tuple((int(out_lo[axis]), int(out_hi[axis])) for axis in range(3)),
        "output_shape_xyz": tuple(int(value) for value in cropped_active.shape),
        "output_origin": origin,
        "crop_iterations": crop_meta["crop_iterations"],
        "final_bbox_ratio_xyz": crop_meta["final_bbox_ratio_xyz"],
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


def _largest_connected_component_mask(active):
    import numpy as np

    try:
        from scipy import ndimage
    except ImportError:
        return active

    labeled, count = ndimage.label(active)
    if count <= 1:
        return active
    sizes = np.bincount(labeled.ravel())
    sizes[0] = 0
    return labeled == int(np.argmax(sizes))


def crop_vtk_images_to_proximal_box_ratio(
    vtk_images,
    vtk_mask,
    *,
    ratio=DEFAULT_FEMUR_EXPERIMENTAL_RATIO,
    proximal_reference_distance_mm=DEFAULT_FEMUR_PROXIMAL_REFERENCE_DISTANCE_MM,
    reference_width=DEFAULT_FEMUR_PROXIMAL_REFERENCE_WIDTH,
    labels=None,
):
    """Crop from the distal side using proximal transverse width as reference."""
    import numpy as np

    from ogo.util.vtk_image import vtk_image_to_numpy

    mask_data = vtk_image_to_numpy(vtk_mask)
    active = _active_crop_mask(mask_data, labels)
    spacing = np.asarray(vtk_mask.GetSpacing(), dtype=np.float64)
    coords = np.argwhere(active)
    lo = coords.min(axis=0).astype(np.int64)
    hi = (coords.max(axis=0) + 1).astype(np.int64)
    size = hi - lo
    proximal_voxels = min(
        int(size[2]),
        max(1, int(round(float(proximal_reference_distance_mm) / float(spacing[2])))),
    )
    proximal_start = int(hi[2]) - proximal_voxels
    proximal = active[:, :, proximal_start : int(hi[2])]
    proximal_coords = np.argwhere(proximal)
    if not len(proximal_coords):
        proximal_coords = coords - np.asarray([0, 0, proximal_start], dtype=np.int64)
    x_width_mm = _percentile_index_width_mm(
        proximal_coords[:, 0] + 0,
        float(spacing[0]),
    )
    y_width_mm = _percentile_index_width_mm(
        proximal_coords[:, 1] + 0,
        float(spacing[1]),
    )
    reference_width = str(reference_width).strip().lower()
    if reference_width == "max_xy":
        reference_axis = 0 if x_width_mm >= y_width_mm else 1
        reference_width_mm = x_width_mm if reference_axis == 0 else y_width_mm
    elif reference_width == "rms_xy":
        reference_axis = None
        reference_width_mm = float(np.sqrt((x_width_mm**2 + y_width_mm**2) / 2.0))
    else:
        raise ValueError("reference_width must be max_xy or rms_xy.")
    target_length_mm = float(ratio) * float(reference_width_mm)
    target_voxels = min(int(size[2]), max(1, int(round(target_length_mm / float(spacing[2])))))
    status = "short" if int(size[2]) <= target_voxels else "cropped"
    keep = np.zeros(active.shape, dtype=bool)
    out_lo = lo.copy()
    out_hi = hi.copy()
    if status == "short":
        keep[tuple(slice(int(lo[axis]), int(hi[axis])) for axis in range(3))] = True
    else:
        out_lo[2] = int(hi[2]) - target_voxels
        keep[
            int(out_lo[0]) : int(out_hi[0]),
            int(out_lo[1]) : int(out_hi[1]),
            int(out_lo[2]) : int(out_hi[2]),
        ] = True
    return _crop_vtk_images_with_keep_mask(
        vtk_images,
        vtk_mask,
        active=active,
        keep=keep,
        meta={
            "enabled": True,
            "method": "proximal_box_ratio",
            "ratio": float(ratio),
            "proximal_reference_distance_mm": float(proximal_reference_distance_mm),
            "status": status,
            "reference_width": reference_width,
            "reference_axis": "rms_xy" if reference_axis is None else "xyz"[reference_axis],
            "reference_width_mm": float(reference_width_mm),
            "proximal_x_width_mm": float(x_width_mm),
            "proximal_y_width_mm": float(y_width_mm),
            "target_length_mm": float(target_length_mm),
            "input_bbox_xyz": tuple((int(lo[axis]), int(hi[axis])) for axis in range(3)),
            "crop_slices_xyz": tuple((int(out_lo[axis]), int(out_hi[axis])) for axis in range(3)),
        },
    )


def crop_vtk_images_to_fixed_proximal_length(
    vtk_images,
    vtk_mask,
    *,
    retained_length_mm=DEFAULT_FEMUR_ROUGH_PRE_ICP_LENGTH_MM,
    labels=None,
):
    """Crop from the distal side to a fixed proximal-distal retained length."""
    import numpy as np

    from ogo.util.vtk_image import vtk_image_to_numpy

    mask_data = vtk_image_to_numpy(vtk_mask)
    active = _active_crop_mask(mask_data, labels)
    spacing = np.asarray(vtk_mask.GetSpacing(), dtype=np.float64)
    coords = np.argwhere(active)
    lo = coords.min(axis=0).astype(np.int64)
    hi = (coords.max(axis=0) + 1).astype(np.int64)
    size = hi - lo
    retained_length_mm = float(retained_length_mm)
    if retained_length_mm <= 0.0:
        raise ValueError("retained_length_mm must be positive.")
    target_voxels = min(
        int(size[2]),
        max(1, int(round(retained_length_mm / float(spacing[2])))),
    )
    status = "short" if int(size[2]) <= target_voxels else "cropped"
    keep = np.zeros(active.shape, dtype=bool)
    out_lo = lo.copy()
    out_hi = hi.copy()
    if status == "short":
        keep[tuple(slice(int(lo[axis]), int(hi[axis])) for axis in range(3))] = True
    else:
        out_lo[2] = int(hi[2]) - target_voxels
        keep[
            int(out_lo[0]) : int(out_hi[0]),
            int(out_lo[1]) : int(out_hi[1]),
            int(out_lo[2]) : int(out_hi[2]),
        ] = True
    cropped_images, crop_face_image, meta = _crop_vtk_images_with_keep_mask(
        vtk_images,
        vtk_mask,
        active=active,
        keep=keep,
        meta={
            "enabled": True,
            "method": "fixed_proximal_length",
            "retained_length_mm": float(target_voxels) * float(spacing[2]),
            "requested_retained_length_mm": retained_length_mm,
            "status": status,
            "input_bbox_xyz": tuple((int(lo[axis]), int(hi[axis])) for axis in range(3)),
            "rough_pre_icp": True,
        },
    )
    return cropped_images, crop_face_image, meta


def crop_vtk_images_to_flat_post_icp_ratio(
    vtk_images,
    vtk_mask,
    *,
    ratio=DEFAULT_FEMUR_EXPERIMENTAL_RATIO,
    labels=None,
):
    """Apply a flat aligned-frame distal crop using z length over y width."""
    import numpy as np

    from ogo.util.vtk_image import vtk_image_to_numpy

    mask_data = vtk_image_to_numpy(vtk_mask)
    active = _active_crop_mask(mask_data, labels)
    spacing = np.asarray(vtk_mask.GetSpacing(), dtype=np.float64)
    coords = np.argwhere(active)
    lo = coords.min(axis=0).astype(np.int64)
    hi = (coords.max(axis=0) + 1).astype(np.int64)
    size = hi - lo
    y_width_mm = float(size[1]) * float(spacing[1])
    z_length_mm = float(size[2]) * float(spacing[2])
    ratio = float(ratio)
    if ratio <= 0.0:
        raise ValueError("ratio must be positive.")
    target_length_mm = min(z_length_mm, ratio * y_width_mm)
    target_voxels = min(
        int(size[2]),
        max(1, int(round(target_length_mm / float(spacing[2])))),
    )
    status = "short" if int(size[2]) <= target_voxels else "cropped"
    keep = np.zeros(active.shape, dtype=bool)
    out_lo = lo.copy()
    out_hi = hi.copy()
    if status == "short":
        keep[tuple(slice(int(lo[axis]), int(hi[axis])) for axis in range(3))] = True
    else:
        out_lo[2] = int(hi[2]) - target_voxels
        keep[
            int(out_lo[0]) : int(out_hi[0]),
            int(out_lo[1]) : int(out_hi[1]),
            int(out_lo[2]) : int(out_hi[2]),
        ] = True
    cropped_images, crop_face_image, meta = _crop_vtk_images_with_keep_mask(
        vtk_images,
        vtk_mask,
        active=active,
        keep=keep,
        meta={
            "enabled": True,
            "method": "post_icp_flat_ratio",
            "ratio": ratio,
            "reference_axis": "y",
            "reference_width_mm": y_width_mm,
            "input_z_length_mm": z_length_mm,
            "target_length_mm": float(target_voxels) * float(spacing[2]),
            "status": status,
            "input_bbox_xyz": tuple((int(lo[axis]), int(hi[axis])) for axis in range(3)),
            "crop_stage": "after ICP on aligned reference grid",
        },
    )
    return cropped_images, crop_face_image, meta


def _active_crop_mask(mask_data, labels):
    import numpy as np

    if labels:
        active = np.isin(mask_data, sorted(int(label) for label in labels))
    else:
        active = mask_data != 0
    if not np.any(active):
        raise ValueError("Cannot crop an empty femur mask.")
    return _largest_connected_component_mask(active)


def _percentile_index_width_mm(indices, spacing, *, low=1.0, high=99.0):
    import numpy as np

    lo, hi = np.percentile(indices.astype(np.float64), [float(low), float(high)])
    return float(hi - lo + 1.0) * float(spacing)


def _crop_vtk_images_with_keep_mask(vtk_images, vtk_mask, *, active, keep, meta):
    import numpy as np
    import vtk

    from ogo.util.vtk_image import vtk_image_to_numpy

    kept_active = active & keep
    if not np.any(kept_active):
        raise ValueError("Femur crop removed all foreground voxels.")
    coords = np.argwhere(kept_active)
    out_lo = coords.min(axis=0).astype(np.int64)
    out_hi = (coords.max(axis=0) + 1).astype(np.int64)
    slices = tuple(slice(int(out_lo[axis]), int(out_hi[axis])) for axis in range(3))
    spacing = np.asarray(vtk_mask.GetSpacing(), dtype=np.float64)
    origin = tuple(
        float(vtk_mask.GetOrigin()[axis]) + float(out_lo[axis]) * float(spacing[axis])
        for axis in range(3)
    )
    cropped_images = []
    for image in vtk_images:
        data = vtk_image_to_numpy(image).copy()
        data[~keep] = 0
        cropped_images.append(_vtk_image_from_array(data[slices], image, origin=origin))
    crop_face = kept_active & _touches_removed_or_background(kept_active, active & ~keep)
    crop_face_image = _vtk_image_from_array(
        crop_face[slices].astype(np.uint8),
        vtk_mask,
        origin=origin,
        vtk_array_type=vtk.VTK_UNSIGNED_CHAR,
    )
    meta = dict(meta)
    meta.update(
        {
            "crop_slices_xyz": tuple((int(out_lo[axis]), int(out_hi[axis])) for axis in range(3)),
            "output_shape_xyz": tuple(int(value) for value in kept_active[slices].shape),
            "output_origin": origin,
            "crop_face_voxels": int(crop_face.sum()),
        }
    )
    return cropped_images, crop_face_image, meta


def _touches_removed_or_background(kept_active, removed_active):
    import numpy as np

    face = np.zeros(kept_active.shape, dtype=bool)
    for axis in range(3):
        before = [slice(None), slice(None), slice(None)]
        after = [slice(None), slice(None), slice(None)]
        before[axis] = slice(1, None)
        after[axis] = slice(None, -1)
        face[tuple(before)] |= kept_active[tuple(before)] & removed_active[tuple(after)]
        face[tuple(after)] |= kept_active[tuple(after)] & removed_active[tuple(before)]
    return face


def _stable_bbox_ratio_crop_bounds(
    *,
    active,
    spacing,
    ratio_xyz,
    crop_from_xyz,
    reference_axes,
    max_iterations=8,
):
    import numpy as np

    coords = np.argwhere(active)
    lo = coords.min(axis=0).astype(np.int64)
    hi = (coords.max(axis=0) + 1).astype(np.int64)
    input_lo = lo.copy()
    input_hi = hi.copy()
    out_lo = lo.copy()
    out_hi = hi.copy()
    crop_surface = None
    input_size = input_hi - input_lo
    input_physical_size = input_size.astype(np.float64) * spacing
    reference_axis = min(reference_axes, key=lambda axis: float(input_physical_size[axis]))
    reference_length_mm = 0.0

    for iteration in range(1, int(max_iterations) + 1):
        retained = active[tuple(slice(int(out_lo[axis]), int(out_hi[axis])) for axis in range(3))]
        retained_coords = np.argwhere(retained)
        if not len(retained_coords):
            raise ValueError("bbox-ratio crop removed all foreground voxels.")
        lo = out_lo + retained_coords.min(axis=0).astype(np.int64)
        hi = out_lo + (retained_coords.max(axis=0) + 1).astype(np.int64)
        size = hi - lo
        physical_size = size.astype(np.float64) * spacing
        reference_axis = min(
            reference_axes,
            key=lambda axis: (float(physical_size[axis]), axis != reference_axis),
        )
        reference_length_mm = float(size[reference_axis]) * float(spacing[reference_axis])

        next_lo = lo.copy()
        next_hi = hi.copy()
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
                crop_surface = {"axis": axis, "side": "min", "normal_sign": -1.0}
            elif mode == "max":
                start = int(lo[axis])
                crop_surface = {"axis": axis, "side": "max", "normal_sign": 1.0}
            else:
                center = 0.5 * (float(lo[axis]) + float(hi[axis]))
                start = int(round(center - 0.5 * float(target_voxels)))
            start = max(int(lo[axis]), min(start, int(hi[axis]) - target_voxels))
            next_lo[axis] = start
            next_hi[axis] = start + target_voxels

        if np.array_equal(next_lo, out_lo) and np.array_equal(next_hi, out_hi):
            break
        out_lo, out_hi = next_lo, next_hi

    if crop_surface is not None:
        axis = int(crop_surface["axis"])
        crop_surface["local_index"] = (
            0 if crop_surface["side"] == "min" else int(out_hi[axis] - out_lo[axis] - 1)
        )

    final_size = out_hi - out_lo
    final_physical = final_size.astype(np.float64) * spacing
    final_ratio = tuple(
        float(final_physical[axis] / reference_length_mm) if reference_length_mm > 0 else 0.0
        for axis in range(3)
    )
    return out_lo, out_hi, crop_surface, {
        "reference_axis": "xyz"[reference_axis],
        "reference_length_mm": reference_length_mm,
        "input_bbox_xyz": tuple((int(input_lo[axis]), int(input_hi[axis])) for axis in range(3)),
        "crop_iterations": iteration,
        "final_bbox_ratio_xyz": final_ratio,
    }


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
    size_bounds_axes=None,
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

    if size_bounds_axes is None:
        size = (
            _axis_extent_from_vector(model_bounds, u_axis) * float(size_fraction[0]),
            _axis_extent_from_vector(model_bounds, v_axis) * float(size_fraction[1]),
        )
    else:
        size = _oriented_bbox_face_size(
            model_bounds,
            center=center,
            u_axis=u_axis,
            v_axis=v_axis,
            varying_axes=size_bounds_axes,
            size_fraction=size_fraction,
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


def _oriented_bbox_face_size(
    model_bounds,
    *,
    center,
    u_axis,
    v_axis,
    varying_axes,
    size_fraction,
):
    """Project a bbox face into an oblique plane's local coordinates."""
    import itertools
    import numpy as np

    axis_tokens = {
        "x": 0,
        "y": 1,
        "z": 2,
        0: 0,
        1: 1,
        2: 2,
    }
    axes = []
    for axis in varying_axes:
        key = str(axis).strip().lower() if isinstance(axis, str) else axis
        if key not in axis_tokens:
            raise ValueError("size_bounds_axes values must be x, y, z, 0, 1, or 2.")
        axes.append(axis_tokens[key])
    if not axes:
        raise ValueError("size_bounds_axes must contain at least one axis.")

    center_arr = np.asarray(center, dtype=np.float64)
    corners = []
    for choices in itertools.product((0, 1), repeat=len(axes)):
        point = center_arr.copy()
        for axis, choice in zip(axes, choices):
            lo, hi = _axis_bounds(model_bounds, axis)
            point[axis] = hi if choice else lo
        corners.append(point)
    relative = np.vstack(corners) - center_arr
    u_values = relative @ np.asarray(u_axis, dtype=np.float64)
    v_values = relative @ np.asarray(v_axis, dtype=np.float64)
    return (
        float(np.ptp(u_values)) * float(size_fraction[0]),
        float(np.ptp(v_values)) * float(size_fraction[1]),
    )


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
    candidate_indices = indices[candidates]
    candidate_distances = distances[candidates]
    candidate_u_values = u_values[candidates]
    candidate_v_values = v_values[candidates]
    candidate_lengths = {
        len(candidate_indices),
        len(candidate_distances),
        len(candidate_u_values),
        len(candidate_v_values),
    }
    if len(candidate_lengths) != 1:
        raise ValueError("candidate indices, distances, and plane coordinates must have the same length.")

    for voxel, distance, u_value, v_value in zip(
        candidate_indices,
        candidate_distances,
        candidate_u_values,
        candidate_v_values,
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
    greater_index, greater_z = _peak_center_z(z, y_max, proximal_indices)

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


# -----------------------------------------------------------------------------
# Hip sideways-fall workflow builder
# -----------------------------------------------------------------------------

#####
# Hip workflow builder
#
# This script sets up the sideways fall FE model on the hip from the density (K2HPO4)
# calibrated image. This script sets up the model for either a left or right femur, as
# specified by the user. The analysis resamples the image to isotropic voxels, transforms
# the image, applies the bone mask and bins the data. It then creates the FE model for
# solving using FAIM (>v8.0, Numerics Solutions Ltd, Calgary, Canada - Steven  Boyd).
#
#####
#
# Andrew Michalski
# University of Calgary
# Biomedical Engineering Graduate Program
# April 29, 2019
# Modified to Py3: March 25, 2020
#####

script_version = 1.0

##
# Import the required modules
import ogo.util.Helper as ogo
import os
import sys
import argparse
import time
from datetime import date
from collections import OrderedDict
import numpy as np
import vtk
import vtkbone

from ogo.fea.boundary import (
    fit_vtk_images_to_physical_bounds,
    generate_projected_material_disk_vtk,
    projected_material_disk_required_bounds,
    should_smooth_resampled_mask,
    smooth_binary_mask_vtk,
)
from ogo.fea.materials import build_femur_material_table
from ogo.fea.model import (
    append_postprocessing_sets,
    directional_face_node_ids_from_voxel_mask,
    find_nodes_on_coordinate_plane,
    interface_node_ids_from_voxel_mask,
    write_model,
)
from ogo.util.echo_arguments import echo_arguments


def remove_extension(filename):
    while True:
        filename, ext = os.path.splitext(filename)
        if not ext:
            break
    return filename


##
# Start script
def sidewaysFallFe(args):
    ogo.message("Start of Script...")

    ##
    # Collect the input arguments
    image = args.calibrated_image
    mask = args.bone_mask
    compartment_mask = args.compartment_mask

    mask_threshold = args.mask_threshold
    iso_resolution = args.iso_resolution
    femur_side = args.femur_side
    poissons_ratio = args.poissons_ratio
    pmma_E = args.pmma_E
    pmma_v = args.pmma_v
    pmma_thick = args.pmma_thick
    pmma_intrusion = args.pmma_intrusion
    pmma_mat_id = args.pmma_mat_id
    fe_displacement = args.fe_displacement
    left_femur_reference = args.left_femur_reference
    right_femur_reference = args.right_femur_reference
    output_file = args.output_file
    pmma_yield_compression = args.pmma_yield_compression
    pmma_yield_tension = args.pmma_yield_tension
    femur_shaft_length = args.femur_shaft_length
    femur_cut_mode = args.femur_cut_mode
    femur_bbox_ratio = args.femur_bbox_ratio
    femur_bbox_crop_from = args.femur_bbox_crop_from
    femur_experimental_crop_ratio = args.femur_experimental_crop_ratio
    femur_proximal_reference_distance = args.femur_proximal_reference_distance
    femur_proximal_reference_width = args.femur_proximal_reference_width
    femur_lesser_trochanter_distal_offset = args.femur_lesser_trochanter_distal_offset
    femur_lesser_trochanter_distal_offset_percent = args.femur_lesser_trochanter_distal_offset_percent
    femur_input_margin = args.femur_input_margin
    cortical_label = args.cortical_label
    trabecular_label = args.trabecular_label
    mask_smoothing_spacing_threshold = args.mask_smoothing_spacing_threshold



    ##
    # Determine image locations and names of files
    image_pathname = os.path.dirname(image)
    image_basename = os.path.basename(image)
    mask_pathname = os.path.dirname(mask)
    mask_basename = os.path.basename(mask)
    script_name = sys.argv[0]


    try:
        N88_fileName = sideways_fall_output_name(output_file, femur_side)
    except ValueError:
        ogo.message("Femur side not recognized. Terminating...")
        sys.exit()


    ##
    # Message the input parameters to the terminal
    ogo.message("Image Path: %s" % image_pathname)
    ogo.message("Image File: %s" % image_basename)
    ogo.message("Mask Path: %s" % mask_pathname)
    ogo.message("Mask File: %s" % mask_basename)
    ogo.message("Isotropic Voxel Size: %8.4f" % iso_resolution)
    if femur_side == 1:
        ogo.message("Femur Side for Model: Left")
    elif femur_side == 2:
        ogo.message("Femur Side for Model: Right")
    else:
        ogo.message("Femur side not recognized. Terminating...")
        sys.exit()
    ogo.message("Trabecular Elastic Function: %s" % (args.elastic_E_func or "default_E"))
    ogo.message("Cortical Elastic Function: %s" % (args.cort_elastic_E_func or args.elastic_E_func or "default_E"))
    ogo.message("Bone Poissons Ratio: %1.1f" % poissons_ratio)
    ogo.message("PMMA Elastic Modulus [MPa]: %8.4f" % pmma_E)
    ogo.message("PMMA Poissons Ration: %1.1f" % pmma_v)
    ogo.message("PMMA Thickness [mm]: %8.4f" % pmma_thick)
    ogo.message("PMMA Intrusion [mm]: %8.4f" % pmma_intrusion)
    ogo.message("PMMA Material ID: %d" % pmma_mat_id)
    ogo.message("Applied Displacement [mm]: %8.4f" % fe_displacement)
    ogo.message("Femur Cut Mode: %s" % femur_cut_mode)
    if femur_cut_mode == "bbox_ratio":
        ogo.message("Femur BBox Ratio [reference, constrained, free]: %s" % str(femur_bbox_ratio))
        ogo.message("Femur BBox Crop From [reference, constrained, free]: %s" % str(femur_bbox_crop_from))
    if femur_cut_mode in {"proximal_box_ratio", "post_icp_flat_ratio"}:
        ogo.message("Femur Proximal Box Crop Ratio: %8.4f" % femur_experimental_crop_ratio)
        ogo.message("Femur Proximal Reference Distance [mm]: %8.4f" % femur_proximal_reference_distance)
        ogo.message("Femur Proximal Reference Width: %s" % femur_proximal_reference_width)
    if femur_cut_mode == "post_icp_flat_ratio":
        ogo.message("Femur Rough Pre-ICP Retained Length [mm]: %8.4f" % femur_shaft_length)
    if compartment_mask is not None:
        ogo.message("Compartment Mask: %s" % compartment_mask)
        ogo.message("Cortical Label: %d" % cortical_label)
        ogo.message("Trabecular Label: %d" % trabecular_label)
    if femur_cut_mode == "fixed_length":
        ogo.message("Retained Proximal Femur Length [mm]: %8.4f" % femur_shaft_length)
    elif femur_lesser_trochanter_distal_offset_percent is not None:
        ogo.message(
            "Lesser Trochanter Distal Cut Offset [%% greater-lesser distance]: %8.4f"
            % femur_lesser_trochanter_distal_offset_percent
        )
    else:
        ogo.message(
            "Lesser Trochanter Distal Cut Offset [mm]: %8.4f"
            % femur_lesser_trochanter_distal_offset
        )
    ogo.message("Input Foreground Safety Margin [mm]: %8.4f" % femur_input_margin)

    ##
    # Read input image
    ogo.message("Reading calibrated image...")
    imageData = ogo.readNii(image)
    input_spacing = imageData.GetSpacing()

    ##
    # Read bone mask
    ogo.message("Reading bone mask...")
    maskData = ogo.readNii(mask)
    maskThres = ogo.maskThreshold(maskData, mask_threshold)
    compartmentData = None
    if compartment_mask is not None:
        ogo.message("Reading trabecular/cortical compartment mask...")
        compartmentData = ogo.readNii(compartment_mask)

    distal_crop_face = None
    bbox_crop_meta = None

    ogo.message("Resampling femur inputs to isotropic spacing before custom crop and ICP...")
    imageData = resample_vtk_image_like_workflow(
        imageData,
        iso_resolution,
        interpolation="bspline",
    )
    maskThres = resample_vtk_image_like_workflow(
        maskThres,
        iso_resolution,
        interpolation="nearest",
    )
    if distal_crop_face is not None:
        distal_crop_face = resample_vtk_image_like_workflow(
            distal_crop_face,
            iso_resolution,
            interpolation="nearest",
        )
    if compartmentData is not None:
        compartmentData = resample_vtk_image_like_workflow(
            compartmentData,
            iso_resolution,
            interpolation="nearest",
        )

    if femur_cut_mode in PRE_ICP_CROP_MODES:
        ogo.message("Applying %s femur crop after isotropic resampling and before ICP..." % femur_cut_mode)
        images_to_crop = [imageData, maskThres] + ([compartmentData] if compartmentData is not None else [])
        if femur_cut_mode == "bbox_ratio":
            cropped_images, distal_crop_face, bbox_crop_meta = crop_vtk_images_to_bbox_ratio(
                images_to_crop,
                maskThres,
                bbox_ratio=femur_bbox_ratio,
                bbox_crop_from=femur_bbox_crop_from,
                labels={1},
            )
        else:
            cropped_images, distal_crop_face, bbox_crop_meta = crop_vtk_images_to_proximal_box_ratio(
                images_to_crop,
                maskThres,
                ratio=femur_experimental_crop_ratio,
                proximal_reference_distance_mm=femur_proximal_reference_distance,
                reference_width=femur_proximal_reference_width,
                labels={1},
            )
        imageData = cropped_images[0]
        maskThres = cropped_images[1]
        if compartmentData is not None:
            compartmentData = cropped_images[2]
        ogo.message(
            "%s crop slices xyz=%s; output shape xyz=%s; crop face voxels=%d."
            % (
                femur_cut_mode,
                bbox_crop_meta["crop_slices_xyz"],
                bbox_crop_meta["output_shape_xyz"],
                bbox_crop_meta["crop_face_voxels"],
            )
        )
    elif femur_cut_mode in ROUGH_PRE_ICP_CROP_MODES:
        ogo.message(
            "Applying fixed-length rough femur crop after isotropic resampling and before ICP..."
        )
        images_to_crop = [imageData, maskThres] + ([compartmentData] if compartmentData is not None else [])
        cropped_images, _rough_crop_face, bbox_crop_meta = crop_vtk_images_to_fixed_proximal_length(
            images_to_crop,
            maskThres,
            retained_length_mm=femur_shaft_length,
            labels={1},
        )
        imageData = cropped_images[0]
        maskThres = cropped_images[1]
        if compartmentData is not None:
            compartmentData = cropped_images[2]
        ogo.message(
            "rough pre-ICP crop slices xyz=%s; output shape xyz=%s; retained length=%8.4f; status=%s."
            % (
                bbox_crop_meta["crop_slices_xyz"],
                bbox_crop_meta["output_shape_xyz"],
                bbox_crop_meta["retained_length_mm"],
                bbox_crop_meta["status"],
            )
        )
    registration_mask = maskThres

    images_to_pad = [imageData, maskThres]
    if distal_crop_face is not None:
        images_to_pad.append(distal_crop_face)
    if compartmentData is not None:
        images_to_pad.append(compartmentData)
    pad_constants = [0] * len(images_to_pad)
    padded_images, padding = pad_vtk_images_to_foreground_margin(
        images_to_pad,
        maskThres,
        margin_mm=femur_input_margin,
        constants=pad_constants,
    )
    imageData = padded_images[0]
    maskThres = padded_images[1]
    next_padded_index = 2
    if distal_crop_face is not None:
        distal_crop_face = padded_images[next_padded_index]
        next_padded_index += 1
    if compartmentData is not None:
        compartmentData = padded_images[next_padded_index]
    if any(padding["lower"]) or any(padding["upper"]):
        ogo.message(
            "Padded isotropic input image extent by lower=%s upper=%s voxels for FE transform safety."
            % (padding["lower"], padding["upper"])
        )
    else:
        ogo.message("Isotropic input image already has sufficient foreground safety margin.")

    ##
    # Use the cropped, isotropic input geometry directly for ICP.
    image_rot = imageData
    mask_rot = maskThres
    distal_crop_face_rot = distal_crop_face
    compartment_rot = compartmentData

    ##
    # Align the input femur with the reference model
    ogo.message("Aligning input with reference model...")
    sample_surface_points = surface_points_from_vtk_mask(
        registration_mask,
        max_points=8000,
        sample_mode="stride",
        sample_offset=0,
    )
    mask_surface = polydata_from_points(sample_surface_points)
    if femur_side == 1:
        ref_poly = ogo.readPolyData(left_femur_reference)
    elif femur_side == 2:
        if os.path.exists(right_femur_reference):
            ref_poly = ogo.readPolyData(right_femur_reference)
        else:
            ogo.message(
                "Right femur reference not found; mirroring left reference in x:",
                right_femur_reference,
            )
            ref_poly = mirror_polydata_x(ogo.readPolyData(left_femur_reference))
    else:
        print("Error: Femur Side not defined. Terminating...")
        sys.exit()

    ref_poly, reference_scale = scale_reference_point_cloud_to_sample(
        ref_poly,
        sample_surface_points,
        max_points=8000,
        reference_sample_mode="linspace",
        min_scale=DEFAULT_FEMUR_REFERENCE_MIN_SCALE,
        max_scale=DEFAULT_FEMUR_REFERENCE_MAX_SCALE,
    )
    ogo.message("Femur reference axis lengths: %s" % str(reference_scale["reference_axis_lengths"]))
    ogo.message("Sample femur axis lengths: %s" % str(reference_scale["sample_axis_lengths"]))
    ogo.message("Femur reference scale factors: %s" % str(reference_scale["scale_factors"]))

    icp_transform = estimate_rigid_icp(
        moving_points=polydata_points(ref_poly),
        fixed_points=sample_surface_points,
        iterations=50,
        tolerance=1.0e-4,
        start_by_matching_centroids_only=True,
        convergence="delta",
        distance_mode="mean",
    )
    icp = point_transform_to_vtk_matrix(
        icp_transform["rotation"],
        icp_transform["translation"],
    )
    ogo.message(
        "ICP reference-to-sample iterations=%d mean_distance=%0.4f"
        % (icp_transform["iterations"], icp_transform["mean_distance"])
    )

    ogo.message("Applying the transformation and isotropic resampling to the image and mask...")
    output_origin, output_size = reference_grid_from_vtk_mask(
        mask_rot,
        icp,
        margin_voxels=int(round(femur_input_margin / max(iso_resolution, 1.0e-6))),
    )
    output_spacing = (float(iso_resolution), float(iso_resolution), float(iso_resolution))
    ogo.message(
        "Reference-frame output grid: origin=%s size=%s spacing=%s"
        % (output_origin, output_size, output_spacing)
    )
    image_trans = transform_resample_vtk_image_to_reference_grid(
        image_rot,
        icp,
        output_origin=output_origin,
        output_size=output_size,
        output_spacing=output_spacing,
        interpolation="bspline",
    )
    mask_trans = transform_resample_vtk_image_to_reference_grid(
        mask_rot,
        icp,
        output_origin=output_origin,
        output_size=output_size,
        output_spacing=output_spacing,
        interpolation="nearest",
    )
    distal_crop_face_trans = (
        transform_resample_vtk_image_to_reference_grid(
            distal_crop_face_rot,
            icp,
            output_origin=output_origin,
            output_size=output_size,
            output_spacing=output_spacing,
            interpolation="nearest",
        )
        if distal_crop_face_rot is not None
        else None
    )
    compartment_trans = (
        transform_resample_vtk_image_to_reference_grid(
            compartment_rot,
            icp,
            output_origin=output_origin,
            output_size=output_size,
            output_spacing=output_spacing,
            interpolation="nearest",
        )
        if compartment_rot is not None
        else None
    )
    smooth_resampled_masks = should_smooth_resampled_mask(
        input_spacing,
        mask_smoothing_spacing_threshold,
    )
    if smooth_resampled_masks:
        ogo.message(
            "smoothing resampled femur mask because input spacing "
            f"{input_spacing} exceeds {mask_smoothing_spacing_threshold} mm..."
        )
        mask_trans = smooth_binary_mask_vtk(mask_trans, close_iter=1, open_iter=1)
    else:
        ogo.message(
            "skipping femur mask smoothing because input spacing "
            f"{input_spacing} is <= {mask_smoothing_spacing_threshold} mm in all dimensions..."
        )

    if femur_cut_mode in PRE_ICP_CROP_MODES:
        ogo.message("Using transformed %s crop face for distal shaft standardization." % femur_cut_mode)
        try:
            mask_z_min, mask_z_max = femur_z_coverage(mask_trans)
        except ValueError as exc:
            ogo.message(str(exc))
            sys.exit(1)
        retained_length_mm = mask_z_max - mask_z_min
        shaft_crop = {
            "cut_mode": femur_cut_mode,
            "cut_plane": "transformed input crop face",
            "pre_icp_crop": bbox_crop_meta,
            "mask_z_min": mask_z_min,
            "mask_z_max": mask_z_max,
            "retained_length_mm": retained_length_mm,
            "warnings": [],
        }
        ogo.message(
            "Transformed %s mask z coverage [%8.4f, %8.4f]; retained z-span=%8.4f"
            % (femur_cut_mode, mask_z_min, mask_z_max, retained_length_mm)
        )
    elif femur_cut_mode == "post_icp_flat_ratio":
        ogo.message("Applying flat post-ICP femur crop in aligned reference coordinates...")
        images_to_crop = [image_trans, mask_trans] + ([compartment_trans] if compartment_trans is not None else [])
        try:
            cropped_images, distal_crop_face_trans, shaft_crop = crop_vtk_images_to_flat_post_icp_ratio(
                images_to_crop,
                mask_trans,
                ratio=femur_experimental_crop_ratio,
                labels={1},
            )
        except ValueError as exc:
            ogo.message(str(exc))
            sys.exit(1)
        image_trans = cropped_images[0]
        mask_trans = cropped_images[1]
        if compartment_trans is not None:
            compartment_trans = cropped_images[2]
        retained_length_mm = shaft_crop["target_length_mm"]
        shaft_crop["pre_icp_crop"] = bbox_crop_meta
        ogo.message(
            "Post-ICP flat crop retained z=%8.4f from y width=%8.4f; slices xyz=%s; status=%s."
            % (
                retained_length_mm,
                shaft_crop["reference_width_mm"],
                shaft_crop["crop_slices_xyz"],
                shaft_crop["status"],
            )
        )
    else:
        ogo.message("Applying flat distal femur crop in reference coordinates...")
        try:
            image_trans, mask_trans, shaft_crop = standardize_femur_shaft_length(
                image_trans,
                mask_trans,
                reference_bounds=ref_poly.GetBounds(),
                retained_length_mm=femur_shaft_length,
                cut_mode=femur_cut_mode,
                lesser_trochanter_distal_offset_mm=femur_lesser_trochanter_distal_offset,
                lesser_trochanter_distal_offset_percent=femur_lesser_trochanter_distal_offset_percent,
            )
        except ValueError as exc:
            ogo.message(str(exc))
            sys.exit(1)
        retained_length_mm = shaft_crop["retained_length_mm"]
        ogo.message(
            "Distal cut z=%8.4f; input mask z coverage [%8.4f, %8.4f]; retained length=%8.4f"
            % (
                shaft_crop["cut_z"],
                shaft_crop["mask_z_min"],
                shaft_crop["mask_z_max"],
                retained_length_mm,
            )
        )
        for warning in shaft_crop.get("warnings", []):
            ogo.message("WARNING: Femur crop QC: %s" % warning)
        if femur_cut_mode == "lesser_trochanter":
            ogo.message(
                "Detected greater trochanter z=%8.4f; lesser trochanter z=%8.4f"
                % (shaft_crop["greater_trochanter_z"], shaft_crop["lesser_trochanter_z"])
            )
        if compartment_trans is not None:
            compartment_trans = flat_crop_vtk_image_below_z(compartment_trans, shaft_crop["cut_z"])

    # ogo.message("Writing out temp images...")
    # ogo.writeNii(image_trans, "temp_image.nii", image_pathname)
    # ogo.writeNii(mask_trans, "temp_mask.nii", image_pathname)
    #
    # left_femur_reference_image = "/Users/andrew.michalski/Development/Ogo/LT_FEMUR_SIDEWAYS_FALL_REF.nii"
    # right_femur_reference_image = "/Users/andrew.michalski/Development/Ogo/RT_FEMUR_SIDEWAYS_FALL_REF.nii"
    # ogo.message("Applying SimpleITK registration of the images...")
    # if femur_side == 1:
    #     ogo.finalRegistration(left_femur_reference_image)
    # elif femur_side == 2:
    #     ogo.finalRegistration(right_femur_reference_image)
    # ogo.message("Reading temp images...")
    # image_trans = ogo.readNii("temp_image2.nii")
    # mask_trans = ogo.readNii("temp_mask2.nii")
    # os.remove("temp_image2.nii")
    # os.remove("temp_mask2.nii")

    # replace any negative values less than -31 to be equivalent to -31.
    # -31 is used as that converts to a minumum elastic modulus value of 0.1 MPa
    # K2HPO4 den = -31 mg/cc => Ash den = 6 mg/cc => E = 0.1 MPa
    image_thres = ogo.bmd_preprocess(image_trans, -31)
    image_ash = ogo.bmd_K2hpo4ToAsh(image_thres)

    cortical_mask = None
    if compartment_trans is not None:
        cortical_mask = cortical_compartment_mask(
            compartment_trans,
            cortical_label=cortical_label,
            trabecular_label=trabecular_label,
        )
        if smooth_resampled_masks:
            ogo.message("smoothing derived cortical compartment mask...")
            cortical_mask = smooth_binary_mask_vtk(cortical_mask, close_iter=1, open_iter=1)
    binned_image, bin_centers = ogo.density2materialID(
        image_ash,
        n_bins=128,
        cort_mask=cortical_mask,
    )

    ##
    # Set up the FE model.
    ogo.message("Setting up the Finite Element Model...")
    # Cast the image to Short to "round" float values to nearest whole number
    cast_image = ogo.cast2short(binned_image)
    cast_mask = ogo.cast2unsignchar(mask_trans)


    # Apply the mask to the bone
    ogo.message("Applying the bone mask to the image...")
    bone_image = ogo.applyMask(cast_image, cast_mask)

    # Ensure connectivity
    ogo.message("Performing Image connectivity...")
    conn = ogo.imageConnectivity(bone_image)


    change = ogo.prepareFiniteElementImage(conn)
    femur_mask_on_model_grid = ogo.prepareFiniteElementImage(cast_mask)
    distal_crop_face_change = (
        ogo.prepareFiniteElementImage(distal_crop_face_trans) if distal_crop_face_trans is not None else None
    )
    if distal_crop_face_change is not None:
        distal_crop_face_change = ogo.applyMask(distal_crop_face_change, ogo.cast2unsignchar(change))
    if femur_cut_mode in PRE_ICP_CROP_MODES:
        ogo.message("Skipping model-grid distal shaft cut; pre-ICP crop face is already transformed.")
    elif femur_cut_mode == "post_icp_flat_ratio":
        ogo.message("Skipping model-grid distal shaft cut; post-ICP crop face is already on the model grid.")
    else:
        model_z_min, model_z_max = femur_z_coverage(change)
        model_cut_z = model_z_max - retained_length_mm
        if model_z_min > model_cut_z:
            if model_z_min - model_cut_z <= 2.0:
                model_cut_z = model_z_min
            else:
                ogo.message(
                    "Femur scan does not include enough distal shaft after model-grid alignment "
                    "(required retained length %.1f mm, available approximately %.1f mm)."
                    % (retained_length_mm, max(0.0, model_z_max - model_z_min))
                )
                sys.exit(1)
        change = flat_crop_vtk_image_below_z(change, model_cut_z)
        femur_mask_on_model_grid = flat_crop_vtk_image_below_z(
            femur_mask_on_model_grid,
            model_cut_z,
        )
        ogo.message(
            "Applied model-grid distal shaft cut z=%8.4f; model mask z coverage [%8.4f, %8.4f]"
            % (model_cut_z, model_z_min, model_z_max)
        )

    # Convert image data to hexahedral elements
    ogo.message("Meshing image data to elements...")
    mesh = ogo.Image2Mesh(change)

    ##
    # Set up the Material Table
    ogo.message("Setting up the Finite Element Material Table...")
    material_table = build_femur_material_table(
        bin_centers,
        n_bins=128,
        elastic_E_func=args.elastic_E_func,
        yield_comp_func=args.yield_comp_func,
        yield_tens_func=args.yield_tens_func,
        cort_elastic_E_func=args.cort_elastic_E_func,
        cort_yield_comp_func=args.cort_yield_comp_func,
        cort_yield_tens_func=args.cort_yield_tens_func,
        cort_poissons_ratio=args.cort_poissons_ratio,
        poissons_ratio=poissons_ratio,
        pmma_mat_id=pmma_mat_id,
        pmma_E=pmma_E,
        pmma_v=pmma_v,
        pmma_yield_tension=pmma_yield_tension,
        pmma_yield_compression=pmma_yield_compression,
        include_cortical=cortical_mask is not None,
    )
    ##
    # Create the preliminary FE model
    ogo.message("Constructing the Finite Element Model...")
    model = ogo.applyTestBase(mesh, material_table)
    model.ComputeBounds()
    mesh_model_bounds = model.GetBounds()
    model_bounds = foreground_voxel_center_bounds(femur_mask_on_model_grid)
    ogo.message("Model Mesh Bounds: %s" % str(mesh_model_bounds))
    ogo.message("Model Foreground Voxel-Center Bounds: %s" % str(model_bounds))

    # Define support vectors for the two PMMA contact fixtures.
    top_support_vector = (0, 1, 0)
    bottom_support_vector = (0, -1, 0)

    ##
    # Create PMMA disks from bbox-relative contact planes. The disk normal
    # points from the contact plane toward the anatomy; the generated PMMA
    # image clears any femur-mask voxels so PMMA never occupies the segmentation.
    femoral_head_plane = bbox_relative_fixture_plane(
        model_bounds,
        center_fraction=FEMORAL_HEAD_FIXTURE_CENTER_FRACTION,
        size_fraction=SIDEWAYS_FALL_FIXTURE_SIZE_FRACTION,
        projection_axis="y",
        shape=SIDEWAYS_FALL_FIXTURE_SHAPE,
    )
    femoral_head_cap_direction = bbox_relative_fixture_direction(
        FEMORAL_HEAD_FIXTURE_CENTER_FRACTION,
        projection_axis="y",
    )
    ogo.message(
        "Femoral Head PMMA Plane: center=%s normal=%s u=%s v=%s size=%s"
        % (
            femoral_head_plane["center"],
            femoral_head_plane["normal"],
            femoral_head_plane["u_axis"],
            femoral_head_plane["v_axis"],
            femoral_head_plane["size"],
        )
    )
    ogo.message("Femoral Head PMMA cap direction: %s" % femoral_head_cap_direction)

    ##
    # Creates the greater trochanter pmma cap
    greater_trochanter_plane = bbox_relative_fixture_plane(
        model_bounds,
        center_fraction=GREATER_TROCHANTER_FIXTURE_CENTER_FRACTION,
        size_fraction=SIDEWAYS_FALL_FIXTURE_SIZE_FRACTION,
        projection_axis="y",
        shape=SIDEWAYS_FALL_FIXTURE_SHAPE,
    )
    greater_trochanter_cap_direction = bbox_relative_fixture_direction(
        GREATER_TROCHANTER_FIXTURE_CENTER_FRACTION,
        projection_axis="y",
    )
    ogo.message(
        "Greater Trochanter PMMA Plane: center=%s normal=%s u=%s v=%s size=%s"
        % (
            greater_trochanter_plane["center"],
            greater_trochanter_plane["normal"],
            greater_trochanter_plane["u_axis"],
            greater_trochanter_plane["v_axis"],
            greater_trochanter_plane["size"],
        )
    )
    ogo.message("Greater Trochanter PMMA cap direction: %s" % greater_trochanter_cap_direction)

    required_bounds = [
        projected_material_disk_required_bounds(
            femur_mask_on_model_grid,
            center=plane["center"],
            normal=plane["normal"],
            u_axis=plane["u_axis"],
            v_axis=plane["v_axis"],
            size=plane["size"],
            shape=plane["shape"],
            thickness=pmma_thick,
            intrusion=pmma_intrusion,
        )
        for plane in (femoral_head_plane, greater_trochanter_plane)
    ]
    required_bounds = [bounds for bounds in required_bounds if bounds is not None]
    if required_bounds:
        bounds_array = np.asarray([model_bounds, *required_bounds], dtype=float)
        desired_bounds = (
            float(bounds_array[:, 0].min()),
            float(bounds_array[:, 1].max()),
            float(bounds_array[:, 2].min()),
            float(bounds_array[:, 3].max()),
            float(bounds_array[:, 4].min()),
            float(bounds_array[:, 5].max()),
        )
        images_to_fit = [change, femur_mask_on_model_grid]
        if distal_crop_face_change is not None:
            images_to_fit.append(distal_crop_face_change)
        fitted_images, contact_fit = fit_vtk_images_to_physical_bounds(
            images_to_fit,
            desired_bounds=desired_bounds,
            constants=[0] * len(images_to_fit),
        )
        change = fitted_images[0]
        femur_mask_on_model_grid = fitted_images[1]
        if distal_crop_face_change is not None:
            distal_crop_face_change = fitted_images[2]
        model_bounds = foreground_voxel_center_bounds(femur_mask_on_model_grid)
        ogo.message(
            "Fit model canvas to projected PMMA contact bounds: output_extent=%s"
            % (contact_fit["output_extent"],)
        )

    femoralHeadPMMA = generate_projected_material_disk_vtk(
        change,
        surface_vtk_image=femur_mask_on_model_grid,
        exclusion_vtk_image=femur_mask_on_model_grid,
        center=femoral_head_plane["center"],
        normal=femoral_head_plane["normal"],
        u_axis=femoral_head_plane["u_axis"],
        v_axis=femoral_head_plane["v_axis"],
        size=femoral_head_plane["size"],
        shape=femoral_head_plane["shape"],
        thickness=pmma_thick,
        intrusion=pmma_intrusion,
        anatomy_constrained=True,
        output_value=pmma_mat_id,
    )
    greaterTrochanterPMMA = generate_projected_material_disk_vtk(
        change,
        surface_vtk_image=femur_mask_on_model_grid,
        exclusion_vtk_image=femur_mask_on_model_grid,
        center=greater_trochanter_plane["center"],
        normal=greater_trochanter_plane["normal"],
        u_axis=greater_trochanter_plane["u_axis"],
        v_axis=greater_trochanter_plane["v_axis"],
        size=greater_trochanter_plane["size"],
        shape=greater_trochanter_plane["shape"],
        thickness=pmma_thick,
        intrusion=pmma_intrusion,
        anatomy_constrained=True,
        output_value=pmma_mat_id,
    )

    ##
    # combines the pmma cap images with the model image.
    ogo.message("Combine PMMA Cap Images with Model Image...")
    combinedImage = ogo.combineImageData_SF(change, femoralHeadPMMA, greaterTrochanterPMMA, pmma_mat_id)

    ##
    # Mesh the final image and create Finite Element Model
    # Ensure connectivity
    ogo.message("Performing Image connectivity...")
    conn2 = ogo.imageConnectivity(combinedImage)

    # Convert image data to hexahedral elements
    ogo.message("Meshing image data to elements...")
    mesh2 = ogo.Image2Mesh(conn2)

    ##
    # Create the final FE model
    ogo.message("Constructing the Finite Element Model...")
    model2 = ogo.applyTestBase(mesh2, material_table)
    model2.ComputeBounds()
    model2_bounds = model2.GetBounds()
    ogo.message("Model 2 Bounds: %s" % str(model2_bounds))

    ##
    # Determine Femoral Head PMMA Cap support nodes.
    ogo.message("Determining Femoral Head PMMA Cap nodes...")
    if femoral_head_cap_direction == "up":
        fh_pmma_bounds = (
            model2_bounds[0],
            model2_bounds[1],
            model2_bounds[3] - 1,
            model2_bounds[3],
            model2_bounds[4],
            model2_bounds[5],
        )
        femoral_head_support_vector = top_support_vector
    else:
        fh_pmma_bounds = (
            model2_bounds[0],
            model2_bounds[1],
            model2_bounds[2],
            model2_bounds[2] + 1,
            model2_bounds[4],
            model2_bounds[5],
        )
        femoral_head_support_vector = bottom_support_vector
    fh_pmma_visible_node_IDS = directional_face_node_ids_from_voxel_mask(
        model2,
        femoralHeadPMMA,
        direction=femoral_head_support_vector,
        name=FEMORAL_HEAD_NODE_SET,
    )

    ogo.message("-- found %d outer-face nodes on Femoral Head PMMA Cap."
        % fh_pmma_visible_node_IDS.GetNumberOfTuples())
    model2.AddNodeSet(fh_pmma_visible_node_IDS)

    ##
    # Determine Greater Trochanter PMMA Cap support nodes
    ogo.message("Determining Greater Trochanter PMMA Cap support nodes...")
    if greater_trochanter_cap_direction == "up":
        gt_pmma_bounds = (
            model2_bounds[0],
            model2_bounds[1],
            model2_bounds[3] - 1,
            model2_bounds[3],
            model2_bounds[4],
            model2_bounds[5],
        )
        greater_trochanter_support_vector = top_support_vector
    else:
        gt_pmma_bounds = (
            model2_bounds[0],
            model2_bounds[1],
            model2_bounds[2],
            model2_bounds[2] + 1,
            model2_bounds[4],
            model2_bounds[5],
        )
        greater_trochanter_support_vector = bottom_support_vector
    gt_pmma_visible_node_IDS = directional_face_node_ids_from_voxel_mask(
        model2,
        greaterTrochanterPMMA,
        direction=greater_trochanter_support_vector,
        name=GREATER_TROCHANTER_NODE_SET,
    )

    ogo.message("-- found %d outer-face nodes on Greater Trochanter PMMA Cap."
        % gt_pmma_visible_node_IDS.GetNumberOfTuples())

    model2.AddNodeSet(gt_pmma_visible_node_IDS)

    ##
    # Distal Femur (df) support nodes. The bbox-ratio crop standardizes shaft
    # length before registration; the support itself remains a bbox-relative
    # recipe plane in the generated model frame.
    ogo.message("Determining distal femur nodes...")
    if femur_cut_mode in PRE_ICP_CROP_MODES:
        if distal_crop_face_change is None:
            ogo.message("%s cut mode requires a transformed distal crop-face mask." % femur_cut_mode)
            sys.exit(1)
        distal_plane = bbox_relative_oriented_contact_plane(
            model_bounds,
            center_fraction=DISTAL_SHAFT_FIXTURE_CENTER_FRACTION,
            size_fraction=DISTAL_SHAFT_FIXTURE_SIZE_FRACTION,
            normal=DISTAL_SHAFT_FIXTURE_NORMAL,
            size_bounds_axes=("x", "y"),
            shape="anatomy",
        )
        ogo.message(
            "Distal Femur bbox-relative oblique support plane: center=%s normal=%s outward=%s size=%s"
            % (
                distal_plane["center"],
                distal_plane["normal"],
                distal_plane["outward_normal"],
                distal_plane["size"],
            )
        )
        distal_surface = projected_crop_face_surface_vtk(
            change,
            distal_plane,
            intrusion=pmma_intrusion,
            output_value=1,
        )
        df_visible_node_IDS = interface_node_ids_from_voxel_mask(
            model2,
            distal_surface,
            change,
            name=DISTAL_FEMUR_NODE_SET,
            direction=distal_plane["normal"],
        )
        if df_visible_node_IDS.GetNumberOfTuples() == 0:
            ogo.message("No distal femur nodes found on the bbox-relative oblique shaft support surface.")
            sys.exit(1)
    elif femur_cut_mode == "post_icp_flat_ratio":
        if distal_crop_face_change is None:
            ogo.message("post_icp_flat_ratio cut mode requires a distal crop-face mask.")
            sys.exit(1)
        distal_plane = crop_face_contact_plane(distal_crop_face_change, change)
        ogo.message(
            "Distal Femur post-ICP flat support plane: center=%s normal=%s outward=%s size=%s"
            % (
                distal_plane["center"],
                distal_plane["normal"],
                distal_plane["outward_normal"],
                distal_plane["size"],
            )
        )
        distal_surface = projected_crop_face_surface_vtk(
            change,
            distal_plane,
            intrusion=pmma_intrusion,
            output_value=1,
        )
        df_visible_node_IDS = interface_node_ids_from_voxel_mask(
            model2,
            distal_surface,
            change,
            name=DISTAL_FEMUR_NODE_SET,
            direction=distal_plane["normal"],
        )
        if df_visible_node_IDS.GetNumberOfTuples() == 0:
            ogo.message("No distal femur nodes found on the post-ICP flat shaft support surface.")
            sys.exit(1)
    else:
        distal_cut_z = model2_bounds[4]
        ogo.message("Distal Femur Cut Plane z: %8.4f" % distal_cut_z)
        df_visible_node_IDS = find_nodes_on_coordinate_plane(model2, "z", distal_cut_z)

    ogo.message("-- found %d distal shaft interface nodes."
        % df_visible_node_IDS.GetNumberOfTuples())

    df_visible_node_IDS.SetName(DISTAL_FEMUR_NODE_SET)
    model2.AddNodeSet(df_visible_node_IDS)

    ##
    # Apply Boundary conditions to PMMA caps at specific sites
    ogo.message("Applying displacement boundary condition to Femoral Head PMMA cap...")
    model2.ApplyBoundaryCondition(
        FEMORAL_HEAD_NODE_SET,
        vtkbone.vtkboneConstraint.SENSE_Y,
        fe_displacement,
        "top_displacement")

    ogo.message("Constraining Greater Trochanter PMMA cap in loading direction...")
    model2.ApplyBoundaryCondition(
        GREATER_TROCHANTER_NODE_SET,
        vtkbone.vtkboneConstraint.SENSE_Y,
        0,
        "bottom_fixed_y_PMMA")

    for SENSE, label in (
        (vtkbone.vtkboneConstraint.SENSE_X, "x"),
        (vtkbone.vtkboneConstraint.SENSE_Z, "z"),
    ):
        ogo.message("Constraining distal femur rigid-body motion...")
        model2.ApplyBoundaryCondition(
            DISTAL_FEMUR_NODE_SET,
            SENSE,
            0,
            f"bottom_fixed_{label}")

    ##
    # Post Processing parameters
    ogo.message("Setting up Post Processing Parameters...")
    append_postprocessing_sets(
        model2,
        SIDEWAYS_FALL_NODE_SETS,
    )

    model2.AppendHistory(
        "Created by %s version %s." % (script_name, script_version))

    ##
    # Write out n88model file
    ogo.message("Writing out n88model file: %s" % N88_fileName)
    write_model(model2, N88_fileName)

    ##
    ogo.message("Done writing n88model.")

    ##
    # End of script
    ogo.message("End of Script.")
    sys.exit()


def main():
     # Setup description
    description='''
This script sets up the sideways fall FE model on the hip from the
        density (K2HPO4) calibrated image. This script sets up the model for either a left or
        right femur, as specified by the user. The analysis resamples the image to isotropic
        voxels, transforms the image, applies the bone mask and bins the data. It then
        creates the FE model for solving using FAIM (v8.1, Numerics Solutions Ltd, Calgary,
        Canada - Steven  Boyd).

        Input: Calibrated K2HPO4 Image (*.nii), Bone Mask (*_MASK.nii)

        Optional Parameters:
        1) Isotropic resample voxel size
        2) Left or Right Femur Bone Mask

        Output: N88 Model (*.n88model)

'''


    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoFEA-hip-builder",
        description=description
    )


    parser.add_argument("calibrated_image",
        help = "*_K2HPO4.nii image file")
    parser.add_argument("bone_mask",
        help = "*_MASK.nii mask image of bone")

    parser.add_argument("--mask_threshold", type = int,
                        default = 1,
                        help = "Set the threshold value to extract the bone of interest from the mask. (Default: %(default)s)")
    parser.add_argument("--iso_resolution", type = float, default = DEFAULT_FEMUR_ISO_RESOLUTION_MM,
                        help = "Set the isotropic voxel size [in mm]. (default: %(default)s [mm])")
    parser.add_argument("--mask_smoothing_spacing_threshold", type=float, default=DEFAULT_FEMUR_MASK_SMOOTHING_SPACING_THRESHOLD_MM,
                        help="Smooth resampled femur masks only when an input spacing dimension exceeds this value. (default: %(default)s [mm])")
    parser.add_argument("--femur_side", type = int, default = 1,
                        help = "Set whether the left or right femur is to be analyzed. 1 = Left; 2 = Right. (default: %(default)s)")
    parser.add_argument("--output_path", type=str, default=None,
                        help="Set output path for the N88 model file. (default: same as input image)")
    parser.add_argument("--poissons_ratio", type=float, default=0.3,
                        help="Sets the Poisson's ratio for the material(s) in the FE model. (default: %(default)s)")
    parser.add_argument("--elastic_E_func", type=str, default=None,
                        help="Function name for trabecular bone Young's modulus. (default: default_E)")
    parser.add_argument("--yield_comp_func", type=str, default=None,
                        help="Function name for trabecular compression yield. (default: none)")
    parser.add_argument("--yield_tens_func", type=str, default=None,
                        help="Function name for trabecular tension yield. (default: none)")
    parser.add_argument("--cort_elastic_E_func", type=str, default=None,
                        help="Function name for cortical bone Young's modulus. Defaults to --elastic_E_func.")
    parser.add_argument("--cort_yield_comp_func", type=str, default=None,
                        help="Function name for cortical compression yield. Defaults to --yield_comp_func.")
    parser.add_argument("--cort_yield_tens_func", type=str, default=None,
                        help="Function name for cortical tension yield. Defaults to --yield_tens_func.")
    parser.add_argument("--cort_poissons_ratio", type=float, default=None,
                        help="Poisson's ratio for cortical bone. Defaults to --poissons_ratio.")
    parser.add_argument("--pmma_E", type=float, default=2500,
                        help="Sets the Elastic Modulus for PMMA caps in the FE model. (default: %(default)s [MPa])")
    parser.add_argument("--pmma_v", type=float, default=0.3,
                        help="Sets the Poisson's ratio for the PMMA material(s) in the FE model. (default: %(default)s)")
    parser.add_argument("--pmma_thick", type=float, default=DEFAULT_PMMA_THICKNESS_MM,
                        help="Sets the fixed thickness for PMMA caps in the FE model. (default: %(default)s [mm])")
    parser.add_argument("--pmma_intrusion", type=float, default=DEFAULT_PMMA_INTRUSION_MM,
                        help="Sets how far anatomy can occupy the fixed PMMA fixture thickness before bone overlap is preserved during material combination. (default: %(default)s [mm])")
    parser.add_argument("--pmma_mat_id", type=int, default=5000,
                        help="Sets the material ID for the PMMA blocks. (default: %(default)s)")
    parser.add_argument("--fe_displacement", type=float, default=DEFAULT_FEMUR_FE_DISPLACEMENT,
                        help="Sets the applied displacement endpoint for the sideways-fall model. The default reports the force at 4%% displacement. (default: %(default)s)")
    parser.add_argument("--femur_shaft_length", type=float, default=DEFAULT_FEMUR_SHAFT_LENGTH_MM,
                        help="Retained proximal femur length [mm] for --femur_cut_mode fixed_length. (default: %(default)s [mm])")
    parser.add_argument("--femur_cut_mode", choices=["bbox_ratio", "proximal_box_ratio", "post_icp_flat_ratio", "lesser_trochanter", "fixed_length"],
                        default=DEFAULT_FEMUR_CUT_MODE,
                        help="Set the distal femur crop. post_icp_flat_ratio uses a fixed rough pre-ICP crop, then a flat aligned-frame final ratio crop after ICP. (default: %(default)s)")
    parser.add_argument("--femur_bbox_ratio", nargs=3, default=DEFAULT_FEMUR_BBOX_RATIO,
                        help="BBox-ratio crop in recipe order: reference constrained free. Use 'none' for a free axis. (default: %(default)s)")
    parser.add_argument("--femur_bbox_crop_from", nargs=3, default=DEFAULT_FEMUR_BBOX_CROP_FROM,
                        help="BBox crop side in recipe order: reference constrained free. Values are min, max, center/null. (default: %(default)s)")
    parser.add_argument("--femur_experimental_crop_ratio", type=float, default=DEFAULT_FEMUR_EXPERIMENTAL_RATIO,
                        help="Retained proximal femur length as a multiple of the proximal transverse max(x, y) width in proximal_box_ratio mode. (default: %(default)s)")
    parser.add_argument("--femur_proximal_reference_distance", type=float, default=DEFAULT_FEMUR_PROXIMAL_REFERENCE_DISTANCE_MM,
                        help="Proximal-most distance [mm] used to estimate transverse reference width in proximal_box_ratio mode. (default: %(default)s)")
    parser.add_argument("--femur_proximal_reference_width", choices=["max_xy", "rms_xy"], default=DEFAULT_FEMUR_PROXIMAL_REFERENCE_WIDTH,
                        help="Transverse width summary for proximal_box_ratio. max_xy uses max(x width, y width); rms_xy uses sqrt((x^2 + y^2) / 2). (default: %(default)s)")
    parser.add_argument("--femur_lesser_trochanter_distal_offset", type=float, default=DEFAULT_LESSER_TROCHANTER_DISTAL_OFFSET_MM,
                        help="Cut this many mm distal to the detected lesser trochanter in lesser_trochanter mode. (default: %(default)s [mm])")
    parser.add_argument("--femur_lesser_trochanter_distal_offset_percent", type=float, default=None,
                        help="Optional percentage of the detected greater-to-lesser trochanter z distance to cut distal to the lesser trochanter. Overrides --femur_lesser_trochanter_distal_offset when set. (default: %(default)s)")
    parser.add_argument("--femur_input_margin", type=float, default=DEFAULT_FEMUR_INPUT_MARGIN_MM,
                        help="Pad the input image/mask as needed so femur foreground has this margin before ICP. (default: %(default)s [mm])")
    parser.add_argument("--compartment_mask", type=str, default=None,
                        help="Optional trabecular/cortical compartment mask aligned with the bone mask. Defaults: cortical=1, trabecular=2.")
    parser.add_argument("--cortical_label", type=int, default=DEFAULT_CORTICAL_LABEL,
                        help="Label value for cortical bone in --compartment_mask. (default: %(default)s)")
    parser.add_argument("--trabecular_label", type=int, default=DEFAULT_TRABECULAR_LABEL,
                        help="Label value for trabecular bone in --compartment_mask. (default: %(default)s)")
    parser.add_argument("--reference_path", type=str, required=False, default=None,
                        help="Path to the reference vtk file for ICP registration. (default: None)")
    parser.add_argument("--pmma_yield_compression", type=float, default=None,
                        help="Sets the yield strength in compression for PMMA material in the FE model. (default: %(default)s [MPa])")
    parser.add_argument("--pmma_yield_tension", type=float, default=None,
                        help="Sets the yield strength in tension for PMMA material in the FE model. (default: %(default)s [MPa])")
    # Parse and display
    args = parser.parse_args()

    # Set default reference paths
    if args.reference_path is None:
        data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "dat")
        args.left_femur_reference = os.path.join(data_dir, "LT_FEMUR_SIDEWAYS_FALL_REF.vtk")
        args.right_femur_reference = os.path.join(data_dir, "RT_FEMUR_SIDEWAYS_FALL_REF.vtk")
    else:
        args.left_femur_reference = os.path.join(args.reference_path, "LT_FEMUR_SIDEWAYS_FALL_REF.vtk")
        args.right_femur_reference = os.path.join(args.reference_path, "RT_FEMUR_SIDEWAYS_FALL_REF.vtk")


    basename = remove_extension(os.path.basename(args.calibrated_image))

    if args.output_path is None:
        output_dir = os.path.dirname(args.calibrated_image)
    else:
        output_dir = args.output_path

    output_dir = os.path.abspath(output_dir)

    args.output_file = os.path.join(output_dir, f"{basename}.n88model")

    print(echo_arguments("ogoFEA-hip-builder", vars(args)))

    # Run program
    sidewaysFallFe(args)

if __name__ == '__main__':
    main()
