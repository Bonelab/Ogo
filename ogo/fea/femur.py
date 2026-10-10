"""Sideways-fall femur models: registration, GT-based shaft crop and supports.

Density uses cubic interpolation and labels use nearest-neighbor. Shared
geometry helpers live in alignment and boundary; the ogoFEA wrapper solves
and reports results. See docs/fea/implementation.md for the code map."""

import ogo.util.Helper as ogo
import os
import sys
import argparse
import json
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
from ogo.fea.image_io import write_vtk_image_with_sitk_geometry
from ogo.fea.hip_registration import estimate_femur_icp
from ogo.fea.shaft_geometry import (
    aligned_distal_scan_boundary,
    capture_distal_scan_face,
    measure_available_shaft,
    verify_model_shaft,
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
POST_ICP_DISTAL_SHAFT_SUPPORT_FRACTION = 0.9
DEFAULT_FEMUR_REFERENCE_MIN_SCALE = (0.8, 0.8, 0.75)
DEFAULT_FEMUR_REFERENCE_MAX_SCALE = (1.2, 1.2, 1.3)
DEFAULT_FEMUR_REGISTRATION_BACKEND = "numpy"
DEFAULT_FEMUR_REGISTRATION_LANDMARKS = 40000
DEFAULT_FEMUR_REGISTRATION_ITERATIONS = 50
DEFAULT_PMMA_THICKNESS_MM = 10.0
DEFAULT_PMMA_INTRUSION_MM = 6.0
DEFAULT_FEMUR_INPUT_MARGIN_MM = DEFAULT_PMMA_THICKNESS_MM + DEFAULT_PMMA_INTRUSION_MM
DEFAULT_FEMUR_SHAFT_LENGTH_MM = 120.0
DEFAULT_FEMUR_ROUGH_PRE_ICP_LENGTH_MM = 120.0
DEFAULT_FEMUR_CUT_MODE = "greater_trochanter_length"
DEFAULT_FEMUR_GREATER_TROCHANTER_DISTAL_LENGTH_MM = 10.0
DEFAULT_FEMUR_GREATER_TROCHANTER_INCLUSION_LENGTH_MM = 0.0
DEFAULT_FEMUR_PISTOIA_CRITICAL_VOLUME_PERCENT = 11.2
DEFAULT_FEMUR_PISTOIA_CRITICAL_STRAIN = 0.009
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
        "critical_volume": DEFAULT_FEMUR_PISTOIA_CRITICAL_VOLUME_PERCENT,
        "critical_strain": DEFAULT_FEMUR_PISTOIA_CRITICAL_STRAIN,
    }


def proximal_sideways_fall_fixture_plane(model_bounds, *, center_fraction):
    """Return a proximal-only sideways-fall PMMA contact plane.

    The shared bbox-relative helper scales y-projected fixtures over the full
    model z extent. That is unsafe for shaft-length sensitivity models because
    the proximal PMMA cap can also contact the distal shaft crop face. Keep the
    contact footprint proximal by using a fixed long-axis footprint anchored at
    the femoral-head end of the model.
    """
    plane = bbox_relative_fixture_plane(
        model_bounds,
        center_fraction=center_fraction,
        size_fraction=SIDEWAYS_FALL_FIXTURE_SIZE_FRACTION,
        projection_axis="y",
        shape=SIDEWAYS_FALL_FIXTURE_SHAPE,
    )
    z_min = float(model_bounds[4])
    z_max = float(model_bounds[5])
    z_span = z_max - z_min
    if z_span <= 0.0:
        raise ValueError("model_bounds must have positive z span.")
    long_axis_mm = min(float(FEMORAL_HEAD_FIXTURE_LONG_AXIS_EXTENSION_MM), z_span)
    x_width_mm = (
        float(model_bounds[1])
        - float(model_bounds[0])
        + float(FEMORAL_HEAD_FIXTURE_WIDTH_EXTENSION_MM)
    )
    center = list(plane["center"])
    center[2] = z_max - long_axis_mm / 2.0
    plane = dict(plane)
    plane["center"] = tuple(float(value) for value in center)
    plane["size"] = (float(long_axis_mm), float(x_width_mm))
    plane["footprint"] = "proximal_only"
    return plane


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


def pistoia_mask_output_path(output_file):
    """Return the model-space Pistoia ROI mask sidecar path for a model."""
    output_path = Path(output_file)
    return output_path.with_name(f"{output_path.stem}_pistoia_mask.nii.gz")


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


def vtk_matrix_from_rows(rows):
    """Return a VTK 4 x 4 matrix from JSON-serializable row values."""
    import vtk

    values = matrix4x4_to_numpy(rows)
    matrix = vtk.vtkMatrix4x4()
    matrix.Identity()
    for row in range(4):
        for col in range(4):
            matrix.SetElement(row, col, float(values[row, col]))
    return matrix


def load_icp_transform(path):
    """Read an ICP transform sidecar JSON file."""
    import json

    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    if data.get("type") != "ogo.femur.icp_transform":
        raise ValueError("ICP transform file has an unexpected type: %s" % data.get("type"))
    return data, vtk_matrix_from_rows(data["matrix"])


def write_icp_transform(path, *, matrix, icp_transform, reference_scale, femur_side, rough_crop):
    """Write an ICP transform sidecar JSON file for repeatable length sweeps."""
    import json

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    matrix_rows = matrix4x4_to_numpy(matrix).tolist()
    data = {
        "type": "ogo.femur.icp_transform",
        "version": 1,
        "femur_side": int(femur_side),
        "matrix": matrix_rows,
        "icp": {
            **{key: value for key, value in icp_transform.items()
               if key not in {"rotation", "translation"}},
            "iterations": int(icp_transform["iterations"]),
            "mean_distance": float(icp_transform["mean_distance"]),
        },
        "reference_scale": reference_scale,
        "rough_pre_icp_crop": rough_crop,
    }
    output.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return data


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


def crop_vtk_images_to_greater_trochanter_length(
    vtk_images,
    vtk_mask,
    *,
    gt_support_vtk,
    distal_boundary=None,
    retained_length_mm=DEFAULT_FEMUR_GREATER_TROCHANTER_DISTAL_LENGTH_MM,
    gt_inclusion_length_mm=DEFAULT_FEMUR_GREATER_TROCHANTER_INCLUSION_LENGTH_MM,
    labels=None,
):
    """Crop a supported model a fixed distance distal to its actual GT disk.

    The inputs must share the model grid and include the generated supports.
    Length is measured between voxel faces, which become mesh-node planes.
    No anatomical landmark or fixture-box surrogate is used for the origin.
    """
    import numpy as np

    from ogo.util.vtk_image import numpy_to_vtk_image, vtk_image_to_numpy

    retained_length_mm = float(retained_length_mm)
    if not np.isfinite(retained_length_mm) or retained_length_mm <= 0.0:
        raise ValueError("retained_length_mm must be positive.")
    if float(gt_inclusion_length_mm) != 0.0:
        raise ValueError("Shaft length starts at the GT disk edge; no additional offset is supported.")

    measurement = measure_available_shaft(
        vtk_mask, gt_support_vtk, requested_length_mm=retained_length_mm,
        distal_boundary=distal_boundary, labels=labels
    )
    active = _active_crop_mask(vtk_image_to_numpy(vtk_mask), labels)
    spacing = np.asarray(vtk_mask.GetSpacing(), dtype=np.float64)
    origin = np.asarray(vtk_mask.GetOrigin(), dtype=np.float64)
    gt_active = (vtk_image_to_numpy(gt_support_vtk) != 0) & active
    extent_z = vtk_mask.GetExtent()[4]
    gt_min_index = int(np.argwhere(gt_active)[:, 2].min())
    gt_support_distal_z = measurement["distal_length_origin_z"]
    available_mm = measurement["available_shaft_length_mm"]
    if not measurement["eligible_for_requested_length"]:
        raise ValueError(
            "Femur scan is too short for the requested greater-trochanter distal shaft length: "
            "requested %.4f mm, available %.4f mm below the GT support distal length origin."
            % (retained_length_mm, max(0.0, available_mm))
        )

    cut_index = int(round(gt_min_index - retained_length_mm / spacing[2]))
    cut_z = float(origin[2] + (extent_z + cut_index - 0.5) * spacing[2])
    keep = np.broadcast_to(np.arange(active.shape[2])[None, None, :] >= cut_index, active.shape)

    cropped_images, crop_face_image, meta = _crop_vtk_images_with_keep_mask(
        vtk_images,
        vtk_mask,
        active=active,
        keep=keep,
        meta={
            **measurement,
            "enabled": True,
            "method": "greater_trochanter_length",
            "retained_length_mm": float(gt_support_distal_z - cut_z),
            "requested_retained_length_mm": retained_length_mm,
            "available_below_gt_support_distal_edge_mm": float(available_mm),
            "greater_trochanter_support_distal_edge_z": float(gt_support_distal_z),
            "distal_length_origin_z": float(gt_support_distal_z),
            "cut_z_mm": float(cut_z),
            "reference_axis": "generated_gt_disk_distal_face_z",
            "status": "cropped",
            "crop_stage": "after ICP and support generation on full-scan model grid",
        },
    )
    # Select the external flat face even when coverage exactly equals the target
    # and no occupied voxel was removed by this crop.
    face_data = np.zeros(cropped_images[0].GetDimensions(), dtype=np.uint8)
    face_index = cut_index - meta["crop_slices_xyz"][2][0]
    if 0 <= face_index < face_data.shape[2]:
        cropped_active = active[tuple(slice(lo, hi) for lo, hi in meta["crop_slices_xyz"])]
        face_data[:, :, face_index] = cropped_active[:, :, face_index]
    crop_face_image = numpy_to_vtk_image(face_data, crop_face_image)
    meta["crop_face_voxels"] = int(np.count_nonzero(face_data))
    if not meta["crop_face_voxels"]:
        raise ValueError("The requested shaft cut has no bone contact surface.")
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
    extent = vtk_mask.GetExtent()
    extent_offset = (extent[0], extent[2], extent[4])
    origin = tuple(
        float(vtk_mask.GetOrigin()[axis]) + float(extent_offset[axis] + out_lo[axis]) * float(spacing[axis])
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


def straight_crop_face_support_surface_vtk(
    crop_face_vtk,
    material_vtk,
    *,
    support_fraction=POST_ICP_DISTAL_SHAFT_SUPPORT_FRACTION,
    output_value=1,
):
    """Return a central straight support patch on a flat distal crop face."""
    import numpy as np
    import vtk

    from ogo.util.vtk_image import numpy_to_vtk_image, vtk_image_to_numpy

    face = vtk_image_to_numpy(crop_face_vtk) != 0
    active = vtk_image_to_numpy(material_vtk) != 0
    out = np.zeros(face.shape, dtype=np.uint8)
    if not np.any(face) or not np.any(active):
        return numpy_to_vtk_image(out, material_vtk, vtk_array_type=vtk.VTK_UNSIGNED_CHAR)
    fraction = float(support_fraction)
    if not (0.0 < fraction <= 1.0):
        raise ValueError("support_fraction must be in (0, 1].")

    face = face & active
    if not np.any(face):
        return numpy_to_vtk_image(out, material_vtk, vtk_array_type=vtk.VTK_UNSIGNED_CHAR)

    coords = np.argwhere(face)
    lo = coords.min(axis=0)
    hi = coords.max(axis=0) + 1
    keep_lo = lo.copy()
    keep_hi = hi.copy()
    for axis in (0, 1):
        width = int(hi[axis] - lo[axis])
        target = max(1, min(width, int(np.floor(width * fraction))))
        margin = max(0, (width - target) // 2)
        keep_lo[axis] = int(lo[axis]) + margin
        keep_hi[axis] = int(keep_lo[axis]) + target

    keep = np.zeros(face.shape, dtype=bool)
    keep[
        int(keep_lo[0]) : int(keep_hi[0]),
        int(keep_lo[1]) : int(keep_hi[1]),
        int(lo[2]) : int(hi[2]),
    ] = True
    out[face & keep] = int(output_value)
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


def remove_extension(filename):
    """Remove all filename suffixes, including compound NIfTI extensions."""
    while True:
        filename, ext = os.path.splitext(filename)
        if not ext:
            break
    return filename


##
# Start script
def _prepare_femur_inputs(
    args,
):
    """Resample full inputs and prepare a registration-only proximal crop."""
    image = args.calibrated_image
    mask = args.bone_mask
    mask_threshold = args.mask_threshold
    iso_resolution = args.iso_resolution
    compartment_mask = args.compartment_mask
    pistoia_mask = args.pistoia_mask
    pistoia_mask_label = args.pistoia_mask_label or []
    femur_shaft_length = args.femur_shaft_length
    femur_input_margin = args.femur_input_margin

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
    native_distal_face = capture_distal_scan_face(maskThres)
    compartmentData = None
    if compartment_mask is not None:
        ogo.message("Reading trabecular/cortical compartment mask...")
        compartmentData = ogo.readNii(compartment_mask)
    pistoiaMaskData = None
    if pistoia_mask is not None:
        ogo.message("Reading Pistoia ROI mask...")
        pistoiaMaskData = ogo.readNii(pistoia_mask)
        if pistoia_mask_label:
            pistoiaMaskData = ogo.maskThreshold(pistoiaMaskData, pistoia_mask_label)

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
    if compartmentData is not None:
        compartmentData = resample_vtk_image_like_workflow(
            compartmentData,
            iso_resolution,
            interpolation="nearest",
        )
    if pistoiaMaskData is not None:
        pistoiaMaskData = resample_vtk_image_like_workflow(
            pistoiaMaskData,
            iso_resolution,
            interpolation="nearest",
        )

    ogo.message("Applying fixed-length rough femur crop after isotropic resampling and before ICP...")
    cropped_images, _rough_crop_face, bbox_crop_meta = crop_vtk_images_to_fixed_proximal_length(
        [imageData, maskThres],
        maskThres,
        retained_length_mm=femur_shaft_length,
        labels={1},
    )
    registration_mask = cropped_images[1]
    ogo.message(
        "rough pre-ICP crop slices xyz=%s; output shape xyz=%s; retained length=%8.4f; status=%s."
        % (
            bbox_crop_meta["crop_slices_xyz"],
            bbox_crop_meta["output_shape_xyz"],
            bbox_crop_meta["retained_length_mm"],
            bbox_crop_meta["status"],
        )
    )

    images_to_pad = [imageData, maskThres]
    if compartmentData is not None:
        images_to_pad.append(compartmentData)
    if pistoiaMaskData is not None:
        images_to_pad.append(pistoiaMaskData)
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
    if compartmentData is not None:
        compartmentData = padded_images[next_padded_index]
        next_padded_index += 1
    if pistoiaMaskData is not None:
        pistoiaMaskData = padded_images[next_padded_index]
    if any(padding["lower"]) or any(padding["upper"]):
        ogo.message(
            "Padded isotropic input image extent by lower=%s upper=%s voxels for FE transform safety."
            % (padding["lower"], padding["upper"])
        )
    else:
        ogo.message("Isotropic input image already has sufficient foreground safety margin.")

    return (
        imageData,
        maskThres,
        compartmentData,
        pistoiaMaskData,
        registration_mask,
        native_distal_face,
        input_spacing,
        bbox_crop_meta,
    )


def _register_femur_inputs(
    args,
    N88_fileName,
    imageData,
    maskThres,
    compartmentData,
    pistoiaMaskData,
    registration_mask,
    native_distal_face,
    input_spacing,
    bbox_crop_meta,
):
    """Estimate or reuse ICP, retaining full anatomy for the GT-relative final crop."""
    iso_resolution = args.iso_resolution
    femur_side = args.femur_side
    left_femur_reference = args.left_femur_reference
    right_femur_reference = args.right_femur_reference
    femur_icp_transform_in = args.femur_icp_transform_in
    femur_icp_transform_out = args.femur_icp_transform_out
    femur_input_margin = args.femur_input_margin
    mask_smoothing_spacing_threshold = args.mask_smoothing_spacing_threshold

    ##
    # Use the cropped, isotropic input geometry directly for ICP.
    image_rot = imageData
    mask_rot = maskThres
    compartment_rot = compartmentData
    pistoia_mask_rot = pistoiaMaskData

    ##
    # Align the input femur with the reference model, or reuse a saved transform.
    if femur_icp_transform_in:
        ogo.message("Loading fixed femur ICP transform...")
        try:
            icp_sidecar, icp = load_icp_transform(femur_icp_transform_in)
        except Exception as exc:
            ogo.message("Unable to load femur ICP transform: %s" % exc)
            sys.exit(1)
        if int(icp_sidecar.get("femur_side", femur_side)) != int(femur_side):
            ogo.message("Femur ICP transform side does not match requested femur side.")
            sys.exit(1)
        icp_transform = icp_sidecar.get("icp", {"iterations": 0, "mean_distance": 0.0})
        reference_scale = icp_sidecar.get("reference_scale", {})
        ogo.message(
            "Using fixed ICP reference-to-sample transform iterations=%d mean_distance=%0.4f"
            % (int(icp_transform.get("iterations", 0)), float(icp_transform.get("mean_distance", 0.0)))
        )
    else:
        ogo.message("Aligning input with reference model...")
        sample_surface_points = surface_points_from_vtk_mask(
            registration_mask,
            max_points=DEFAULT_FEMUR_REGISTRATION_LANDMARKS,
            sample_mode="linspace",
            sample_offset=0,
        )
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
            max_points=DEFAULT_FEMUR_REGISTRATION_LANDMARKS,
            reference_sample_mode="linspace",
            min_scale=DEFAULT_FEMUR_REFERENCE_MIN_SCALE,
            max_scale=DEFAULT_FEMUR_REFERENCE_MAX_SCALE,
        )
        ogo.message("Femur reference axis lengths: %s" % str(reference_scale["reference_axis_lengths"]))
        ogo.message("Sample femur axis lengths: %s" % str(reference_scale["sample_axis_lengths"]))
        ogo.message("Femur reference scale factors: %s" % str(reference_scale["scale_factors"]))

        icp_transform = estimate_femur_icp(
            moving_points=polydata_points(ref_poly),
            fixed_points=sample_surface_points,
            iterations=DEFAULT_FEMUR_REGISTRATION_ITERATIONS,
        )
        icp = point_transform_to_vtk_matrix(
            icp_transform["rotation"],
            icp_transform["translation"],
        )
        ogo.message(
            "ICP reference-to-sample iterations=%d mean_distance=%0.4f"
            % (icp_transform["iterations"], icp_transform["mean_distance"])
        )
        if not femur_icp_transform_out:
            femur_icp_transform_out = str(N88_fileName).replace(".n88model", "_icp.json")
        if femur_icp_transform_out:
            try:
                write_icp_transform(
                    femur_icp_transform_out,
                    matrix=icp,
                    icp_transform=icp_transform,
                    reference_scale=reference_scale,
                    femur_side=femur_side,
                    rough_crop=bbox_crop_meta,
                )
            except Exception as exc:
                ogo.message("Unable to write femur ICP transform: %s" % exc)
                sys.exit(1)
            ogo.message("Wrote femur ICP transform: %s" % femur_icp_transform_out)

    ogo.message("Applying the transformation and isotropic resampling to the image and mask...")
    output_origin, output_size = reference_grid_from_vtk_mask(
        mask_rot,
        icp,
        margin_voxels=int(round(femur_input_margin / max(iso_resolution, 1.0e-6))),
    )
    output_spacing = (float(iso_resolution), float(iso_resolution), float(iso_resolution))
    distal_scan_boundary = aligned_distal_scan_boundary(
        native_distal_face, icp, resampling_spacing=output_spacing,
        output_spacing=output_spacing,
    )
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
    pistoia_mask_trans = (
        transform_resample_vtk_image_to_reference_grid(
            pistoia_mask_rot,
            icp,
            output_origin=output_origin,
            output_size=output_size,
            output_spacing=output_spacing,
            interpolation="nearest",
        )
        if pistoia_mask_rot is not None
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

    return (
        image_trans,
        mask_trans,
        compartment_trans,
        pistoia_mask_trans,
        smooth_resampled_masks,
        distal_scan_boundary,
        icp_transform,
        reference_scale,
    )


def sidewaysFallFe(args):
    """Build a sideways-fall model and its shaft geometry and QC sidecars."""
    ogo.message("Start of Script...")

    ##
    # Collect the input arguments
    image = args.calibrated_image
    mask = args.bone_mask
    compartment_mask = args.compartment_mask
    pistoia_mask = args.pistoia_mask
    pistoia_mask_label = args.pistoia_mask_label or []

    iso_resolution = args.iso_resolution
    femur_side = args.femur_side
    poissons_ratio = args.poissons_ratio
    pmma_E = args.pmma_E
    pmma_v = args.pmma_v
    pmma_thick = args.pmma_thick
    pmma_intrusion = args.pmma_intrusion
    pmma_mat_id = args.pmma_mat_id
    fe_displacement = args.fe_displacement
    output_file = args.output_file
    pmma_yield_compression = args.pmma_yield_compression
    pmma_yield_tension = args.pmma_yield_tension
    femur_shaft_length = args.femur_shaft_length
    femur_greater_trochanter_distal_length = args.femur_greater_trochanter_distal_length
    femur_greater_trochanter_inclusion_length = args.femur_greater_trochanter_inclusion_length
    femur_input_margin = args.femur_input_margin
    femur_icp_transform_in = args.femur_icp_transform_in
    femur_icp_transform_out = args.femur_icp_transform_out
    cortical_label = args.cortical_label
    trabecular_label = args.trabecular_label


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
    ogo.message("Femur Crop Recipe: registration-only rough crop, then a complete flat shaft distal to the GT disk edge")
    if femur_icp_transform_in:
        ogo.message("Femur ICP Transform In: %s" % femur_icp_transform_in)
    if femur_icp_transform_out:
        ogo.message("Femur ICP Transform Out: %s" % femur_icp_transform_out)
    if args.femur_cut_mode != DEFAULT_FEMUR_CUT_MODE:
        raise ValueError("Only the GT-disk-relative complete flat shaft recipe is supported.")
    if femur_greater_trochanter_inclusion_length != 0:
        raise ValueError("The shaft length origin must be the GT disk distal edge, without an extra offset.")
    ogo.message("Femur Rough Pre-ICP Retained Length [mm]: %8.4f" % femur_shaft_length)
    ogo.message("Femur Shaft Length Distal To GT Disk [mm]: %8.4f" % femur_greater_trochanter_distal_length)
    if compartment_mask is not None:
        ogo.message("Compartment Mask: %s" % compartment_mask)
        ogo.message("Cortical Label: %d" % cortical_label)
        ogo.message("Trabecular Label: %d" % trabecular_label)
    if pistoia_mask is not None:
        ogo.message("Pistoia ROI Mask: %s" % pistoia_mask)
        if pistoia_mask_label:
            ogo.message("Pistoia ROI Label(s): %s" % pistoia_mask_label)
    ogo.message("Input Foreground Safety Margin [mm]: %8.4f" % femur_input_margin)

    (
        imageData,
        maskThres,
        compartmentData,
        pistoiaMaskData,
        registration_mask,
        native_distal_face,
        input_spacing,
        bbox_crop_meta,
    ) = _prepare_femur_inputs(
        args,
    )

    (
        image_trans,
        mask_trans,
        compartment_trans,
        pistoia_mask_trans,
        smooth_resampled_masks,
        distal_scan_boundary,
        icp_transform,
        reference_scale,
    ) = _register_femur_inputs(
        args,
        N88_fileName,
        imageData,
        maskThres,
        compartmentData,
        pistoiaMaskData,
        registration_mask,
        native_distal_face,
        input_spacing,
        bbox_crop_meta,
    )

    # Keep the aligned anatomy intact until the actual support disks exist.


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
    pistoia_mask_on_model_grid = (
        ogo.prepareFiniteElementImage(ogo.cast2unsignchar(pistoia_mask_trans))
        if pistoia_mask_trans is not None
        else None
    )


    # Apply the mask to the bone
    ogo.message("Applying the bone mask to the image...")
    bone_image = ogo.applyMask(cast_image, cast_mask)

    # Ensure connectivity
    ogo.message("Performing Image connectivity...")
    conn = ogo.imageConnectivity(bone_image)


    change = ogo.prepareFiniteElementImage(conn)
    femur_mask_on_model_grid = ogo.prepareFiniteElementImage(cast_mask)
    distal_crop_face_change = None

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
    femoral_head_plane = proximal_sideways_fall_fixture_plane(
        model_bounds,
        center_fraction=FEMORAL_HEAD_FIXTURE_CENTER_FRACTION,
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
    greater_trochanter_plane = proximal_sideways_fall_fixture_plane(
        model_bounds,
        center_fraction=GREATER_TROCHANTER_FIXTURE_CENTER_FRACTION,
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
        if pistoia_mask_on_model_grid is not None:
            images_to_fit.append(pistoia_mask_on_model_grid)
        fitted_images, contact_fit = fit_vtk_images_to_physical_bounds(
            images_to_fit,
            desired_bounds=desired_bounds,
            constants=[0] * len(images_to_fit),
        )
        change = fitted_images[0]
        femur_mask_on_model_grid = fitted_images[1]
        next_fit_index = 2
        if pistoia_mask_on_model_grid is not None:
            pistoia_mask_on_model_grid = fitted_images[next_fit_index]
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
        keep_largest_component=True,
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
        keep_largest_component=True,
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

    # Freeze proximal supports before shortening the shaft. Measuring their
    # generated voxel faces avoids mistaking the larger fixture box for a disk.
    greaterTrochanterPMMA = ogo.applyMask(greaterTrochanterPMMA, ogo.cast2unsignchar(conn2))
    femoralHeadPMMA = ogo.applyMask(femoralHeadPMMA, ogo.cast2unsignchar(conn2))
    images_to_crop = [conn2, change, femur_mask_on_model_grid, femoralHeadPMMA, greaterTrochanterPMMA]
    if pistoia_mask_on_model_grid is not None:
        images_to_crop.append(pistoia_mask_on_model_grid)
    shaft_crop = {}
    try:
        shaft_crop = measure_available_shaft(
            conn2, greaterTrochanterPMMA,
            requested_length_mm=femur_greater_trochanter_distal_length,
            distal_boundary=distal_scan_boundary,
        )
        shaft_crop.update(pre_icp_crop=bbox_crop_meta,
                          proximal_supports_frozen_before_crop=True,
                          registration={key: value for key, value in icp_transform.items()
                                        if key not in {"rotation", "translation"}},
                          reference_scale=reference_scale)
        with open(str(N88_fileName).replace(".n88model", "_shaft_geometry.json"), "w") as stream:
            json.dump(shaft_crop, stream, indent=2)
        cropped_images, distal_crop_face_change, shaft_crop = crop_vtk_images_to_greater_trochanter_length(
            images_to_crop,
            conn2,
            gt_support_vtk=greaterTrochanterPMMA,
            distal_boundary=distal_scan_boundary,
            retained_length_mm=femur_greater_trochanter_distal_length,
            gt_inclusion_length_mm=femur_greater_trochanter_inclusion_length,
        )
    except ValueError as exc:
        # Preserve a visual explanation of coverage failure without inventing
        # a distal boundary or sending an ineligible model to the solver.
        if shaft_crop.get("status") == "too_short":
            from ogo.fea.qc_render import try_export_model_qc

            available_model = ogo.applyTestBase(ogo.Image2Mesh(conn2), material_table)
            from ogo.fea.validation import measure_model, write_measurements

            write_measurements(N88_fileName, measure_model(
                available_model, "hip", shaft=shaft_crop, ineligible=True))
            shaft_crop["exports"] = {"qc_3d": try_export_model_qc(
                N88_fileName, "hip", model=available_model, ineligible=True)}
            with open(str(N88_fileName).replace(".n88model", "_shaft_geometry.json"), "w") as stream:
                json.dump(shaft_crop, stream, indent=2)
        ogo.message(str(exc))
        sys.exit(1)
    conn2, change, femur_mask_on_model_grid, femoralHeadPMMA, greaterTrochanterPMMA = cropped_images[:5]
    if pistoia_mask_on_model_grid is not None:
        pistoia_mask_on_model_grid = cropped_images[5]
    retained_length_mm = shaft_crop["retained_length_mm"]
    shaft_crop["pre_icp_crop"] = bbox_crop_meta
    shaft_crop["proximal_supports_frozen_before_crop"] = True
    shaft_crop["registration"] = {
        key: value for key, value in icp_transform.items()
        if key not in {"rotation", "translation"}
    }
    shaft_crop["reference_scale"] = reference_scale
    ogo.message(
        "Generated GT-disk-edge crop: shaft=%.4f mm, GT distal face z=%.4f, "
        "shaft cut face z=%.4f, available=%.4f mm."
        % (retained_length_mm, shaft_crop["distal_length_origin_z"],
           shaft_crop["cut_z_mm"], shaft_crop["available_below_gt_support_distal_edge_mm"])
    )
    with open(str(N88_fileName).replace(".n88model", "_shaft_geometry.json"), "w") as stream:
        json.dump(shaft_crop, stream, indent=2)

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
    # Distal Femur (df) support nodes on the final flat post-ICP crop face.
    ogo.message("Determining distal femur nodes...")
    if distal_crop_face_change is None:
        ogo.message("Final femur crop requires a distal crop-face mask.")
        sys.exit(1)
    distal_support_direction = (0.0, 0.0, -1.0)
    ogo.message(
        "Distal Femur post-ICP straight support patch: fraction=%8.4f direction=%s"
        % (POST_ICP_DISTAL_SHAFT_SUPPORT_FRACTION, distal_support_direction)
    )
    distal_surface = straight_crop_face_support_surface_vtk(
        distal_crop_face_change,
        change,
        support_fraction=POST_ICP_DISTAL_SHAFT_SUPPORT_FRACTION,
        output_value=1,
    )
    df_visible_node_IDS = directional_face_node_ids_from_voxel_mask(
        model2,
        distal_surface,
        name=DISTAL_FEMUR_NODE_SET,
        direction=distal_support_direction,
    )
    if df_visible_node_IDS.GetNumberOfTuples() == 0:
        ogo.message("No distal femur nodes found on the post-ICP straight shaft support patch.")
        sys.exit(1)

    ogo.message("-- found %d distal shaft outer-face nodes."
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

    if pistoia_mask_on_model_grid is not None:
        mask_sidecar = pistoia_mask_output_path(N88_fileName)
        ogo.message("Writing model-space Pistoia ROI mask: %s" % mask_sidecar)
        write_vtk_image_with_sitk_geometry(
            ogo.cast2unsignchar(pistoia_mask_on_model_grid),
            mask_sidecar,
        )

    ##
    # Write out n88model file
    ogo.message("Writing out n88model file: %s" % N88_fileName)
    shaft_crop = verify_model_shaft(model2, shaft_crop)
    from ogo.fea.validation import measure_model, write_measurements

    write_measurements(N88_fileName, measure_model(
        model2, "hip", shaft=shaft_crop, registration_rotation=icp_transform["rotation"]))
    write_model(model2, N88_fileName)
    saved_model = vtkbone.vtkboneN88ModelReader()
    saved_model.SetFileName(str(N88_fileName))
    saved_model.Update()
    shaft_crop = verify_model_shaft(saved_model.GetOutput(), shaft_crop)

    from ogo.fea.qc_render import try_export_model_qc

    shaft_crop["exports"] = {"qc_3d": try_export_model_qc(N88_fileName, "hip", model=model2)}
    with open(str(N88_fileName).replace(".n88model", "_shaft_geometry.json"), "w") as stream:
        json.dump(shaft_crop, stream, indent=2)

    ##
    ogo.message("Done writing n88model.")

    ##
    # End of script
    ogo.message("End of Script.")
    sys.exit()


def build_parser():
    """Define the femur builder's arguments and defaults."""
    description = """Build a sideways-fall femur model from calibrated K2HPO4 density\n    and a femur label mask. Estimate ICP from a proximal crop, then retain a\n    complete flat shaft measured from the distal GT support edge. Use ogoFEA\n    hip to also solve and report results."""


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
                        help="Registration-only proximal crop length [mm]; the solved model uses the full scan before its final GT-relative crop. (default: %(default)s)")
    parser.add_argument("--femur_cut_mode", choices=[DEFAULT_FEMUR_CUT_MODE],
                        default=DEFAULT_FEMUR_CUT_MODE,
                        help=argparse.SUPPRESS)
    parser.add_argument("--femur_greater_trochanter_distal_length", type=float, default=DEFAULT_FEMUR_GREATER_TROCHANTER_DISTAL_LENGTH_MM,
                        help="Retained complete flat shaft length [mm] from the actual GT disk distal edge. Insufficient coverage fails generation. (default: %(default)s)")
    parser.add_argument("--femur_greater_trochanter_inclusion_length", type=float, default=DEFAULT_FEMUR_GREATER_TROCHANTER_INCLUSION_LENGTH_MM,
                        help=argparse.SUPPRESS)
    parser.add_argument("--femur_input_margin", type=float, default=DEFAULT_FEMUR_INPUT_MARGIN_MM,
                        help="Pad the input image/mask as needed so femur foreground has this margin before ICP. (default: %(default)s [mm])")
    parser.add_argument("--femur_icp_transform_in", type=str, default=None,
                        help="Optional femur ICP transform JSON to reuse instead of estimating ICP. Intended for fixed-transform length sweeps. (default: %(default)s)")
    parser.add_argument("--femur_icp_transform_out", type=str, default=None,
                        help="Optional path to write the estimated femur ICP transform JSON. (default: %(default)s)")
    parser.add_argument("--compartment_mask", type=str, default=None,
                        help="Optional trabecular/cortical compartment mask aligned with the bone mask. Defaults: cortical=1, trabecular=2.")
    parser.add_argument("--pistoia_mask", type=str, default=None,
                        help="Optional ROI mask aligned with the input image for model-space masked Pistoia reporting.")
    parser.add_argument("--pistoia_mask_label", type=int, action="append", default=None,
                        help="Label to keep from --pistoia_mask before masked Pistoia. Repeat for a multi-label ROI.")
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
    return parser


def main(argv=None):
    """Generate a femur model from explicit arguments or the command line."""
    args = build_parser().parse_args(argv)

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
