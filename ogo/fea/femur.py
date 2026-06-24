"""Femur-specific defaults and helpers for sideways-fall FE generation.

Default femur workflow:
1. Read the calibrated density image and whole-femur mask. The public wrapper
   chooses left/right from ``--side`` and the lower-level script uses
   ``femur_side=1`` for left and ``2`` for right.
2. Pre-rotate the femur by side to provide a stable starting orientation for
   ICP. This is a VTK reslice in the native input grid.
3. Run ICP to the bundled side-specific femur reference, then apply the ICP
   transform and set the final output spacing to 1.0 x 1.0 x 1.0 mm in one
   shared VTK reslice helper. Image data use cubic interpolation; bone and
   compartment labels use nearest-neighbor interpolation.
4. Smooth the transformed femur mask with one binary close/open pass only when
   at least one input spacing dimension is coarser than 2 mm. If a compartment
   mask is supplied, the derived cortical binary mask follows the same rule.
5. Standardize the distal shaft by cutting on a flat model-grid z plane. The
   default cut mode detects the lesser-trochanter cross-section peak and keeps
   the femur to 50 mm distal to that landmark. If the scan does not include the
   required distal field of view, model generation fails instead of silently
   using a shorter shaft. ``fixed_length`` mode is available for debugging.
6. Generate two geometric PMMA fixtures: a round femoral-head loading fixture
   and a box-like greater-trochanter contact fixture. Defaults are 6 mm PMMA
   thickness and 6 mm PMMA intrusion toward bone before material merging. The
   femoral-head footprint is widened by 10 mm and lengthened by 80 mm so the
   anatomical cropping determines the actual contact.
7. Apply sideways-fall boundary conditions: prescribed displacement at the
   femoral-head PMMA cap toward the greater trochanter, loading-direction
   constraint at the greater-trochanter PMMA cap, and distal shaft constraints
   to remove rigid-body motion.
8. Build materials with the same shared bone/PMMA material-table helper used by
   the spine workflow. If no compartment mask is supplied, femur bone is one
   trabecular-style region. If supplied, cortical=1 and trabecular=2 by default.
"""


LEFT_FEMUR = 1
RIGHT_FEMUR = 2

DEFAULT_FEMUR_ISO_RESOLUTION_MM = 1.0
DEFAULT_FEMUR_FE_DISPLACEMENT = -4.0
DEFAULT_FEMUR_TARGET_DISPLACEMENT_PERCENT = 4.0
DEFAULT_FEMUR_MASK_SMOOTHING_SPACING_THRESHOLD_MM = 2.0
FEMORAL_HEAD_FIXTURE_WIDTH_EXTENSION_MM = 10.0
FEMORAL_HEAD_FIXTURE_LONG_AXIS_EXTENSION_MM = 80.0
DEFAULT_PMMA_THICKNESS_MM = 6.0
DEFAULT_PMMA_INTRUSION_MM = 6.0
DEFAULT_FEMUR_SHAFT_LENGTH_MM = 100.0
DEFAULT_FEMUR_CUT_MODE = "lesser_trochanter"
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
    """Return the legacy output suffix for a femur side."""
    if femur_side == LEFT_FEMUR:
        return "_LT_FEMUR_SF"
    if femur_side == RIGHT_FEMUR:
        return "_RT_FEMUR_SF"
    raise ValueError("femur_side must be 1 for left or 2 for right.")


def side_rotation(femur_side):
    """Return the pre-alignment z rotation for a femur side."""
    if femur_side == LEFT_FEMUR:
        return 90
    if femur_side == RIGHT_FEMUR:
        return -90
    raise ValueError("femur_side must be 1 for left or 2 for right.")


def sideways_fall_output_name(output_file, femur_side):
    """Apply the legacy left/right suffix to a sideways-fall output path."""
    suffix = side_suffix(femur_side)
    if str(output_file).endswith(".n88model"):
        return str(output_file).replace(".n88model", f"{suffix}.n88model")
    return f"{output_file}{suffix}.n88model"


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


def detect_lesser_trochanter_cut_z(
    vtk_mask,
    *,
    distal_offset_mm=0.0,
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

    lesser_index = int(distal_indices[np.argmax(area[distal_indices])])
    lesser_z = float(z[lesser_index])
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
):
    """Apply a flat distal crop to aligned femur image and mask.

    In the default mode, the cut plane is detected from the lesser-trochanter
    cross-section profile. A fixed retained length can be requested explicitly
    for debugging or legacy comparisons.
    """
    cut_mode = str(cut_mode).lower()
    if cut_mode == "lesser_trochanter":
        meta = detect_lesser_trochanter_cut_z(
            mask_vtk,
            distal_offset_mm=lesser_trochanter_distal_offset_mm,
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
