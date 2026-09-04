"""Shared boundary-condition geometry helpers for FE workflows."""

from ogo.util.vtk_image import numpy_to_vtk_image
from ogo.util.vtk_image import vtk_image_to_numpy


def axis_index(axis):
    """Return the array axis index for a named anatomical/image axis."""
    axis_map = {"x": 0, "y": 1, "z": 2}
    try:
        return axis_map[str(axis).lower()]
    except KeyError as exc:
        raise ValueError("axis must be 'x', 'y', or 'z'.") from exc


def _axis_bounds(bounds, axis):
    index = 2 * int(axis)
    return float(bounds[index]), float(bounds[index + 1])


def _axis_extent_from_vector(bounds, vector):
    import numpy as np

    spans = np.asarray(
        [
            float(bounds[1]) - float(bounds[0]),
            float(bounds[3]) - float(bounds[2]),
            float(bounds[5]) - float(bounds[4]),
        ],
        dtype=float,
    )
    return float(np.sum(np.abs(np.asarray(vector, dtype=float)) * spans))


def _contact_plane_axes(projection_axis, normal_sign):
    """Return stable in-plane axes for a coordinate-axis contact plane."""
    axis = axis_index(projection_axis)
    sign = 1.0 if float(normal_sign) >= 0.0 else -1.0
    # Workflow editor markups serialize nearly-cardinal RAS axes with tiny
    # floating-point components. The projected-disk algorithm buckets columns in
    # plane coordinates, so preserving that authored basis avoids off-by-one
    # footprint differences at the contact edge.
    ras_quarter_turn = 6.123233995736766e-17
    if axis == 0:
        return (0.0, 0.0, sign), (0.0, -1.0, 0.0)
    if axis == 1:
        if sign > 0.0:
            return (
                1.208164965010145e-16,
                ras_quarter_turn,
                -1.0,
            ), (
                -1.0,
                ras_quarter_turn,
                -1.208164965010145e-16,
            )
        return (
            9.331099404642518e-17,
            ras_quarter_turn,
            1.0,
        ), (
            -1.0,
            -ras_quarter_turn,
            9.331099404642518e-17,
        )
    return (1.0, 0.0, 0.0), (0.0, sign, 0.0)


def bbox_relative_contact_bounds(
    model_bounds,
    *,
    center_fraction,
    size_fraction,
    projection_axis="y",
    shape="rectangle",
):
    """Return physical contact ROI bounds from bbox-relative plane values.

    ``center_fraction`` is authored in x/y/z order against voxel-center bounds.
    The two ``size_fraction`` values scale the lateral axes of the projection
    plane in coordinate order. The projection axis itself spans the full input
    bounds so downstream geometry can find the nearest supported anatomy.
    """
    if len(model_bounds) != 6:
        raise ValueError("model_bounds must contain x/y/z min/max values.")
    if len(center_fraction) != 3:
        raise ValueError("center_fraction must contain x/y/z fractions.")
    if len(size_fraction) != 2:
        raise ValueError("size_fraction must contain the two in-plane scale factors.")

    projection_index = axis_index(projection_axis)
    out = [float(value) for value in model_bounds]
    lateral_axes = [axis for axis in (0, 1, 2) if axis != projection_index]
    lateral_lengths = []
    for size_value, axis in zip(size_fraction, lateral_axes):
        lo, hi = _axis_bounds(model_bounds, axis)
        span = hi - lo
        if span <= 0.0:
            raise ValueError("model_bounds must have positive span on every lateral axis.")
        length = span * float(size_value)
        if length <= 0.0:
            raise ValueError("size_fraction values must be positive.")
        lateral_lengths.append(length)

    if str(shape).strip().lower() == "square":
        lateral_lengths = [min(lateral_lengths)] * len(lateral_lengths)

    for length, axis in zip(lateral_lengths, lateral_axes):
        lo, hi = _axis_bounds(model_bounds, axis)
        span = hi - lo
        center = lo + float(center_fraction[axis]) * span
        out[2 * axis] = center - length / 2.0
        out[2 * axis + 1] = center + length / 2.0

    lo, hi = _axis_bounds(model_bounds, projection_index)
    out[2 * projection_index] = lo
    out[2 * projection_index + 1] = hi
    return tuple(out)


def bbox_relative_contact_direction(center_fraction, *, projection_axis="y"):
    """Return the cap side implied by a bbox-relative plane fraction."""
    if len(center_fraction) != 3:
        raise ValueError("center_fraction must contain x/y/z fractions.")
    return "up" if float(center_fraction[axis_index(projection_axis)]) >= 0.5 else "down"


def bbox_relative_contact_plane(
    model_bounds,
    *,
    center_fraction,
    size_fraction,
    projection_axis="y",
    shape="square",
):
    """Return a physical contact-plane definition from bbox-relative values."""
    if len(model_bounds) != 6:
        raise ValueError("model_bounds must contain x/y/z min/max values.")
    if len(center_fraction) != 3:
        raise ValueError("center_fraction must contain x/y/z fractions.")
    if len(size_fraction) != 2:
        raise ValueError("size_fraction must contain the two in-plane scale factors.")

    projection_index = axis_index(projection_axis)
    center = []
    for axis in range(3):
        lo, hi = _axis_bounds(model_bounds, axis)
        center.append(lo + float(center_fraction[axis]) * (hi - lo))

    normal_sign = -1.0 if float(center_fraction[projection_index]) >= 0.5 else 1.0
    normal = [0.0, 0.0, 0.0]
    normal[projection_index] = normal_sign
    u_axis, v_axis = _contact_plane_axes(projection_axis, normal_sign)
    size = (
        _axis_extent_from_vector(model_bounds, u_axis) * float(size_fraction[0]),
        _axis_extent_from_vector(model_bounds, v_axis) * float(size_fraction[1]),
    )
    return {
        "center": tuple(float(value) for value in center),
        "normal": tuple(float(value) for value in normal),
        "u_axis": tuple(float(value) for value in u_axis),
        "v_axis": tuple(float(value) for value in v_axis),
        "size": tuple(float(value) for value in size),
        "shape": str(shape).strip().lower(),
    }


def bounds_with_reference_extent(current_bounds, reference_bounds):
    """Return bounds centered on ``current_bounds`` with ``reference_bounds`` extents."""
    import numpy as np

    current = np.asarray(current_bounds, dtype=float)
    reference = np.asarray(reference_bounds, dtype=float)
    if current.shape != (6,) or reference.shape != (6,):
        raise ValueError("bounds must contain x/y/z min/max values.")
    out = []
    for axis in range(3):
        current_lo, current_hi = current[2 * axis], current[2 * axis + 1]
        reference_lo, reference_hi = reference[2 * axis], reference[2 * axis + 1]
        center = (current_lo + current_hi) / 2.0
        extent = reference_hi - reference_lo
        out.extend([center - extent / 2.0, center + extent / 2.0])
    return tuple(float(value) for value in out)


def foreground_voxel_center_bounds_from_mask(mask, *, origin, spacing):
    """Return x/y/z physical bounds of nonzero voxel centers."""
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
    return foreground_voxel_center_bounds_from_mask(
        vtk_image_to_numpy(vtk_image),
        origin=vtk_image.GetOrigin(),
        spacing=vtk_image.GetSpacing(),
    )


def largest_connected_component(mask):
    """Keep the largest connected component in a binary NumPy mask."""
    import numpy as np
    from scipy.ndimage import label

    if not np.any(mask):
        return mask.astype(bool)
    labels, _ = label(mask)
    sizes = np.bincount(labels.ravel())
    sizes[0] = 0
    return labels == int(np.argmax(sizes))


def should_smooth_resampled_mask(spacing, threshold_mm):
    """Return True when input spacing is coarse enough to justify mask smoothing."""
    return any(float(value) > float(threshold_mm) for value in spacing)


def smooth_binary_mask_vtk(mask_vtk, *, close_iter=1, open_iter=1):
    """Apply light binary closing/opening to a VTK mask and keep the main component."""
    import numpy as np
    import vtk
    from scipy.ndimage import binary_closing, binary_dilation, binary_erosion

    mask = vtk_image_to_numpy(mask_vtk, processing_order=True) > 0
    if close_iter > 0:
        mask = binary_closing(mask, structure=np.ones((3, 3, 3), dtype=bool), iterations=int(close_iter))
    if open_iter > 0:
        mask = binary_erosion(mask, iterations=int(open_iter))
        mask = binary_dilation(mask, iterations=int(open_iter))
    mask = largest_connected_component(mask)
    return numpy_to_vtk_image(
        mask.astype(np.uint8),
        mask_vtk,
        vtk.VTK_UNSIGNED_CHAR,
        processing_order=True,
    )


def smooth_label_mask_vtk(label_vtk, *, sigma_mm=1.0, threshold=0.5):
    """Smooth a multi-label mask by smoothing each label probability separately."""
    import numpy as np
    import SimpleITK as sitk
    import vtk

    labels_xyz = np.asarray(vtk_image_to_numpy(label_vtk))
    unique_labels = [int(value) for value in np.unique(labels_xyz) if int(value) != 0]
    if not unique_labels:
        return numpy_to_vtk_image(
            labels_xyz.astype(np.uint8, copy=False),
            label_vtk,
            vtk.VTK_UNSIGNED_CHAR,
        )

    labels_zyx = np.transpose(labels_xyz, (2, 1, 0))
    scores = []
    for label_value in unique_labels:
        image = sitk.GetImageFromArray((labels_zyx == label_value).astype(np.float32))
        image.SetSpacing(tuple(float(value) for value in label_vtk.GetSpacing()))
        smoothed = sitk.SmoothingRecursiveGaussian(image, float(sigma_mm))
        scores.append(sitk.GetArrayFromImage(smoothed))
    stacked = np.stack(scores, axis=0)
    best_index = np.argmax(stacked, axis=0)
    best_score = np.take_along_axis(
        stacked,
        best_index[np.newaxis, ...],
        axis=0,
    )[0]
    out_zyx = np.zeros(labels_zyx.shape, dtype=labels_xyz.dtype)
    label_values = np.asarray(unique_labels, dtype=labels_xyz.dtype)
    active = best_score >= float(threshold)
    out_zyx[active] = label_values[best_index[active]]
    out_xyz = np.transpose(out_zyx, (2, 1, 0))
    return numpy_to_vtk_image(
        out_xyz.astype(np.uint8, copy=False),
        label_vtk,
        vtk.VTK_UNSIGNED_CHAR,
    )


def resample_vtk_image_to_spacing(vtk_image, target_spacing_mm, *, interpolation="nearest"):
    """Resample a VTK image to isotropic spacing while preserving physical origin."""
    import numpy as np
    import SimpleITK as sitk
    import vtk
    from vtk.util.numpy_support import numpy_to_vtk

    target = float(target_spacing_mm)
    target_spacing = (target, target, target)
    array_xyz = vtk_image_to_numpy(vtk_image)
    image = sitk.GetImageFromArray(np.transpose(array_xyz, (2, 1, 0)))
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

    out_xyz = np.transpose(sitk.GetArrayFromImage(resampler.Execute(image)), (2, 1, 0))
    vtk_type = vtk.VTK_UNSIGNED_CHAR if interpolation == "nearest" else vtk_image.GetScalarType()
    scalars = numpy_to_vtk(out_xyz.ravel(order="F"), deep=True, array_type=vtk_type)
    out = vtk.vtkImageData()
    out.SetDimensions(out_xyz.shape)
    out.SetOrigin(*vtk_image.GetOrigin())
    out.SetSpacing(*target_spacing)
    out.GetPointData().SetScalars(scalars)
    return out


def _bounding_box(mask):
    import numpy as np

    coords = np.array(np.where(mask))
    if coords.size == 0:
        raise ValueError("Cannot generate a cap from an empty mask.")
    min_coords = coords.min(axis=1)
    max_coords = coords.max(axis=1) + 1
    return tuple(slice(min_coords[i], max_coords[i]) for i in range(mask.ndim))


def _expanded_bbox(mask, extrusion_axis, direction, thickness):
    bb = _bounding_box(mask)
    bb_list = list(bb)
    if direction == "up":
        bb_list[extrusion_axis] = slice(
            bb[extrusion_axis].start,
            min(bb[extrusion_axis].stop + thickness, mask.shape[extrusion_axis]),
        )
    elif direction == "down":
        bb_list[extrusion_axis] = slice(
            max(bb[extrusion_axis].start - thickness, 0),
            bb[extrusion_axis].stop,
        )
    else:
        raise ValueError("direction must be 'up' or 'down'.")
    return tuple(bb_list), bb


def _contact_surface(mask, extrusion_axis, direction):
    import numpy as np

    surface = np.zeros_like(mask, dtype=bool)
    dims = mask.shape
    if extrusion_axis == 0:
        for j in range(dims[1]):
            for k in range(dims[2]):
                hit = np.flatnonzero(mask[:, j, k])
                if hit.size:
                    surface[int(hit.max() if direction == "up" else hit.min()), j, k] = True
    elif extrusion_axis == 1:
        for i in range(dims[0]):
            for k in range(dims[2]):
                hit = np.flatnonzero(mask[i, :, k])
                if hit.size:
                    surface[i, int(hit.max() if direction == "up" else hit.min()), k] = True
    else:
        for i in range(dims[0]):
            for j in range(dims[1]):
                hit = np.flatnonzero(mask[i, j, :])
                if hit.size:
                    surface[i, j, int(hit.max() if direction == "up" else hit.min())] = True
    return surface


def _shift_surface(surface, extrusion_axis, direction, thickness):
    import numpy as np

    cap = np.zeros_like(surface, dtype=bool)
    axis_len = int(surface.shape[extrusion_axis])
    step_sign = 1 if direction == "up" else -1
    for step in range(1, max(int(thickness), 0) + 1):
        shifted = np.zeros_like(surface, dtype=bool)
        if step_sign > 0:
            src = slice(0, axis_len - step)
            dst = slice(step, axis_len)
        else:
            src = slice(step, axis_len)
            dst = slice(0, axis_len - step)
        if extrusion_axis == 0:
            shifted[dst, :, :] = surface[src, :, :]
        elif extrusion_axis == 1:
            shifted[:, dst, :] = surface[:, src, :]
        else:
            shifted[:, :, dst] = surface[:, :, src]
        cap |= shifted
    return cap


def _fit_cap(mask, cap, extrusion_axis, direction, thickness):
    import numpy as np

    flat_cap = np.zeros_like(cap, dtype=bool)
    cap = largest_connected_component(cap)
    if not np.any(cap):
        return flat_cap

    thickness_records = []
    if extrusion_axis == 0:
        outer = int(np.max(np.where(cap)[0])) if direction == "up" else int(np.min(np.where(cap)[0]))
        for j in range(mask.shape[1]):
            for k in range(mask.shape[2]):
                hit = np.flatnonzero(mask[:, j, k])
                cap_hit = np.flatnonzero(cap[:, j, k])
                if hit.size == 0 or cap_hit.size == 0:
                    continue
                contact = int(hit.max() if direction == "up" else hit.min())
                start, stop = (contact + 1, outer + 1) if direction == "up" else (outer, contact)
                if start < stop:
                    flat_cap[start:stop, j, k] = True
                    thickness_records.append(((j, k), int(stop - start)))
    elif extrusion_axis == 1:
        outer = int(np.max(np.where(cap)[1])) if direction == "up" else int(np.min(np.where(cap)[1]))
        for i in range(mask.shape[0]):
            for k in range(mask.shape[2]):
                hit = np.flatnonzero(mask[i, :, k])
                cap_hit = np.flatnonzero(cap[i, :, k])
                if hit.size == 0 or cap_hit.size == 0:
                    continue
                contact = int(hit.max() if direction == "up" else hit.min())
                start, stop = (contact + 1, outer + 1) if direction == "up" else (outer, contact)
                if start < stop:
                    flat_cap[i, start:stop, k] = True
                    thickness_records.append(((i, k), int(stop - start)))
    else:
        outer = int(np.max(np.where(cap)[2])) if direction == "up" else int(np.min(np.where(cap)[2]))
        for i in range(mask.shape[0]):
            for j in range(mask.shape[1]):
                hit = np.flatnonzero(mask[i, j, :])
                cap_hit = np.flatnonzero(cap[i, j, :])
                if hit.size == 0 or cap_hit.size == 0:
                    continue
                contact = int(hit.max() if direction == "up" else hit.min())
                start, stop = (contact + 1, outer + 1) if direction == "up" else (outer, contact)
                if start < stop:
                    flat_cap[i, j, start:stop] = True
                    thickness_records.append(((i, j), int(stop - start)))

    if thickness_records:
        values = np.asarray([thick for _, thick in thickness_records], dtype=float)
        upper = max(float(thickness * 2), 2.0 * float(np.mean(values)))
        for coords, thick in thickness_records:
            if float(thick) <= upper:
                continue
            if extrusion_axis == 0:
                j, k = coords
                flat_cap[:, j, k] = False
            elif extrusion_axis == 1:
                i, k = coords
                flat_cap[i, :, k] = False
            else:
                i, j = coords
                flat_cap[i, j, :] = False
    return flat_cap


def _projected_footprint(surface, extrusion_axis, shape):
    import numpy as np
    from scipy.ndimage import binary_fill_holes

    footprint = surface.any(axis=extrusion_axis)
    if not np.any(footprint):
        return footprint

    coords = np.array(np.where(footprint))
    mins = coords.min(axis=1)
    maxs = coords.max(axis=1)
    if shape == "box":
        out = np.zeros_like(footprint, dtype=bool)
        out[tuple(slice(mins[i], maxs[i] + 1) for i in range(2))] = True
        return out
    if shape == "round":
        grids = np.ogrid[tuple(slice(0, size) for size in footprint.shape)]
        center = (mins + maxs) / 2.0
        radii = np.maximum((maxs - mins + 1) / 2.0, 0.5)
        distance = sum(((grids[i] - center[i]) / radii[i]) ** 2 for i in range(2))
        return distance <= 1.0
    if shape == "fit":
        return binary_fill_holes(footprint)
    raise ValueError("shape must be 'fit', 'box', or 'round'.")


def _selector_for_material_column(lateral, extrusion_axis):
    """Build an index selector for one column running along the extrusion axis."""
    selector = [slice(None), slice(None), slice(None)]
    lateral_iter = iter(lateral)
    for axis in range(3):
        if axis != extrusion_axis:
            selector[axis] = next(lateral_iter)
    return selector


def _largest_planar_component(mask):
    """Keep the largest connected region in a 2D projected footprint."""
    import numpy as np

    values = np.asarray(mask, dtype=bool)
    visited = np.zeros(values.shape, dtype=bool)
    best = []
    for start_array in np.argwhere(values):
        start = tuple(int(v) for v in start_array)
        if visited[start]:
            continue
        component = []
        queue = [start]
        visited[start] = True
        while queue:
            coord = queue.pop(0)
            component.append(coord)
            for axis in range(2):
                for offset in (-1, 1):
                    neighbor = [coord[0], coord[1]]
                    neighbor[axis] += offset
                    if neighbor[axis] < 0 or neighbor[axis] >= values.shape[axis]:
                        continue
                    token = tuple(neighbor)
                    if visited[token] or not values[token]:
                        continue
                    visited[token] = True
                    queue.append(token)
        if len(component) > len(best):
            best = component
    out = np.zeros(values.shape, dtype=bool)
    if best:
        coords = tuple(np.asarray(best, dtype=np.int64).T)
        out[coords] = True
    return out


def _fill_short_boolean_gaps_in_line(line, max_gap):
    """Fill short False runs when they are bracketed by True values."""
    import numpy as np

    out = np.asarray(line, dtype=bool).copy()
    false_runs = np.flatnonzero(~out)
    if false_runs.size == 0:
        return out
    start = 0
    while start < false_runs.size:
        stop = start + 1
        while stop < false_runs.size and false_runs[stop] == false_runs[stop - 1] + 1:
            stop += 1
        run = false_runs[start:stop]
        if (
            run.size <= int(max_gap)
            and int(run[0]) > 0
            and int(run[-1]) < out.size - 1
            and out[int(run[0]) - 1]
            and out[int(run[-1]) + 1]
        ):
            out[run] = True
        start = stop
    return out


def _fill_short_boolean_gaps_in_footprint(values, max_gap=2):
    """Fill small row-wise and column-wise holes in a projected footprint."""
    import numpy as np

    out = np.asarray(values, dtype=bool).copy()
    for row in range(out.shape[0]):
        out[row, :] = _fill_short_boolean_gaps_in_line(out[row, :], max_gap)
    for col in range(out.shape[1]):
        out[:, col] = _fill_short_boolean_gaps_in_line(out[:, col], max_gap)
    return out


def _fill_short_boolean_gaps_in_volume(values, max_gap=2):
    """Fill short bracketed gaps along each axis in a 3D mask."""
    import numpy as np

    out = np.asarray(values, dtype=bool).copy()
    for y in range(out.shape[1]):
        for z in range(out.shape[2]):
            out[:, y, z] = _fill_short_boolean_gaps_in_line(out[:, y, z], max_gap)
    for x in range(out.shape[0]):
        for z in range(out.shape[2]):
            out[x, :, z] = _fill_short_boolean_gaps_in_line(out[x, :, z], max_gap)
    for x in range(out.shape[0]):
        for y in range(out.shape[1]):
            out[x, y, :] = _fill_short_boolean_gaps_in_line(out[x, y, :], max_gap)
    return out


def _clean_projected_footprint(mask):
    """Clean a candidate cap footprint before it is extruded into a disk."""
    import numpy as np
    from scipy.ndimage import binary_fill_holes

    values = np.asarray(mask, dtype=bool)
    values = _fill_short_boolean_gaps_in_footprint(values, max_gap=2)
    values = binary_fill_holes(values).astype(bool)
    return _largest_planar_component(values)


def _apply_requested_footprint_shape(mask, *, shape):
    """Convert a projected support mask to the requested disk footprint."""
    import numpy as np

    values = np.asarray(mask, dtype=bool)
    if not np.any(values):
        return values
    mode = str(shape).strip().lower()
    if mode in {"", "anatomy", "fit", "projected", "mask"}:
        return values
    coords = np.argwhere(values)
    lo = coords.min(axis=0)
    hi = coords.max(axis=0)
    yy, xx = np.indices(values.shape, dtype=np.float64)
    center = (lo + hi) / 2.0
    half = np.maximum((hi - lo + 1) / 2.0, 0.5)
    if mode in {"box", "rectangle", "rectangular", "square"}:
        shaped = np.ones(values.shape, dtype=bool)
    elif mode in {"round", "circle", "circular"}:
        norm_y = (yy - center[0]) / half[0]
        norm_x = (xx - center[1]) / half[1]
        shaped = (norm_y * norm_y + norm_x * norm_x) <= 1.0
    elif mode in {"hex", "hexagon", "hexagonal"}:
        norm_y = np.abs((yy - center[0]) / half[0])
        norm_x = np.abs((xx - center[1]) / half[1])
        shaped = (norm_x <= 1.0) & (norm_y <= 1.0) & (norm_x + norm_y / 2.0 <= 1.0)
    else:
        raise ValueError("disk shape must be one of anatomy, rectangle, square, round, or hex")
    return shaped


def _unit_vector(vector, name):
    """Return a unit vector from a physical-space direction."""
    import numpy as np

    values = np.asarray(vector, dtype=float)
    if values.shape != (3,):
        raise ValueError(f"{name} must contain three values.")
    norm = float(np.linalg.norm(values))
    if norm <= 1.0e-12:
        raise ValueError(f"{name} must have nonzero length.")
    return values / norm


def _plane_axes_from_normal(normal):
    """Return a stable orthonormal in-plane basis for a plane normal."""
    import numpy as np

    helper = np.asarray([1.0, 0.0, 0.0], dtype=float)
    if abs(float(np.dot(normal, helper))) > 0.85:
        helper = np.asarray([0.0, 1.0, 0.0], dtype=float)
    u_axis = _unit_vector(np.cross(normal, helper), "u_axis")
    v_axis = _unit_vector(np.cross(normal, u_axis), "v_axis")
    return u_axis, v_axis


def pad_vtk_images_to_physical_bounds(vtk_images, *, desired_bounds, constants=None):
    """Pad aligned VTK images so the requested physical x/y/z bounds fit."""
    import numpy as np
    import vtk

    images = list(vtk_images)
    if not images:
        return [], {"lower": (0, 0, 0), "upper": (0, 0, 0)}
    constants = [0] * len(images) if constants is None else list(constants)
    if len(constants) != len(images):
        raise ValueError("constants must match vtk_images length.")
    if len(desired_bounds) != 6:
        raise ValueError("desired_bounds must contain x/y/z min/max values.")

    reference = images[0]
    spacing = np.asarray(reference.GetSpacing(), dtype=float)
    if np.any(spacing <= 0.0):
        raise ValueError("VTK image spacing must be positive.")
    extent = reference.GetExtent()
    current_min = np.asarray(
        [
            reference.GetOrigin()[0] + float(extent[0]) * spacing[0],
            reference.GetOrigin()[1] + float(extent[2]) * spacing[1],
            reference.GetOrigin()[2] + float(extent[4]) * spacing[2],
        ],
        dtype=float,
    )
    current_max = np.asarray(
        [
            reference.GetOrigin()[0] + float(extent[1]) * spacing[0],
            reference.GetOrigin()[1] + float(extent[3]) * spacing[1],
            reference.GetOrigin()[2] + float(extent[5]) * spacing[2],
        ],
        dtype=float,
    )
    desired = np.asarray(desired_bounds, dtype=float)
    desired_min = desired[[0, 2, 4]]
    desired_max = desired[[1, 3, 5]]
    lower = np.maximum(0, np.ceil((current_min - desired_min) / spacing).astype(int))
    upper = np.maximum(0, np.ceil((desired_max - current_max) / spacing).astype(int))

    if not np.any(lower) and not np.any(upper):
        return images, {
            "lower": tuple(int(value) for value in lower),
            "upper": tuple(int(value) for value in upper),
        }

    output_extent = (
        int(extent[0] - lower[0]),
        int(extent[1] + upper[0]),
        int(extent[2] - lower[1]),
        int(extent[3] + upper[1]),
        int(extent[4] - lower[2]),
        int(extent[5] + upper[2]),
    )
    if len(images) != len(constants):
        raise ValueError("images and constants must contain the same number of items.")

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
        image_origin = image.GetOrigin()
        image_spacing = out.GetSpacing()
        out.SetOrigin(
            float(image_origin[0]) + float(output_extent[0]) * float(image_spacing[0]),
            float(image_origin[1]) + float(output_extent[2]) * float(image_spacing[1]),
            float(image_origin[2]) + float(output_extent[4]) * float(image_spacing[2]),
        )
        out.SetExtent(0, dims_out[0] - 1, 0, dims_out[1] - 1, 0, dims_out[2] - 1)
        padded.append(out)

    return padded, {
        "lower": tuple(int(value) for value in lower),
        "upper": tuple(int(value) for value in upper),
    }


def fit_vtk_images_to_physical_bounds(vtk_images, *, desired_bounds, constants=None):
    """Crop or pad aligned VTK images to the requested physical bounds.

    The output extent is chosen from voxel-center coordinates, then reset to a
    zero-based extent with the origin shifted so physical coordinates are
    preserved. Values outside the source extent are filled from ``constants``.
    """
    import numpy as np
    import vtk

    images = list(vtk_images)
    if not images:
        return [], {"output_extent": (0, -1, 0, -1, 0, -1)}
    constants = [0] * len(images) if constants is None else list(constants)
    if len(constants) != len(images):
        raise ValueError("constants must match vtk_images length.")
    if len(desired_bounds) != 6:
        raise ValueError("desired_bounds must contain x/y/z min/max values.")

    reference = images[0]
    spacing = np.asarray(reference.GetSpacing(), dtype=float)
    if np.any(spacing <= 0.0):
        raise ValueError("VTK image spacing must be positive.")
    origin = np.asarray(reference.GetOrigin(), dtype=float)
    desired = np.asarray(desired_bounds, dtype=float)
    desired_min = desired[[0, 2, 4]]
    desired_max = desired[[1, 3, 5]]
    lower_index = np.floor((desired_min - origin) / spacing).astype(int)
    upper_index = np.ceil((desired_max - origin) / spacing).astype(int)
    output_extent = (
        int(lower_index[0]),
        int(upper_index[0]),
        int(lower_index[1]),
        int(upper_index[1]),
        int(lower_index[2]),
        int(upper_index[2]),
    )

    fitted = []
    for image, constant in zip(images, constants):
        pad = vtk.vtkImageConstantPad()
        pad.SetInputData(image)
        pad.SetOutputWholeExtent(output_extent)
        pad.SetConstant(float(constant))
        pad.Update()
        out = vtk.vtkImageData()
        out.DeepCopy(pad.GetOutput())
        dims_out = out.GetDimensions()
        image_origin = image.GetOrigin()
        image_spacing = out.GetSpacing()
        out.SetOrigin(
            float(image_origin[0]) + float(output_extent[0]) * float(image_spacing[0]),
            float(image_origin[1]) + float(output_extent[2]) * float(image_spacing[1]),
            float(image_origin[2]) + float(output_extent[4]) * float(image_spacing[2]),
        )
        out.SetExtent(0, dims_out[0] - 1, 0, dims_out[1] - 1, 0, dims_out[2] - 1)
        fitted.append(out)

    return fitted, {"output_extent": output_extent}


def projected_material_disk_required_bounds(
    surface_vtk_image,
    *,
    center,
    normal,
    u_axis=None,
    v_axis=None,
    size=(24.0, 24.0),
    shape="anatomy",
    thickness=3.0,
    intrusion=2.0,
    stable_surface=False,
    stable_surface_fraction=0.55,
    stable_surface_min_area_fraction=0.12,
    stable_surface_max_depth=None,
):
    """Return conservative physical bounds needed for a plane-projected disk."""
    import numpy as np

    active = vtk_image_to_numpy(surface_vtk_image) != 0
    if not np.any(active):
        return None

    spacing = tuple(float(value) for value in surface_vtk_image.GetSpacing())
    origin = tuple(float(value) for value in surface_vtk_image.GetOrigin())
    extent = surface_vtk_image.GetExtent()
    extent_start = (extent[0], extent[2], extent[4])
    center = np.asarray(center, dtype=float)
    if center.shape != (3,):
        raise ValueError("center must contain three values.")
    normal = _unit_vector(normal, "normal")
    if u_axis is None or v_axis is None:
        u_axis, v_axis = _plane_axes_from_normal(normal)
    else:
        u_axis = _unit_vector(u_axis, "u_axis")
        v_axis = _unit_vector(v_axis, "v_axis")
    if len(size) != 2:
        raise ValueError("size must contain two in-plane lengths.")

    half_u = max(float(size[0]) / 2.0, 0.0)
    half_v = max(float(size[1]) / 2.0, 0.0)
    active_indices = np.argwhere(active)
    active_points = _physical_points_from_indices(
        active_indices,
        spacing=spacing,
        origin=origin,
        extent_start=extent_start,
    )
    rel = active_points - center
    distance = rel @ normal
    u = rel @ u_axis
    v = rel @ v_axis
    tolerance = max(min(spacing) * 0.75, 1.0e-6)
    inside = _inside_projected_shape(
        shape,
        u,
        v,
        half_u,
        half_v,
        tolerance=tolerance,
    )
    forward = distance >= -tolerance
    if not np.any(inside & forward):
        return None

    first_distance, distance_by_key = _surface_distance_by_projected_bucket(
        active_indices[inside & forward],
        distance[inside & forward],
        u[inside & forward],
        v[inside & forward],
        spacing=spacing,
    )
    surface_distance = first_distance
    if stable_surface and distance_by_key:
        stable_contact = _stable_surface_contact_from_bucket_distances(
            distance_by_key,
            stable_surface_fraction=stable_surface_fraction,
            stable_surface_min_area_fraction=stable_surface_min_area_fraction,
            stable_surface_max_depth=stable_surface_max_depth,
        )
        surface_distance = float(stable_contact["stable_distance"])
    cap_inner_distance = surface_distance + max(float(intrusion), 0.0)
    cap_outer_distance = cap_inner_distance - max(float(thickness), 0.0)
    d_min = min(cap_inner_distance, cap_outer_distance)
    d_max = max(cap_inner_distance, cap_outer_distance)
    corners = []
    for u_value in (-half_u, half_u):
        for v_value in (-half_v, half_v):
            for d_value in (d_min, d_max):
                corners.append(center + u_value * u_axis + v_value * v_axis + d_value * normal)
    points = np.vstack(corners)
    margin = np.asarray(spacing, dtype=float)
    lower = points.min(axis=0) - margin
    upper = points.max(axis=0) + margin
    return (
        float(lower[0]),
        float(upper[0]),
        float(lower[1]),
        float(upper[1]),
        float(lower[2]),
        float(upper[2]),
    )


def _physical_points_from_indices(indices, *, spacing, origin, extent_start=(0, 0, 0)):
    """Convert x/y/z array indices to physical image coordinates."""
    import numpy as np

    idx = np.asarray(indices, dtype=float)
    return np.asarray(origin, dtype=float) + (
        idx + np.asarray(extent_start, dtype=float)
    ) * np.asarray(spacing, dtype=float)


def _projected_bucket_key(u, v, *, spacing):
    resolution = max(min(float(value) for value in spacing), 1.0e-6)
    return int(round(float(u) / resolution)), int(round(float(v) / resolution))


def _inside_projected_shape(shape, u, v, half_u, half_v, *, tolerance):
    import numpy as np

    token = str(shape).strip().lower()
    if token in {"anatomy", "rectangle", "rectangular", "box"}:
        return (np.abs(u) <= half_u + tolerance) & (np.abs(v) <= half_v + tolerance)
    if token == "square":
        half = min(float(half_u), float(half_v))
        return (np.abs(u) <= half + tolerance) & (np.abs(v) <= half + tolerance)
    if token in {"round", "circle", "circular", "oval"}:
        half_u = max(float(half_u) + tolerance, 1.0e-9)
        half_v = max(float(half_v) + tolerance, 1.0e-9)
        return ((u / half_u) ** 2 + (v / half_v) ** 2) <= 1.0
    if token in {"hex", "hexagon", "hexagonal"}:
        half = max(min(float(half_u), float(half_v)), 1.0e-9)
        uu = u / half
        vv = v / half
        hex_tol = tolerance / half
        return (
            (np.abs(uu) <= 1.0 + hex_tol)
            & (np.abs(0.5 * uu + 0.8660254 * vv) <= 1.0 + hex_tol)
            & (np.abs(0.5 * uu - 0.8660254 * vv) <= 1.0 + hex_tol)
        )
    raise ValueError("disk shape must be one of anatomy, rectangle, square, round, or hex.")


def _surface_distance_by_projected_bucket(
    indices,
    distance,
    u,
    v,
    *,
    spacing,
):
    import numpy as np

    best = {}
    lengths = {len(indices), len(distance), len(u), len(v)}
    if len(lengths) != 1:
        raise ValueError("indices, distance, u, and v must have the same length.")

    for voxel, dist, uu, vv in zip(indices, distance, u, v):
        if float(dist) < 0.0:
            continue
        key = _projected_bucket_key(float(uu), float(vv), spacing=spacing)
        current = best.get(key)
        if current is None or float(dist) < current[0]:
            best[key] = (float(dist), np.asarray(voxel, dtype=np.int64))
    if not best:
        return 0.0, {}
    first_distance = min(value[0] for value in best.values())
    return float(first_distance), {key: float(value[0]) for key, value in best.items()}


def _clean_footprint_keys(keys):
    """Return the largest filled component from projected footprint keys."""
    import numpy as np

    if not keys:
        return set()
    coords = np.asarray(list(keys), dtype=int)
    lower = coords.min(axis=0)
    upper = coords.max(axis=0)
    footprint = np.zeros(tuple((upper - lower + 1).tolist()), dtype=bool)
    shifted = coords - lower
    footprint[shifted[:, 0], shifted[:, 1]] = True
    cleaned = _clean_projected_footprint(footprint)
    cleaned_coords = np.argwhere(cleaned) + lower
    return {tuple(int(value) for value in coord) for coord in cleaned_coords}


def _stable_surface_contact_from_bucket_distances(
    distance_by_key,
    *,
    stable_surface_fraction=0.55,
    stable_surface_min_area_fraction=0.12,
    stable_surface_max_depth=None,
):
    """Find the first broad projected surface behind isolated protrusions."""
    import numpy as np

    if not distance_by_key:
        return {
            "first_distance": 0.0,
            "stable_distance": 0.0,
            "stable_depth": 0.0,
            "first_area": 0,
            "stable_area": 0,
            "peak_area": 0,
            "stable_keys": set(),
        }

    distances = np.asarray(list(distance_by_key.values()), dtype=float)
    first_distance = float(np.min(distances))
    full_area = int(len(distance_by_key))
    first_keys = {
        key
        for key, distance in distance_by_key.items()
        if float(distance) <= first_distance + 1.0e-9
    }
    first_area = len(_clean_footprint_keys(first_keys))

    if stable_surface_max_depth is None:
        max_distance = float(np.max(distances))
    else:
        max_distance = first_distance + max(float(stable_surface_max_depth), 0.0)
    thresholds = np.unique(distances[distances <= max_distance + 1.0e-9])
    if thresholds.size == 0:
        thresholds = np.asarray([first_distance], dtype=float)

    surfaces = []
    peak_area = 0
    for threshold in thresholds:
        keys = {
            key
            for key, distance in distance_by_key.items()
            if float(distance) <= float(threshold) + 1.0e-9
        }
        clean_keys = _clean_footprint_keys(keys)
        area = len(clean_keys)
        peak_area = max(peak_area, area)
        surfaces.append((float(threshold), area, clean_keys))

    target_area = max(
        1,
        int(np.ceil(max(float(stable_surface_fraction), 0.0) * peak_area)),
        int(np.ceil(max(float(stable_surface_min_area_fraction), 0.0) * full_area)),
    )
    stable_distance, stable_area, stable_keys = surfaces[-1]
    for threshold, area, keys in surfaces:
        if area >= target_area:
            stable_distance = threshold
            stable_area = area
            stable_keys = keys
            break

    return {
        "first_distance": first_distance,
        "stable_distance": float(stable_distance),
        "stable_depth": float(stable_distance - first_distance),
        "first_area": int(first_area),
        "stable_area": int(stable_area),
        "peak_area": int(peak_area),
        "stable_keys": stable_keys,
    }


def projected_stable_surface_contact(
    active_mask,
    *,
    spacing,
    origin,
    center,
    normal,
    u_axis=None,
    v_axis=None,
    size=(24.0, 24.0),
    shape="square",
    extent_start=(0, 0, 0),
    stable_surface_fraction=0.55,
    stable_surface_min_area_fraction=0.12,
    stable_surface_max_depth=None,
):
    """Return stable projected contact-surface metrics for a physical plane."""
    import numpy as np

    active = np.asarray(active_mask, dtype=bool)
    if active.ndim != 3:
        raise ValueError("active_mask must be a 3D x/y/z array.")
    if not np.any(active):
        return _stable_surface_contact_from_bucket_distances({})

    spacing = tuple(float(value) for value in spacing)
    origin = tuple(float(value) for value in origin)
    center = np.asarray(center, dtype=float)
    if center.shape != (3,):
        raise ValueError("center must contain three values.")
    normal = _unit_vector(normal, "normal")
    if u_axis is None or v_axis is None:
        u_axis, v_axis = _plane_axes_from_normal(normal)
    else:
        u_axis = _unit_vector(u_axis, "u_axis")
        v_axis = _unit_vector(v_axis, "v_axis")
    if len(size) != 2:
        raise ValueError("size must contain two in-plane lengths.")
    half_u = max(float(size[0]) / 2.0, 0.0)
    half_v = max(float(size[1]) / 2.0, 0.0)

    indices = np.argwhere(active)
    points = _physical_points_from_indices(
        indices,
        spacing=spacing,
        origin=origin,
        extent_start=extent_start,
    )
    rel = points - center
    distance = rel @ normal
    u = rel @ u_axis
    v = rel @ v_axis
    tolerance = max(min(spacing) * 0.75, 1.0e-6)
    inside = _inside_projected_shape(
        shape,
        u,
        v,
        half_u,
        half_v,
        tolerance=tolerance,
    )
    forward = distance >= -tolerance
    candidate = inside & forward
    if not np.any(candidate):
        return _stable_surface_contact_from_bucket_distances({})
    _, distance_by_key = _surface_distance_by_projected_bucket(
        indices[candidate],
        distance[candidate],
        u[candidate],
        v[candidate],
        spacing=spacing,
    )
    return _stable_surface_contact_from_bucket_distances(
        distance_by_key,
        stable_surface_fraction=stable_surface_fraction,
        stable_surface_min_area_fraction=stable_surface_min_area_fraction,
        stable_surface_max_depth=stable_surface_max_depth,
    )


def generate_projected_material_disk_mask(
    active_mask,
    *,
    spacing,
    origin,
    center,
    normal,
    u_axis=None,
    v_axis=None,
    size=(24.0, 24.0),
    shape="square",
    thickness=3.0,
    intrusion=2.0,
    anatomy_constrained=True,
    stable_surface=False,
    stable_surface_fraction=0.55,
    stable_surface_min_area_fraction=0.12,
    stable_surface_max_depth=None,
    close_gaps_3d=False,
    close_gaps_3d_max_gap=2,
    material_mask=None,
    extent_start=(0, 0, 0),
):
    """Generate a physical-space material disk from an interactive contact plane.

    The plane normal points from the authored plane toward the anatomy. The
    disk has a fixed total ``thickness`` in physical units. ``intrusion`` lets
    the local anatomy occupy that much of the fixed disk thickness; generated
    disk voxels are then cleared anywhere material already exists, so bone is
    preserved.
    """
    import numpy as np

    active = np.asarray(active_mask, dtype=bool)
    if active.ndim != 3:
        raise ValueError("active_mask must be a 3D x/y/z array.")
    material = active if material_mask is None else np.asarray(material_mask) > 0
    if material.shape != active.shape:
        raise ValueError("material_mask must match active_mask shape.")
    if not np.any(active):
        return np.zeros(active.shape, dtype=bool)

    spacing = tuple(float(value) for value in spacing)
    origin = tuple(float(value) for value in origin)
    center = np.asarray(center, dtype=float)
    if center.shape != (3,):
        raise ValueError("center must contain three values.")
    normal = _unit_vector(normal, "normal")
    if u_axis is None or v_axis is None:
        u_axis, v_axis = _plane_axes_from_normal(normal)
    else:
        u_axis = _unit_vector(u_axis, "u_axis")
        v_axis = _unit_vector(v_axis, "v_axis")
    if len(size) != 2:
        raise ValueError("size must contain two in-plane lengths.")
    half_u = max(float(size[0]) / 2.0, 0.5)
    half_v = max(float(size[1]) / 2.0, 0.5)
    thickness = max(float(thickness), 0.0)
    intrusion = max(float(intrusion), 0.0)
    if thickness <= 0.0:
        return np.zeros(active.shape, dtype=bool)

    active_indices = np.argwhere(active)
    active_points = _physical_points_from_indices(
        active_indices,
        spacing=spacing,
        origin=origin,
        extent_start=extent_start,
    )
    active_rel = active_points - center
    active_distance = active_rel @ normal
    active_u = active_rel @ u_axis
    active_v = active_rel @ v_axis
    tolerance = max(min(spacing) * 0.75, 1.0e-6)
    active_inside = _inside_projected_shape(
        shape,
        active_u,
        active_v,
        half_u,
        half_v,
        tolerance=tolerance,
    )
    forward = active_distance >= -tolerance
    candidate = active_inside & forward
    if not np.any(candidate):
        return np.zeros(active.shape, dtype=bool)

    first_distance, distance_by_key = _surface_distance_by_projected_bucket(
        active_indices[candidate],
        active_distance[candidate],
        active_u[candidate],
        active_v[candidate],
        spacing=spacing,
    )
    if not distance_by_key:
        return np.zeros(active.shape, dtype=bool)

    support_distance = first_distance
    stable_keys = None
    if stable_surface:
        stable_contact = _stable_surface_contact_from_bucket_distances(
            distance_by_key,
            stable_surface_fraction=stable_surface_fraction,
            stable_surface_min_area_fraction=stable_surface_min_area_fraction,
            stable_surface_max_depth=stable_surface_max_depth,
        )
        support_distance = float(stable_contact["stable_distance"])
        footprint_min = support_distance - max(2.0 * min(spacing), tolerance)
        footprint_max = support_distance + intrusion + tolerance
        stable_keys = _clean_footprint_keys(
            {
                key
                for key, distance in distance_by_key.items()
                if footprint_min <= float(distance) <= footprint_max
            }
        )
        if not stable_keys:
            return np.zeros(active.shape, dtype=bool)

    if anatomy_constrained:
        depth_limit = support_distance + intrusion + tolerance
        distance_by_key = {
            key: distance
            for key, distance in distance_by_key.items()
            if distance <= depth_limit and (stable_keys is None or key in stable_keys)
        }
        if not distance_by_key:
            return np.zeros(active.shape, dtype=bool)

    cap_inner_distance = support_distance + intrusion
    cap_outer_distance = cap_inner_distance - thickness
    full_indices = np.argwhere(np.ones(active.shape, dtype=bool))
    full_points = _physical_points_from_indices(
        full_indices,
        spacing=spacing,
        origin=origin,
        extent_start=extent_start,
    )
    rel = full_points - center
    distance = rel @ normal
    u = rel @ u_axis
    v = rel @ v_axis
    inside = _inside_projected_shape(
        shape,
        u,
        v,
        half_u,
        half_v,
        tolerance=tolerance,
    )
    empty = ~material[tuple(full_indices.T)]

    if anatomy_constrained:
        if len(u) != len(v):
            raise ValueError("u and v must have the same length.")
        keys = [_projected_bucket_key(float(uu), float(vv), spacing=spacing) for uu, vv in zip(u, v)]
        local_surface = np.asarray(
            [distance_by_key.get(key, np.nan) for key in keys],
            dtype=float,
        )
        bucket_mask = np.isfinite(local_surface)
        flat_outer = cap_outer_distance
        if flat_outer >= support_distance:
            flat_outer = support_distance - thickness
        if stable_surface:
            depth_ok = (
                (distance >= flat_outer - tolerance)
                & (distance <= cap_inner_distance + tolerance)
            )
        else:
            local_min = np.minimum(flat_outer, local_surface)
            local_max = np.maximum(flat_outer, local_surface)
            depth_ok = (distance >= local_min - tolerance) & (distance <= local_max + tolerance)
    else:
        bucket_mask = np.ones(full_indices.shape[0], dtype=bool)
        depth_ok = (
            (distance >= cap_outer_distance - tolerance)
            & (distance <= cap_inner_distance + tolerance)
        )

    keep = inside & depth_ok & empty & bucket_mask
    disk = np.zeros(active.shape, dtype=bool)
    if np.any(keep):
        disk[tuple(full_indices[keep].T)] = True
    if close_gaps_3d and np.any(disk):
        disk = _fill_short_boolean_gaps_in_volume(
            disk,
            max_gap=int(close_gaps_3d_max_gap),
        )
        disk &= ~material
    return disk


def _generate_projected_bone_caps_from_mask(mask, *, extrusion_axis, thickness, intrusion_depth, shape):
    """Generate inferior and superior bone caps from projected surface footprints.

    The builder looks at each material column, records the first and last bone
    voxel, then builds a flat cap from columns close enough to the global
    inferior/superior surfaces.
    """
    import numpy as np

    thickness = max(1, int(thickness))
    intrusion = max(
        0,
        int(round(thickness * 2.5)) if intrusion_depth is None else int(intrusion_depth),
    )
    inferior = np.zeros(mask.shape, dtype=bool)
    superior = np.zeros(mask.shape, dtype=bool)
    lateral_shape = tuple(mask.shape[idx] for idx in range(3) if idx != extrusion_axis)
    inferior_surface_by_column = np.full(lateral_shape, -1, dtype=np.int32)
    superior_surface_by_column = np.full(lateral_shape, -1, dtype=np.int32)
    for lateral in np.ndindex(lateral_shape):
        selector = _selector_for_material_column(lateral, extrusion_axis)
        occupied = np.flatnonzero(mask[tuple(selector)])
        if occupied.size == 0:
            continue
        inferior_surface_by_column[lateral] = int(occupied.min())
        superior_surface_by_column[lateral] = int(occupied.max())
    has_bone_in_column = inferior_surface_by_column >= 0
    if not np.any(has_bone_in_column):
        return inferior, superior

    global_inferior_surface = int(np.min(inferior_surface_by_column[has_bone_in_column]))
    global_superior_surface = int(np.max(superior_surface_by_column[has_bone_in_column]))
    inferior_depth_limit = min(
        global_inferior_surface + intrusion,
        int(np.max(inferior_surface_by_column[has_bone_in_column])),
    )
    superior_depth_limit = max(
        global_superior_surface - intrusion,
        int(np.min(superior_surface_by_column[has_bone_in_column])),
    )
    inferior_footprint = _apply_requested_footprint_shape(
        _clean_projected_footprint(
            has_bone_in_column & (inferior_surface_by_column <= inferior_depth_limit)
        ),
        shape=shape,
    )
    superior_footprint = _apply_requested_footprint_shape(
        _clean_projected_footprint(
            has_bone_in_column & (superior_surface_by_column >= superior_depth_limit)
        ),
        shape=shape,
    )

    inferior_start = max(0, global_inferior_surface + intrusion - thickness)
    inferior_stop = min(mask.shape[extrusion_axis], global_inferior_surface + intrusion)
    superior_start = max(0, global_superior_surface - intrusion + 1)
    superior_stop = min(mask.shape[extrusion_axis], superior_start + thickness)

    for lateral in np.ndindex(lateral_shape):
        if inferior_footprint[lateral]:
            selector = _selector_for_material_column(lateral, extrusion_axis)
            selector[extrusion_axis] = slice(inferior_start, inferior_stop)
            inferior[tuple(selector)] = True
        if superior_footprint[lateral]:
            selector = _selector_for_material_column(lateral, extrusion_axis)
            selector[extrusion_axis] = slice(superior_start, superior_stop)
            superior[tuple(selector)] = True
    return inferior, superior


def _generate_axis_aligned_anatomy_cap_from_mask(mask, *, extrusion_axis, direction, thickness, intrusion_depth):
    """Generate an axis-aligned cap from the local anatomy surface.

    The cap has a fixed total thickness. ``intrusion_depth`` does not mean
    "overwrite this much bone". It means that the anatomy is allowed to occupy
    that much of the fixed cap thickness before the PMMA part starts. After the
    cap is built, any overlap with bone is cleared so bone material is preserved.
    """
    import numpy as np

    thickness = max(1, int(thickness))
    intrusion = max(0, int(intrusion_depth))
    cap = np.zeros(mask.shape, dtype=bool)
    lateral_shape = tuple(mask.shape[idx] for idx in range(3) if idx != extrusion_axis)
    surface_position_by_column = np.full(lateral_shape, -1, dtype=np.int32)

    for lateral in np.ndindex(lateral_shape):
        selector = _selector_for_material_column(lateral, extrusion_axis)
        occupied = np.flatnonzero(mask[tuple(selector)])
        if occupied.size == 0:
            continue
        surface_position_by_column[lateral] = int(occupied.max() if direction == "up" else occupied.min())

    has_bone_in_column = surface_position_by_column >= 0
    if not np.any(has_bone_in_column):
        return cap

    # Keep columns whose surface is close enough to the outermost contact
    # surface. Intrusion is the amount of the fixed cap thickness that anatomy
    # may occupy; deeper columns would make the PMMA thicker than requested.
    maximum_supported_surface_depth = intrusion
    if direction == "up":
        flat_contact_position = int(np.max(surface_position_by_column[has_bone_in_column]))
        included_columns = has_bone_in_column & (
            surface_position_by_column >= flat_contact_position - maximum_supported_surface_depth
        )
        outer_cap_position = flat_contact_position - intrusion + thickness
        for lateral in np.argwhere(included_columns):
            lateral = tuple(int(value) for value in lateral)
            bone_surface_position = int(surface_position_by_column[lateral])
            selector = _selector_for_material_column(lateral, extrusion_axis)
            selector[extrusion_axis] = slice(
                max(bone_surface_position + 1, 0),
                min(outer_cap_position + 1, mask.shape[extrusion_axis]),
            )
            cap[tuple(selector)] = True
    elif direction == "down":
        flat_contact_position = int(np.min(surface_position_by_column[has_bone_in_column]))
        included_columns = has_bone_in_column & (
            surface_position_by_column <= flat_contact_position + maximum_supported_surface_depth
        )
        outer_cap_position = flat_contact_position + intrusion - thickness
        for lateral in np.argwhere(included_columns):
            lateral = tuple(int(value) for value in lateral)
            bone_surface_position = int(surface_position_by_column[lateral])
            selector = _selector_for_material_column(lateral, extrusion_axis)
            selector[extrusion_axis] = slice(
                max(outer_cap_position, 0),
                min(bone_surface_position, mask.shape[extrusion_axis]),
            )
            cap[tuple(selector)] = True
    else:
        raise ValueError("direction must be 'up' or 'down'.")

    cap[mask] = False
    return cap


def _shaped_contact_cap(mask, surface, extrusion_axis, direction, thickness, shape):
    import numpy as np

    footprint = _projected_footprint(surface, extrusion_axis, shape)
    cap = np.zeros_like(mask, dtype=bool)
    hit = np.where(mask)[extrusion_axis]
    outer = int(hit.max() + thickness if direction == "up" else hit.min() - thickness)
    outer = int(np.clip(outer, 0, mask.shape[extrusion_axis] - 1))

    footprint_coords = np.array(np.where(footprint)).T
    for coord in footprint_coords:
        if extrusion_axis == 0:
            j, k = coord
            ray_hit = np.flatnonzero(mask[:, j, k])
        elif extrusion_axis == 1:
            i, k = coord
            ray_hit = np.flatnonzero(mask[i, :, k])
        else:
            i, j = coord
            ray_hit = np.flatnonzero(mask[i, j, :])
        if ray_hit.size == 0:
            continue
        contact = int(ray_hit.max() if direction == "up" else ray_hit.min())
        start, stop = (contact + 1, outer + 1) if direction == "up" else (outer, contact)
        if start >= stop:
            continue
        if extrusion_axis == 0:
            cap[start:stop, j, k] = True
        elif extrusion_axis == 1:
            cap[i, start:stop, k] = True
        else:
            cap[i, j, start:stop] = True
    return cap


def generate_bone_cap_mask(mask, *, axis="x", direction="up", thickness=5, shape="fit", intrusion=None):
    """Generate a one-sided cap mask from a binary bone/contact mask.

    Without ``intrusion`` this uses the surface-projection path. With
    ``intrusion`` it uses the fixed-thickness cap convention:

    * ``thickness`` is the total requested cap thickness.
    * ``intrusion`` is how far the anatomy is allowed to occupy that thickness.
    * generated PMMA cap voxels are removed anywhere they overlap bone.

    The spine workflow calls this with ``shape="anatomy"`` to build a flat
    outside face whose inner side follows nearby body anatomy.
    """
    import numpy as np

    shape = str(shape).lower()
    thickness = int(thickness)
    extrusion_axis = axis_index(axis)
    if intrusion is not None:
        work_bb, _ = _expanded_bbox(mask, extrusion_axis, direction, thickness + int(intrusion))
        cropped = mask[work_bb].astype(bool)
        if shape in {"anatomy", "interactive", "interactive_anatomy"}:
            cap = _generate_axis_aligned_anatomy_cap_from_mask(
                cropped,
                extrusion_axis=extrusion_axis,
                direction=direction,
                thickness=thickness,
                intrusion_depth=int(intrusion),
            )
            output = np.zeros_like(mask, dtype=bool)
            output[work_bb] = cap
            output[mask] = False
            return output
        inferior, superior = _generate_projected_bone_caps_from_mask(
            cropped,
            extrusion_axis=extrusion_axis,
            thickness=thickness,
            intrusion_depth=int(intrusion),
            shape=shape,
        )
        cap = superior if direction == "up" else inferior
        output = np.zeros_like(mask, dtype=bool)
        output[work_bb] = cap
        output[mask] = False
        return largest_connected_component(output)

    work_bb, _ = _expanded_bbox(mask, extrusion_axis, direction, thickness)
    cropped = mask[work_bb].astype(bool)
    surface = _contact_surface(cropped, extrusion_axis, direction)

    if shape == "fit":
        cap = _fit_cap(cropped, _shift_surface(surface, extrusion_axis, direction, thickness), extrusion_axis, direction, thickness)
    else:
        cap = _shaped_contact_cap(cropped, surface, extrusion_axis, direction, thickness, shape)

    output = np.zeros_like(mask, dtype=bool)
    output[work_bb] = cap
    output[mask] = False
    return largest_connected_component(output)


def generate_bone_cap_vtk(
    labelled_vtk_image,
    *,
    label_value,
    axis="x",
    direction="up",
    thickness=5,
    shape="fit",
    intrusion=None,
    output_value=1,
):
    """Generate a VTK PMMA/support cap from a labelled bone mask."""
    import numpy as np

    labels = vtk_image_to_numpy(labelled_vtk_image, processing_order=True)
    cap = generate_bone_cap_mask(
        labels == label_value,
        axis=axis,
        direction=direction,
        thickness=thickness,
        shape=shape,
        intrusion=intrusion,
    )
    return numpy_to_vtk_image(
        (cap.astype(np.uint16) * int(output_value)),
        labelled_vtk_image,
        processing_order=True,
    )


def generate_projected_material_disk_vtk(
    material_vtk_image,
    *,
    surface_vtk_image=None,
    exclusion_vtk_image=None,
    center,
    normal,
    u_axis=None,
    v_axis=None,
    size=(24.0, 24.0),
    shape="square",
    thickness=3.0,
    intrusion=2.0,
    anatomy_constrained=True,
    stable_surface=False,
    stable_surface_fraction=0.55,
    stable_surface_min_area_fraction=0.12,
    stable_surface_max_depth=None,
    close_gaps_3d=False,
    close_gaps_3d_max_gap=2,
    output_value=1,
):
    """Generate a VTK material disk from a physical-space contact plane.

    ``surface_vtk_image`` can be supplied when the anatomical surface should be
    traced from a mask rather than from the binned material image. Nonzero
    voxels in ``exclusion_vtk_image`` are treated as protected anatomy that the
    generated disk must never occupy.
    """
    import numpy as np
    import vtk

    material = vtk_image_to_numpy(material_vtk_image)
    surface = (
        material
        if surface_vtk_image is None
        else vtk_image_to_numpy(surface_vtk_image)
    )
    exclusion = material > 0
    if exclusion_vtk_image is not None:
        exclusion = exclusion | (vtk_image_to_numpy(exclusion_vtk_image) > 0)
    if surface.shape != material.shape:
        raise ValueError("surface_vtk_image must match material_vtk_image shape.")
    if exclusion.shape != material.shape:
        raise ValueError("exclusion_vtk_image must match material_vtk_image shape.")
    extent = material_vtk_image.GetExtent()
    disk = generate_projected_material_disk_mask(
        surface > 0,
        spacing=material_vtk_image.GetSpacing(),
        origin=material_vtk_image.GetOrigin(),
        center=center,
        normal=normal,
        u_axis=u_axis,
        v_axis=v_axis,
        size=size,
        shape=shape,
        thickness=thickness,
        intrusion=intrusion,
        anatomy_constrained=anatomy_constrained,
        stable_surface=stable_surface,
        stable_surface_fraction=stable_surface_fraction,
        stable_surface_min_area_fraction=stable_surface_min_area_fraction,
        stable_surface_max_depth=stable_surface_max_depth,
        close_gaps_3d=close_gaps_3d,
        close_gaps_3d_max_gap=close_gaps_3d_max_gap,
        material_mask=exclusion,
        extent_start=(extent[0], extent[2], extent[4]),
    )
    return numpy_to_vtk_image(
        disk.astype(np.uint16) * int(output_value),
        material_vtk_image,
        vtk_array_type=vtk.VTK_UNSIGNED_SHORT,
    )
