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


def _projected_contact_within_depth(mask, extrusion_axis, direction, depth):
    """Project only columns with bone near the flat fixture contact plane."""
    import numpy as np

    hit = np.where(mask)[extrusion_axis]
    contact = int(hit.max() if direction == "up" else hit.min())
    depth = max(int(depth), 1)
    if direction == "up":
        support_slice = slice(max(contact - depth + 1, 0), contact + 1)
    elif direction == "down":
        support_slice = slice(contact, min(contact + depth, mask.shape[extrusion_axis]))
    else:
        raise ValueError("direction must be 'up' or 'down'.")

    if extrusion_axis == 0:
        return mask[support_slice, :, :].any(axis=0)
    if extrusion_axis == 1:
        return mask[:, support_slice, :].any(axis=1)
    return mask[:, :, support_slice].any(axis=2)


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


def _full_fixture_cap(
    mask,
    surface,
    extrusion_axis,
    direction,
    thickness,
    shape,
    intrusion=0,
    crop_to_contact=False,
):
    import numpy as np

    footprint = _projected_footprint(surface, extrusion_axis, shape)
    if crop_to_contact:
        footprint &= _projected_contact_within_depth(mask, extrusion_axis, direction, intrusion)
    cap = np.zeros_like(mask, dtype=bool)
    hit = np.where(mask)[extrusion_axis]
    contact = int(hit.min() if direction == "down" else hit.max())
    intrusion = max(int(intrusion), 0)
    if direction == "up":
        slab = slice(
            max(contact - intrusion + 1, 0),
            min(contact + thickness + 1, mask.shape[extrusion_axis]),
        )
    else:
        slab = slice(
            max(contact - thickness, 0),
            min(contact + intrusion, mask.shape[extrusion_axis]),
        )
    if extrusion_axis == 0:
        cap[slab, :, :] = footprint
    elif extrusion_axis == 1:
        cap[:, slab, :] = footprint[:, None, :]
    else:
        cap[:, :, slab] = footprint[:, :, None]
    return cap


def generate_bone_cap_mask(mask, *, axis="x", direction="up", thickness=5, shape="fit"):
    """Generate a one-sided cap mask from a binary bone/contact mask."""
    import numpy as np

    shape = str(shape).lower()
    thickness = int(thickness)
    extrusion_axis = axis_index(axis)
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


def generate_fixture_cap_mask(
    mask,
    *,
    axis="x",
    direction="up",
    thickness=5,
    shape="box",
    intrusion=0,
    crop_to_contact=False,
):
    """Generate a full geometric fixture cap from a binary ROI/contact mask.

    Unlike a bone-fit cap, fixture caps keep their full square/round footprint.
    They may overlap the bone mask; the later image-combine step preserves bone
    material IDs where PMMA and bone overlap.
    """
    import numpy as np

    shape = str(shape).lower()
    if shape not in {"box", "round"}:
        raise ValueError("fixture cap shape must be 'box' or 'round'.")
    thickness = int(thickness)
    extrusion_axis = axis_index(axis)
    work_bb, _ = _expanded_bbox(mask, extrusion_axis, direction, thickness)
    cropped = mask[work_bb].astype(bool)
    surface = _contact_surface(cropped, extrusion_axis, direction)
    cap = _full_fixture_cap(
        cropped,
        surface,
        extrusion_axis,
        direction,
        thickness,
        shape,
        intrusion=intrusion,
        crop_to_contact=crop_to_contact,
    )

    output = np.zeros_like(mask, dtype=bool)
    output[work_bb] = cap
    return largest_connected_component(output)


def label_foreground_in_bounds_vtk(vtk_image, bounds, *, label_value=1):
    """Label nonzero voxels inside physical bounds for contact-cap generation."""
    import numpy as np
    import vtk

    dims = vtk_image.GetDimensions()
    extent = vtk_image.GetExtent()
    spacing = vtk_image.GetSpacing()
    origin = vtk_image.GetOrigin()
    data = vtk_image_to_numpy(vtk_image)

    x = origin[0] + (np.arange(dims[0]) + extent[0]) * spacing[0]
    y = origin[1] + (np.arange(dims[1]) + extent[2]) * spacing[1]
    z = origin[2] + (np.arange(dims[2]) + extent[4]) * spacing[2]
    in_bounds = (
        (x[:, None, None] >= bounds[0])
        & (x[:, None, None] <= bounds[1])
        & (y[None, :, None] >= bounds[2])
        & (y[None, :, None] <= bounds[3])
        & (z[None, None, :] >= bounds[4])
        & (z[None, None, :] <= bounds[5])
    )

    output_data = np.where((data != 0) & in_bounds, label_value, 0).astype(np.uint16)
    return numpy_to_vtk_image(output_data, vtk_image, vtk_array_type=vtk.VTK_UNSIGNED_SHORT)


def generate_bone_cap_vtk(
    labelled_vtk_image,
    *,
    label_value,
    axis="x",
    direction="up",
    thickness=5,
    shape="fit",
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
    )
    return numpy_to_vtk_image(
        (cap.astype(np.uint16) * int(output_value)),
        labelled_vtk_image,
        processing_order=True,
    )


def generate_fixture_cap_vtk(
    labelled_vtk_image,
    *,
    label_value,
    axis="x",
    direction="up",
    thickness=5,
    shape="box",
    intrusion=0,
    crop_to_contact=False,
    output_value=1,
):
    """Generate a full geometric VTK fixture cap from a labelled ROI mask."""
    import numpy as np

    labels = vtk_image_to_numpy(labelled_vtk_image, processing_order=True)
    cap = generate_fixture_cap_mask(
        labels == label_value,
        axis=axis,
        direction=direction,
        thickness=thickness,
        shape=shape,
        intrusion=intrusion,
        crop_to_contact=crop_to_contact,
    )
    return numpy_to_vtk_image(
        (cap.astype(np.uint16) * int(output_value)),
        labelled_vtk_image,
        processing_order=True,
    )
