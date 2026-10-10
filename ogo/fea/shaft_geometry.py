"""Shaft coverage and final-model checks relative to the generated GT disk."""

import numpy as np
from vtk.util.numpy_support import vtk_to_numpy

from ogo.util.vtk_image import vtk_image_to_numpy


def capture_distal_scan_face(vtk_mask, *, labels=None):
    """Keep the native distal voxel-face footprint before cropping or padding.

    Background padding may precede the first occupied slice. Use that slice,
    not the NIfTI bounding box or the separate 120-mm registration crop.
    Coordinates are in Ogo's canonical RAS input frame.
    """
    data = vtk_image_to_numpy(vtk_mask)
    active = np.isin(data, list(labels)) if labels is not None else data != 0
    occupied_z = np.flatnonzero(np.any(active, axis=(0, 1)))
    if not len(occupied_z):
        raise ValueError("Cannot capture the distal scan boundary from an empty mask.")
    first_z = int(occupied_z[0])
    xy = np.argwhere(active[:, :, first_z])
    extent = np.asarray(vtk_mask.GetExtent()).reshape(3, 2)[:, 0]
    origin = np.asarray(vtk_mask.GetOrigin())
    spacing = np.asarray(vtk_mask.GetSpacing())
    z_face = float(origin[2] + (extent[2] + first_z - .5) * spacing[2])
    corners = []
    for offset in ((-.5, -.5), (.5, -.5), (.5, .5), (-.5, .5)):
        corners.append(np.column_stack((
            origin[:2] + (xy + extent[:2] + offset) * spacing[:2],
            np.full(len(xy), z_face),
        )))
    return {"points_xyz": np.unique(np.concatenate(corners), axis=0),
            "native_face_z_mm": z_face,
            "native_distal_slice_index": int(extent[2] + first_z),
            "native_face_voxels": int(len(xy)),
            "native_face_at_image_boundary": first_z == 0}


def aligned_distal_scan_boundary(native_face, reference_to_input, *,
                                 resampling_spacing, output_spacing):
    """Carry the full scan end through the exact rigid ICP and bound resampling.

    ICP maps model/reference coordinates to input coordinates, so invert it.
    The clearance is the projected half-voxel radius of the first resampling
    grid plus half a final Z voxel. This bounds the two nearest-neighbor
    mask resamplings without a cohort-specific distance or percentile.
    """
    if hasattr(reference_to_input, "GetElement"):
        matrix = np.array([[reference_to_input.GetElement(i, j) for j in range(4)]
                           for i in range(4)])
    else:
        matrix = np.asarray(reference_to_input, dtype=float)
    if (matrix.shape != (4, 4) or not np.all(np.isfinite(matrix))
            or not np.allclose(matrix[3], [0, 0, 0, 1])
            or not np.allclose(matrix[:3, :3].T @ matrix[:3, :3], np.eye(3), atol=1e-6)
            or not np.isclose(np.linalg.det(matrix[:3, :3]), 1, atol=1e-6)):
        raise ValueError("The distal scan boundary requires a finite rigid ICP matrix.")
    inverse = np.linalg.inv(matrix)
    if inverse[2, 2] <= 0:
        raise ValueError("The distal input boundary is not inferior after ICP.")
    first_spacing = np.asarray(resampling_spacing, dtype=float)
    final_spacing = np.asarray(output_spacing, dtype=float)
    if (first_spacing.shape != (3,) or final_spacing.shape != (3,)
            or not np.all(np.isfinite([first_spacing, final_spacing]))
            or np.any(first_spacing <= 0) or np.any(final_spacing <= 0)):
        raise ValueError("Boundary resampling spacing must contain three positive values.")
    points = np.asarray(native_face["points_xyz"], dtype=float)
    if (points.ndim != 2 or points.shape[1] != 3 or not len(points)
            or not np.all(np.isfinite(points))):
        raise ValueError("The distal scan boundary has no finite physical face points.")
    aligned_z = points @ inverse[2, :3] + inverse[2, 3]
    clearance = .5 * np.dot(np.abs(inverse[2, :3]), first_spacing) + .5 * final_spacing[2]
    return {**{key: value for key, value in native_face.items() if key != "points_xyz"},
            "source": "full native femur distal voxel face before resampling and rough crop",
            "aligned_face_z_min_mm": float(aligned_z.min()),
            "aligned_face_z_max_mm": float(aligned_z.max()),
            "resampling_clearance_mm": float(clearance),
            "required_flat_face_z_mm": float(aligned_z.max() + clearance)}


def measure_available_shaft(vtk_mask, gt_support_vtk, *, requested_length_mm,
                            distal_boundary=None, labels=None):
    """Measure the full aligned model before cropping, including short scans.

    All coordinates refer to voxel faces. Coverage ends at the first flat
    grid face that clears the transformed original scan end, not the lowest
    corner of an oblique truncation. The GT origin is its connected disk edge.
    """
    requested = float(requested_length_mm)
    if not np.isfinite(requested) or requested <= 0:
        raise ValueError("requested_length_mm must be positive.")
    if distal_boundary is None:
        raise ValueError("The transformed original distal scan boundary is required.")
    required_face = float(distal_boundary["required_flat_face_z_mm"])
    if not np.isfinite(required_face):
        raise ValueError("The distal scan boundary must be finite.")
    spacing = np.asarray(vtk_mask.GetSpacing(), dtype=float)
    origin = np.asarray(vtk_mask.GetOrigin(), dtype=float)
    if (gt_support_vtk.GetExtent() != vtk_mask.GetExtent()
            or not np.allclose(gt_support_vtk.GetOrigin(), origin)
            or not np.allclose(gt_support_vtk.GetSpacing(), spacing)):
        raise ValueError("GT support and model must share the same grid.")
    data = vtk_image_to_numpy(vtk_mask)
    active = np.isin(data, list(labels)) if labels is not None else data != 0
    if not np.any(active):
        raise ValueError("Cannot measure shaft length from an empty model.")
    disk = (vtk_image_to_numpy(gt_support_vtk) != 0) & active
    if not np.any(disk):
        raise ValueError("Cannot define shaft length from an empty GT support.")
    gt_index = int(np.argwhere(disk)[:, 2].min())
    distal_index = int(np.argwhere(active)[:, 2].min())
    extent_z = vtk_mask.GetExtent()[4]
    gt_face = float(origin[2] + (extent_z + gt_index - 0.5) * spacing[2])
    distal_face = float(origin[2] + (extent_z + distal_index - 0.5) * spacing[2])
    first_grid_face = origin[2] + (extent_z - .5) * spacing[2]
    safe_index = int(np.ceil((required_face - first_grid_face) / spacing[2] - 1e-9))
    safe_face = float(max(distal_face, first_grid_face + safe_index * spacing[2]))
    available = max(0.0, gt_face - safe_face)
    eligible = available + 1e-6 >= requested
    return {"type": "ogo.femur.shaft_geometry", "version": 3,
            "definition": "generated GT disk minimum Z face minus complete-section safe flat face Z",
            "available_shaft_length_mm": available,
            "raw_available_shaft_length_mm": gt_face - distal_face,
            "requested_shaft_length_mm": requested,
            "measured_shaft_length_mm": None,
            "eligible_for_requested_length": bool(eligible),
            "status": "measured" if eligible else "too_short",
            "full_model_distal_face_z_mm": distal_face,
            "safe_distal_face_z_mm": safe_face,
            "oblique_trim_loss_mm": safe_face - distal_face,
            "distal_scan_boundary": dict(distal_boundary),
            "complete_distal_section_verified": False,
            "greater_trochanter_support_distal_edge_z": gt_face,
            "distal_length_origin_z": gt_face,
            "available_below_gt_support_distal_edge_mm": available}


def verify_model_shaft(model, metadata):
    """Check the actual FE node sets before solving; return a verified sidecar."""
    points = vtk_to_numpy(model.GetPoints().GetData())
    coordinates = {}
    for key, name in (("gt", "Greater_Trochanter_PMMA_Nodes"),
                      ("distal", "Distal_Femur_Nodes")):
        ids = model.GetNodeSet(name)
        if ids is None or ids.GetNumberOfTuples() == 0:
            raise ValueError("Cannot verify shaft length: missing %s." % name)
        coordinates[key] = points[vtk_to_numpy(ids).astype(int)]
    gt_z = float(coordinates["gt"][:, 2].min())
    distal_z = float(np.median(coordinates["distal"][:, 2]))
    span = float(np.ptp(coordinates["distal"][:, 2]))
    length = gt_z - distal_z
    retained = float(metadata["retained_length_mm"])
    if (abs(length - retained) > 0.01 or span > 0.01
            or abs(gt_z - metadata["distal_length_origin_z"]) > 0.01
            or abs(distal_z - metadata["cut_z_mm"]) > 0.01):
        raise ValueError("Saved shaft BC nodes do not match the generated GT-disk crop.")
    if distal_z < float(metadata["safe_distal_face_z_mm"]) - .01:
        raise ValueError("Saved shaft cut does not clear the original distal scan boundary.")
    return dict(metadata, status="verified", measured_shaft_length_mm=length,
                complete_distal_section_verified=True,
                gt_bc_distal_z_mm=gt_z, distal_bc_median_z_mm=distal_z,
                distal_bc_z_span_mm=span, gt_bc_nodes=len(coordinates["gt"]),
                distal_bc_nodes=len(coordinates["distal"]))
