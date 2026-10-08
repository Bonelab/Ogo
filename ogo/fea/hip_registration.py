"""Deterministic rigid hip registration in the canonical superior-Z frame.

Reference scaling, when used, is performed by the caller. This fit is rigid:
it never scales the participant's anatomy or its solved voxel dimensions.
"""

import numpy as np
from scipy.spatial import cKDTree

from ogo.fea.alignment import estimate_rigid_icp


def estimate_femur_icp(*, moving_points, fixed_points, iterations=50):
    """Select the best non-inverted fit from four axial starting rotations.

    Bidirectional mean surface distance avoids choosing a fit that matches only
    one part of the femur. Positive superior-axis projection excludes SI flips;
    it does not replace visual alignment QC. KD-tree queries use one worker so
    cohort-level parallelism remains bounded.
    """
    moving = np.asarray(moving_points, dtype=float)
    fixed = np.asarray(fixed_points, dtype=float)
    if (moving.ndim != 2 or fixed.ndim != 2 or moving.shape[1] != 3
            or fixed.shape[1] != 3 or min(len(moving), len(fixed)) < 3
            or not np.all(np.isfinite(moving)) or not np.all(np.isfinite(fixed))):
        raise ValueError("ICP requires finite (n, 3) surfaces with at least three points.")
    fixed_tree = cKDTree(fixed)
    trials = []
    best = None
    for angle in (0, 90, 180, 270):
        theta = np.deg2rad(angle)
        rotation = np.array([[np.cos(theta), -np.sin(theta), 0],
                             [np.sin(theta), np.cos(theta), 0], [0, 0, 1.]])
        fit = estimate_rigid_icp(
            moving_points=moving, fixed_points=fixed, iterations=iterations,
            initial_transform={"rotation": rotation,
                               "translation": fixed.mean(axis=0) - moving.mean(axis=0) @ rotation.T},
            nearest_workers=1,
        )
        transformed = moving @ fit["rotation"].T + fit["translation"]
        forward = float(fixed_tree.query(transformed, workers=1)[0].mean())
        reverse = float(cKDTree(transformed).query(fixed, workers=1)[0].mean())
        score = (forward + reverse) / 2
        superior_projection = float(fit["rotation"][2, 2])
        valid = bool(superior_projection > 0 and np.isfinite(score))
        trials.append({"start_deg": angle, "forward_distance_mm": forward,
                       "reverse_distance_mm": reverse, "symmetric_distance_mm": score,
                       "superior_axis_projection": superior_projection,
                       "orientation_valid": valid, "iterations": int(fit["iterations"])})
        if valid and (best is None or score < best["symmetric_distance_mm"]):
            best = dict(fit, mean_distance=forward, symmetric_distance_mm=score,
                        selected_start_deg=angle, superior_axis_projection=superior_projection)
    if best is None:
        raise ValueError("No femur ICP candidate preserved superior orientation.")
    return dict(best, strategy="axial_multistart_rigid", candidates=trials,
                reference_points=len(moving), sample_points=len(fixed), visual_qc_status="pending")
