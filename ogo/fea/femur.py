"""Femur-specific helpers for sideways-fall FE model generation."""


LEFT_FEMUR = 1
RIGHT_FEMUR = 2


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
