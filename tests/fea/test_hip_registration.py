"""Hip registration must cover both ends and retain superior orientation."""

import numpy as np
import pytest

from ogo.fea import alignment


def test_public_hip_defaults_are_the_locked_ten_mm_recipe():
    from ogo.fea import femur

    assert femur.DEFAULT_FEMUR_CUT_MODE == "greater_trochanter_length"
    assert femur.DEFAULT_FEMUR_GREATER_TROCHANTER_DISTAL_LENGTH_MM == 10
    assert femur.DEFAULT_FEMUR_PISTOIA_CRITICAL_VOLUME_PERCENT == 11.2
    assert femur.DEFAULT_FEMUR_REGISTRATION_LANDMARKS == 40000


def test_stride_sampling_does_not_truncate_the_proximal_end():
    points = np.column_stack([np.zeros(15106), np.zeros(15106), np.arange(15106)])
    sampled = alignment.sample_points(points, max_points=8000, mode="stride")
    assert sampled[:, 2].max() >= 15104
    assert sampled[:, 2].min() <= 1
    assert len(sampled) <= 8000


def test_multistart_recovers_axial_rotation_without_scaling_native_geometry():
    from ogo.fea.hip_registration import estimate_femur_icp

    rng = np.random.default_rng(42)
    points = rng.normal(size=(800, 3)) * [5, 12, 30]
    points[:200] += [15, 20, 35]
    rotation = np.diag([-1.0, -1.0, 1.0])
    fixed = points @ rotation.T + [9, 6, 2]
    fit = estimate_femur_icp(moving_points=points, fixed_points=fixed, iterations=50)
    assert fit["symmetric_distance_mm"] < 1e-8
    assert np.linalg.det(fit["rotation"]) == pytest.approx(1)
    assert fit["rotation"][2, 2] > 0
    assert len(fit["candidates"]) == 4
    again = estimate_femur_icp(moving_points=points, fixed_points=fixed, iterations=50)
    np.testing.assert_allclose(fit["rotation"], again["rotation"])


def test_superior_inversion_cannot_win_candidate_selection(monkeypatch):
    from ogo.fea import hip_registration

    calls = iter([0, 1, 2, 3])

    def fit(**kwargs):
        number = next(calls)
        rotation = np.diag([1., -1., -1.]) if number == 0 else np.eye(3)
        return {"rotation": rotation, "translation": np.zeros(3),
                "iterations": 1, "mean_distance": 0.0}

    monkeypatch.setattr(hip_registration, "estimate_rigid_icp", fit)
    points = np.array([[0., 0., 0.], [1., 0., 1.], [0., 1., 2.]])
    result = hip_registration.estimate_femur_icp(moving_points=points, fixed_points=points)
    assert not result["candidates"][0]["orientation_valid"]
    assert result["selected_start_deg"] != 0
