"""Verify the study example never substitutes neighboring vertebral levels."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip('SimpleITK')
source = Path(__file__).resolve().parents[2] / 'examples/fea/cohort_recovery/recover_spines.py'
spec = importlib.util.spec_from_file_location('recover_spines', source)
recovery = importlib.util.module_from_spec(spec)
spec.loader.exec_module(recovery)


def test_only_l1_body_and_process_are_retained():
    level = np.array([20, 20, 19, 21])
    parts = np.array([2, 1, 2, 1])
    np.testing.assert_array_equal(recovery.isolate_l1(level, parts, np.eye(4), np.eye(4)), [2, 1, 0, 0])


def test_mismatched_affine_is_rejected():
    shifted = np.eye(4)
    shifted[0, 3] = 2
    with pytest.raises(ValueError, match='grids'):
        recovery.isolate_l1(np.array([20, 20]), np.array([2, 1]), np.eye(4), shifted)


def test_missing_l1_and_unknown_part_are_rejected():
    for level, parts in [([19, 21], [2, 1]), ([20, 20], [2, 3])]:
        with pytest.raises(ValueError):
            recovery.isolate_l1(np.array(level), np.array(parts), np.eye(4), np.eye(4))
