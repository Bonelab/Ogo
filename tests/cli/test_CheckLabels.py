"""Scientific label outcomes before changing the command interface."""

import numpy as np
import pytest
import SimpleITK as sitk


@pytest.mark.parametrize("options,status,report", [
    ([], 0, "Labels found: [0, 1, 2, 91]"),
    (["--require", "1", "2"], 0, "All required labels are present."),
    (["--require", "1", "9"], 2, "Missing required labels: [9]"),
    (
        ["--include-description"], 0,
        "Labels found: 0 (Clear Label), 1 (Femur Right), 2 (Femur Left), 91 (Adipose)",
    ),
    (
        ["--require", "1", "9", "--include-description"], 2,
        "Missing required labels: 9 (L2)",
    ),
])
def test_label_outcomes(run_installed_command, tmp_path, options, status, report):
    data = np.array([[[0, 1], [2, 91]], [[1, 0], [2, 0]]], dtype=np.uint8)
    source = tmp_path / "segmentation.nii.gz"
    sitk.WriteImage(sitk.GetImageFromArray(data), str(source))
    original = source.read_bytes()

    result = run_installed_command("ogoCheckLabels", source, *options)

    assert result.returncode == status, result.stdout + result.stderr
    # Streams and elapsed/severity prefixes will intentionally change later.
    assert report in result.stdout + result.stderr
    assert source.read_bytes() == original
