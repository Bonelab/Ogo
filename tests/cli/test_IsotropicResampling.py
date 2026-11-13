"""Analytical voxel and geometry checks for the existing resampling command."""

import numpy as np
import pytest
import SimpleITK as sitk


def test_mask_samples_and_physical_geometry(run_installed_command, tmp_path):
    data = np.zeros((4, 4, 6), dtype=np.uint8)
    data[:, 1:3, 1:5] = 7
    image = sitk.GetImageFromArray(data)
    image.SetSpacing((1.0, 1.0, 2.0))
    image.SetOrigin((11.0, -7.0, 3.0))
    image.SetDirection((0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0))
    source, target = tmp_path / "mask.nii.gz", tmp_path / "resampled.nii.gz"
    sitk.WriteImage(image, str(source))
    original = source.read_bytes()

    result = run_installed_command(
        "ogoIsotropicResampling", source, target,
        "--iso_resolution", "2", "--image_type", "mask",
    )

    assert result.returncode == 0, result.stdout + result.stderr
    output = sitk.ReadImage(str(target))
    assert output.GetSize() == (3, 2, 4)
    assert output.GetSpacing() == pytest.approx((2.0, 2.0, 2.0))
    assert output.GetOrigin() == pytest.approx(image.GetOrigin())
    assert output.GetDirection() == pytest.approx(image.GetDirection())
    assert output.GetPixelID() == sitk.sitkUInt8
    # Target positions select every second x/y voxel; z spacing is unchanged.
    np.testing.assert_array_equal(sitk.GetArrayFromImage(output), data[:, ::2, ::2])
    assert source.read_bytes() == original


def test_constant_ct_at_existing_resolution(run_installed_command, tmp_path):
    image = sitk.GetImageFromArray(np.full((6, 6, 6), 123.25, dtype=np.float32))
    image.SetOrigin((4.0, 5.0, -6.0))
    source, target = tmp_path / "ct.nii.gz", tmp_path / "ct-resampled.nii.gz"
    sitk.WriteImage(image, str(source))
    original = source.read_bytes()

    result = run_installed_command(
        "ogoIsotropicResampling", source, target,
        "--iso_resolution", "1", "--image_type", "ct",
    )

    assert result.returncode == 0, result.stdout + result.stderr
    output = sitk.ReadImage(str(target))
    assert output.GetSize() == image.GetSize()
    assert output.GetSpacing() == pytest.approx(image.GetSpacing())
    assert output.GetOrigin() == pytest.approx(image.GetOrigin())
    assert output.GetDirection() == pytest.approx(image.GetDirection())
    assert output.GetPixelID() == sitk.sitkFloat32
    np.testing.assert_allclose(sitk.GetArrayFromImage(output), 123.25, rtol=0, atol=1e-5)
    assert source.read_bytes() == original
