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

    # Exercise the installed pilot command.
    result = run_installed_command(
        "ogo", "resample-image", source, target,
        "--iso-resolution", "2", "--image-type", "mask",
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout == ""
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
        "ogo", "resample-image", source, target,
        "--iso-resolution", "1", "--image-type", "ct",
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout == ""
    output = sitk.ReadImage(str(target))
    assert output.GetSize() == image.GetSize()
    assert output.GetSpacing() == pytest.approx(image.GetSpacing())
    assert output.GetOrigin() == pytest.approx(image.GetOrigin())
    assert output.GetDirection() == pytest.approx(image.GetDirection())
    assert output.GetPixelID() == sitk.sitkFloat32
    np.testing.assert_allclose(sitk.GetArrayFromImage(output), 123.25, rtol=0, atol=1e-5)
    assert source.read_bytes() == original


@pytest.fixture
def source_image(tmp_path):
    source = tmp_path / "source.nii.gz"
    image = sitk.GetImageFromArray(np.full((6, 6, 6), 17, dtype=np.uint8))
    sitk.WriteImage(image, str(source))
    return source


@pytest.mark.parametrize("resolution", [0, -1, float("nan"), float("inf"), -float("inf")])
def test_invalid_resolution_does_not_write(run_installed_command, source_image, tmp_path, resolution):
    from ogo.cli.IsotropicResampling import resample_image_file
    target = tmp_path / "output.nii.gz"
    with pytest.raises(ValueError, match="resolution"):
        resample_image_file(source_image, target, iso_resolution=resolution)
    result = run_installed_command("ogo", "resample-image", source_image, target,
                                   f"--iso-resolution={resolution}")
    assert result.returncode == 1
    assert result.stdout == ""
    assert "[ERROR]" in result.stderr
    assert not target.exists()


@pytest.mark.parametrize("mode", ["c", "", "invalid"])
def test_invalid_direct_mode_is_rejected(source_image, tmp_path, mode):
    from ogo.cli.IsotropicResampling import resample_image_file
    from ogo.util.Helper import isotropicResampling
    target = tmp_path / "output.nii.gz"
    with pytest.raises(ValueError, match="image type"):
        resample_image_file(source_image, target, image_type=mode)
    with pytest.raises(ValueError, match="image type"):
        isotropicResampling(None, 1.0, mode)
    assert not target.exists()


def test_existing_output_requires_permission_in_python_and_cli(
    run_installed_command, source_image, tmp_path, monkeypatch,
):
    from ogo.cli.IsotropicResampling import resample_image_file
    target = tmp_path / "output.nii.gz"
    target.write_bytes(b"existing result")
    original = source_image.read_bytes()
    def no_input(*args):
        pytest.fail("Resampling must never prompt")
    monkeypatch.setattr("builtins.input", no_input)
    with pytest.raises(FileExistsError, match="overwrite"):
        resample_image_file(source_image, target)
    result = run_installed_command("ogo", "resample-image", source_image, target)
    assert result.returncode == 1
    assert result.stdout == ""
    assert "--overwrite" in result.stderr
    assert target.read_bytes() == b"existing result"
    assert source_image.read_bytes() == original


@pytest.mark.parametrize("interface", ["python", "cli"])
def test_explicit_overwrite_replaces_output(run_installed_command, source_image, tmp_path, interface):
    from ogo.cli.IsotropicResampling import resample_image_file
    target = tmp_path / "output.nii.gz"
    target.write_bytes(b"existing result")
    if interface == "python":
        assert resample_image_file(source_image, target, image_type="mask", overwrite=True) == str(target)
    else:
        result = run_installed_command("ogo", "--quiet", "resample-image", source_image,
                                       target, "--image-type", "mask", "-o")
        assert result.returncode == 0, result.stderr
        assert result.stdout == result.stderr == ""
    output = sitk.ReadImage(str(target))
    assert output.GetPixelID() == sitk.sitkUInt8
    np.testing.assert_array_equal(sitk.GetArrayFromImage(output), 17)


@pytest.mark.parametrize("failure", ["missing", "corrupt", "write"])
def test_file_failures_return_one(run_installed_command, source_image, tmp_path, failure):
    source = source_image
    target = tmp_path / "output.nii.gz"
    if failure == "missing":
        source = tmp_path / "absent.nii.gz"
    elif failure == "corrupt":
        source = tmp_path / "corrupt.nii.gz"
        source.write_bytes(b"invalid image")
    else:
        target = tmp_path / "absent-directory" / "output.nii.gz"
    original = source.read_bytes() if source.is_file() else None
    result = run_installed_command("ogo", "resample-image", source, target)
    assert result.returncode == 1
    assert result.stdout == ""
    assert "[ERROR]" in result.stderr
    assert not target.exists()
    if original is not None:
        assert source.read_bytes() == original


@pytest.mark.parametrize("options", [["--iso_resolution", "1"], ["--image_type", "mask"],
                                     ["--image-type", "invalid"]])
def test_old_or_invalid_options_are_usage_errors(run_installed_command, source_image, tmp_path, options):
    target = tmp_path / "output.nii.gz"
    result = run_installed_command("ogo", "resample-image", source_image, target, *options)
    assert result.returncode == 2
    assert result.stdout == ""
    assert not target.exists()


@pytest.mark.parametrize("kind", ["directory", "dangling-symlink"])
def test_other_existing_destinations_are_refused(source_image, tmp_path, kind):
    from ogo.cli.IsotropicResampling import resample_image_file
    target = tmp_path / "output.nii.gz"
    if kind == "directory":
        target.mkdir()
    else:
        target.symlink_to(tmp_path / "absent.nii.gz")
    with pytest.raises(FileExistsError):
        resample_image_file(source_image, target)
    assert target.is_dir() if kind == "directory" else target.is_symlink()
