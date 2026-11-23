"""Preserved label operations and the installed CheckLabels interface."""

import numpy as np
import pytest
import SimpleITK as sitk

from ogo.cli.CheckLabels import add_label_descriptions, find_missing_labels, read_labels


@pytest.fixture
def segmentation(tmp_path):
    data = np.array([[[0, 1], [2, 91]], [[1, 0], [2, 0]]], dtype=np.uint8)
    source = tmp_path / "segmentation.nii.gz"
    sitk.WriteImage(sitk.GetImageFromArray(data), str(source))
    return source


@pytest.mark.parametrize("options,status,report", [
    ([], 0, "Labels found: [0, 1, 2, 91]"),
    (["--require"], 0, "Labels found: [0, 1, 2, 91]"),
    (["--require", "1", "2"], 0, "All required labels are present."),
    (["--require", "0", "2", "2"], 0, "All required labels are present."),
    (["--require", "1", "9"], 2, "Missing required labels: [9]"),
    (["--require", "9", "2", "8", "9"], 2, "Missing required labels: [9, 8, 9]"),
    (
        ["--include-description"], 0,
        "Labels found: 0 (Clear Label), 1 (Femur Right), 2 (Femur Left), 91 (Adipose)",
    ),
    (
        ["--require", "1", "9", "--include-description"], 2,
        "Missing required labels: 9 (L2)",
    ),
    (
        ["--require", "999", "9", "999", "--include-description"], 2,
        "Missing required labels: 999 (unknown label), 9 (L2), 999 (unknown label)",
    ),
])
def test_label_outcomes(run_installed_command, segmentation, options, status, report):
    original = segmentation.read_bytes()

    result = run_installed_command("ogo", "check-labels", segmentation, *options)

    assert result.returncode == status, result.stdout + result.stderr
    if status == 0:
        assert result.stdout == report + "\n"
        assert report not in result.stderr
    else:
        assert result.stdout == ""
        assert "[ERROR] " + report in result.stderr
    assert "[INFO] Reading input image:" in result.stderr
    assert segmentation.read_bytes() == original


def test_read_labels_returns_inventory_without_output(segmentation, capsys):
    original = segmentation.read_bytes()
    assert read_labels(segmentation) == [0, 1, 2, 91]
    assert capsys.readouterr() == ("", "")
    assert segmentation.read_bytes() == original


def test_read_labels_preserves_uint8_conversion(tmp_path):
    data = np.array([[[0, 1, 255, 256, 257, 511]]], dtype=np.uint16)
    source = tmp_path / "uint16.nii.gz"
    sitk.WriteImage(sitk.GetImageFromArray(data), str(source))
    original = source.read_bytes()

    # Retain historical conversion rather than silently supporting wider IDs.
    assert read_labels(source) == [0, 1, 255]
    assert source.read_bytes() == original


@pytest.mark.parametrize("present,required,expected", [
    ([0, 1, 2], [9, 2, 8, 9], [9, 8, 9]),
    ([0, 1, 2], [0, 2, 2], []),
    ([0, 1, 2], [], []),
    ([], [9, 9], [9, 9]),
])
def test_missing_labels_are_returned_without_side_effects(present, required, expected, capsys):
    original_present, original_required = present[:], required[:]
    assert find_missing_labels(present, required) == expected
    assert present == original_present
    assert required == original_required
    assert capsys.readouterr() == ("", "")


def test_description_fallback_preserves_order_and_duplicates():
    assert add_label_descriptions([999, 0, 999]) == [
        "999 (unknown label)", "0 (Clear Label)", "999 (unknown label)",
    ]


@pytest.mark.parametrize("kind,error_type", [
    ("missing", FileNotFoundError),
    ("directory", FileNotFoundError),
    ("corrupt", RuntimeError),
])
def test_input_failures_raise_directly_and_return_status_one(
    run_installed_command, tmp_path, kind, error_type, capsys,
):
    source = tmp_path / "invalid.nii.gz"
    if kind == "directory":
        source.mkdir()
    elif kind == "corrupt":
        source.write_bytes(b"not an image")

    with pytest.raises(error_type):
        read_labels(source)
    assert capsys.readouterr() == ("", "")

    result = run_installed_command("ogo", "check-labels", source)
    assert result.returncode == 1, result.stdout + result.stderr
    assert result.stdout == ""
    assert "[ERROR]" in result.stderr
    assert str(source) in result.stderr
    assert "Traceback" not in result.stderr
    if kind == "corrupt":
        assert source.read_bytes() == b"not an image"


def test_unusual_suffix_warns_but_readable_image_succeeds(run_installed_command, tmp_path):
    source = tmp_path / "segmentation.mha"
    data = np.array([[[0, 2, 1]]], dtype=np.uint8)
    sitk.WriteImage(sitk.GetImageFromArray(data), str(source))
    original = source.read_bytes()

    result = run_installed_command("ogo", "--quiet", "check-labels", source)
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout == "Labels found: [0, 1, 2]\n"
    assert "[WARNING] Input image does not end with .nii or .nii.gz" in result.stderr
    assert "[INFO]" not in result.stderr
    assert source.read_bytes() == original


@pytest.mark.parametrize("options,status,stdout,diagnostic", [
    ([], 0, "Labels found: [0, 1, 2, 91]\n", None),
    (["--require", "1", "2"], 0, "All required labels are present.\n", None),
    (["--require", "9"], 2, "", "Missing required labels: [9]"),
])
def test_quiet_retains_results_and_failures(
    run_installed_command, segmentation, options, status, stdout, diagnostic,
):
    result = run_installed_command("ogo", "--quiet", "check-labels", segmentation, *options)
    assert result.returncode == status, result.stdout + result.stderr
    assert result.stdout == stdout
    if diagnostic is None:
        assert result.stderr == ""
    else:
        assert "[ERROR] " + diagnostic in result.stderr
        assert "[INFO]" not in result.stderr


def test_verbose_explains_an_unreadable_image(run_installed_command, tmp_path):
    source = tmp_path / "corrupt.nii"
    source.write_bytes(b"not an image")
    result = run_installed_command("ogo", "--verbose", "check-labels", source)
    assert result.returncode == 1
    assert result.stdout == ""
    assert "[ERROR] Cannot read input image:" in result.stderr
    assert "[DEBUG] Image read failed" in result.stderr
    assert "Traceback" in result.stderr


def test_check_labels_does_not_import_helper_or_vtk(run_python, segmentation):
    code = """
import importlib.abc
import sys

class BlockUnrelatedImports(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in {'vtk', 'vtkmodules', 'vtkbone', 'scipy'} or fullname == 'ogo.util.Helper':
            raise ModuleNotFoundError('Unrelated import blocked: ' + fullname)

sys.meta_path.insert(0, BlockUnrelatedImports())
from ogo.cli.cli import main
raise SystemExit(main(sys.argv[1:]))
"""
    result = run_python("-c", code, "--quiet", "check-labels", segmentation)
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout == "Labels found: [0, 1, 2, 91]\n"
    assert result.stderr == ""
