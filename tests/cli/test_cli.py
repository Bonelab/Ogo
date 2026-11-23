"""Exercise the installed root dispatcher and its deliberately small catalog."""

import sys
from types import SimpleNamespace

import pytest

from ogo import __version__
from ogo.cli import cli


@pytest.mark.parametrize("arguments", [[], ["--help"]])
def test_root_help_lists_the_pilot(run_installed_command, arguments):
    result = run_installed_command("ogo", *arguments)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "print-labels" in result.stdout
    assert "Print standard label definitions" in result.stdout
    assert result.stderr == ""


def test_version_reports_the_installed_version(run_installed_command):
    result = run_installed_command("ogo", "--version")
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout.strip() == f"ogo {__version__}"
    assert result.stderr == ""


def test_selected_command_handles_its_own_help(run_installed_command):
    result = run_installed_command("ogo", "print-labels", "--help")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "usage: ogo print-labels" in result.stdout
    assert "--echo {labels,colors16,colors128}" in result.stdout
    assert result.stderr == ""


@pytest.mark.parametrize("arguments", [
    ["unknown-command"],
    ["check-labels"],
    ["resample-image"],
    ["--unknown-option"],
    ["print-labels", "--unknown-option"],
    ["print-labels", "--echo", "invalid"],
    ["print-labels", "--echo"],
])
def test_invalid_requests_are_usage_errors(run_installed_command, arguments):
    result = run_installed_command("ogo", *arguments)
    assert result.returncode == 2
    assert result.stdout == ""
    assert "usage:" in result.stderr
    assert "error:" in result.stderr


def test_dispatch_preserves_arguments_and_failure_status(monkeypatch):
    received = []

    def handler(argv):
        received.extend(argv)
        return 7

    monkeypatch.setitem(cli.COMMANDS, "example", ("ogo.cli.example", "Example for testing"))
    monkeypatch.setitem(sys.modules, "ogo.cli.example", SimpleNamespace(main=handler))
    assert cli.main(["example", "input.nii.gz", "--option", "value"]) == 7
    assert received == ["input.nii.gz", "--option", "value"]


@pytest.mark.parametrize("arguments", [["--help"], ["print-labels", "--echo", "labels"]])
def test_help_and_labels_work_when_scientific_imports_are_unavailable(run_python, arguments):
    code = """
import importlib.abc
import sys

class BlockScientificImports(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in {
            'numpy', 'scipy', 'SimpleITK', 'vtk', 'vtkmodules', 'vtkbone'
        } or fullname == 'ogo.util.Helper':
            raise ModuleNotFoundError('Scientific import blocked by test: ' + fullname)

sys.meta_path.insert(0, BlockScientificImports())
from ogo.cli.cli import main
raise SystemExit(main(sys.argv[1:]))
"""
    result = run_python("-c", code, *arguments)
    assert result.returncode == 0, result.stdout + result.stderr
    expected = "print-labels" if arguments == ["--help"] else '"Femur Right"'
    assert expected in result.stdout
    assert result.stderr == ""
