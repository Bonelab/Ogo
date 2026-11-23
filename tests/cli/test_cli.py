"""Exercise the installed root dispatcher and its deliberately small catalog."""

import io
import logging
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
    assert "check-labels" in result.stdout
    assert "resample-image" not in result.stdout
    assert result.stderr == ""


def test_version_reports_the_installed_version(run_installed_command):
    result = run_installed_command("ogo", "--version")
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout.strip() == f"ogo {__version__}"
    assert result.stderr == ""


@pytest.mark.parametrize("command,option", [
    ("print-labels", "--echo {labels,colors16,colors128}"),
    ("check-labels", "--include-description"),
])
def test_selected_command_handles_its_own_help(run_installed_command, command, option):
    result = run_installed_command("ogo", command, "--help")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "usage: ogo " + command in result.stdout
    assert option in result.stdout
    assert result.stderr == ""


@pytest.mark.parametrize("arguments", [
    ["unknown-command"],
    ["check-labels"],
    ["resample-image"],
    ["--unknown-option"],
    ["--quiet", "--verbose", "print-labels"],
    ["print-labels", "--quiet"],
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


@pytest.mark.parametrize("options,visible", [
    ([], {"INFO", "WARNING", "ERROR"}),
    (["--verbose"], {"DEBUG", "INFO", "WARNING", "ERROR"}),
    (["--quiet"], {"WARNING", "ERROR"}),
])
def test_verbosity_controls_diagnostics_only(monkeypatch, capsys, options, visible):
    def handler(argv):
        logger = logging.getLogger("ogo.cli.example")
        for level in ("DEBUG", "INFO", "WARNING", "ERROR"):
            logger.log(getattr(logging, level), "Diagnostic %s", level)
        print("Requested result")
        return 0

    monkeypatch.setitem(cli.COMMANDS, "example", ("ogo.cli.example", "Example for testing"))
    monkeypatch.setitem(sys.modules, "ogo.cli.example", SimpleNamespace(main=handler))
    assert cli.main([*options, "example"]) == 0
    output = capsys.readouterr()
    assert output.out == "Requested result\n"
    for level in ("DEBUG", "INFO", "WARNING", "ERROR"):
        assert (f"[{level}] Diagnostic {level}" in output.err) == (level in visible)


def test_logging_timing_continuations_and_repeated_calls_preserve_handlers(monkeypatch, capsys):
    logger = logging.getLogger("ogo")
    application_output = io.StringIO()
    application_handler = logging.StreamHandler(application_output)
    monkeypatch.setattr(logger, "handlers", [*logger.handlers, application_handler])
    previous = (logger.handlers[:], logger.level, logger.propagate)
    root = logging.getLogger()
    previous_root = (root.handlers[:], root.level)
    ticks = iter([100.0, 100.25, 200.0, 200.5])
    monkeypatch.setattr(cli.time, "monotonic", lambda: next(ticks))

    def handler(argv):
        logging.getLogger("ogo.cli.example").info("Reading input\nSize: 1")
        return 0

    monkeypatch.setitem(cli.COMMANDS, "example", ("ogo.cli.example", "Example for testing"))
    monkeypatch.setitem(sys.modules, "ogo.cli.example", SimpleNamespace(main=handler))
    for _ in range(2):
        assert cli.main(["example"]) == 0
        assert (logger.handlers, logger.level, logger.propagate) == previous
        assert (root.handlers, root.level) == previous_root

    lines = capsys.readouterr().err.splitlines()
    assert lines == [
        "    0.25s [INFO] Reading input", " " * 17 + "Size: 1",
        "    0.50s [INFO] Reading input", " " * 17 + "Size: 1",
    ]
    assert application_output.getvalue() == "Reading input\nSize: 1\n" * 2


@pytest.mark.parametrize("error", [RuntimeError("unexpected failure"), SystemExit(2)])
def test_logging_is_restored_when_a_handler_raises(monkeypatch, error):
    logger = logging.getLogger("ogo")
    previous = (logger.handlers[:], logger.level, logger.propagate)

    def handler(argv):
        raise error

    monkeypatch.setitem(cli.COMMANDS, "example", ("ogo.cli.example", "Example for testing"))
    monkeypatch.setitem(sys.modules, "ogo.cli.example", SimpleNamespace(main=handler))
    with pytest.raises(type(error)):
        cli.main(["example"])
    assert (logger.handlers, logger.level, logger.propagate) == previous


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
