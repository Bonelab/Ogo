import ast
from importlib import metadata
import os
from pathlib import Path
import subprocess
import sysconfig

from packaging.requirements import Requirement
from packaging.version import Version


def test_setup_py_does_not_bootstrap_setup_dependencies():
    """Build isolation supplies dependencies without setuptools bootstrapping."""
    repo_root = Path(__file__).resolve().parents[1]
    setup_py = repo_root / "setup.py"
    tree = ast.parse(setup_py.read_text(encoding="utf-8"), filename=str(setup_py))

    for node in ast.walk(tree):
        if not isinstance(node, ast.keyword):
            continue
        assert node.arg != "setup_requires", (
            "setup.py should not use setup_requires; that setuptools bootstrap path "
            "is fragile in isolated builds."
        )


def test_build_system_declares_pbr_dependency():
    """Editable installs need pbr before setup.py metadata is evaluated."""
    repo_root = Path(__file__).resolve().parents[1]
    pyproject_toml = repo_root / "pyproject.toml"

    assert pyproject_toml.exists(), "pyproject.toml should declare build requirements."
    text = pyproject_toml.read_text(encoding="utf-8")
    assert "[build-system]" in text
    assert '"pbr' in text


def test_qc_image_dependencies_are_declared():
    requirements = Path(__file__).resolve().parents[1] / "requirements.txt"
    names = {
        Requirement(line).name.lower()
        for line in requirements.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    }
    assert {"nibabel", "pillow", "netcdf4"} <= names


def test_editable_install_uses_checkout_and_matching_version():
    """CI must exercise this checkout rather than another installed Ogo."""
    import ogo

    root = Path(__file__).resolve().parents[1]
    assert Path(ogo.__file__).resolve() == root / "ogo" / "__init__.py"
    assert Version(ogo.__version__) == Version(metadata.version("ogo"))


def test_installed_label_command_runs_outside_checkout(tmp_path):
    """Exercise a real installed command, not just its help or parser."""
    executable = "ogo.exe" if os.name == "nt" else "ogo"
    command = Path(sysconfig.get_path("scripts")) / executable
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    env.pop("PYTHONHOME", None)
    result = subprocess.run(
        [str(command), "print-labels", "--echo", "labels"], cwd=tmp_path, env=env,
        capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert '"Femur Right"' in result.stdout
    assert result.stderr == ""


def test_installed_package_registers_only_the_root_command():
    entry_points = {
        entry.name: entry.value
        for entry in metadata.distribution("ogo").entry_points
        if entry.group == "console_scripts"
    }
    assert entry_points == {"ogo": "ogo.cli.cli:main"}
