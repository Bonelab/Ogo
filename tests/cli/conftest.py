"""Outside-checkout checks for installed commands and unmigrated modules."""

import os
from pathlib import Path
import subprocess
import sys
import sysconfig

import pytest


def _run(arguments, directory):
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    env.pop("PYTHONHOME", None)
    return subprocess.run(
        list(map(str, arguments)),
        cwd=directory,
        env=env,
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        timeout=60,
    )


@pytest.fixture(scope="session")
def installed_ogo_identity(tmp_path_factory):
    """Do not let an inherited installation silently supply baseline results."""
    result = _run(
        [sys.executable, "-c", "import ogo; print(ogo.__file__)"],
        tmp_path_factory.mktemp("installed-identity"),
    )
    assert result.returncode == 0, result.stdout + result.stderr
    expected = Path(__file__).resolve().parents[2] / "ogo" / "__init__.py"
    assert Path(result.stdout.strip()).resolve() == expected, (
        f"Installed commands must use {expected}, got {result.stdout.strip()}"
    )


@pytest.fixture
def run_installed_command(tmp_path, installed_ogo_identity):
    """Run this environment's executable outside the checkout, without stdin."""
    def run(name, *arguments):
        suffix = ".exe" if os.name == "nt" else ""
        executable = Path(sysconfig.get_path("scripts")) / (name + suffix)
        return _run([executable, *arguments], tmp_path)

    return run


@pytest.fixture
def run_python(tmp_path, installed_ogo_identity):
    """Run Python for import checks and tests of source awaiting migration."""
    def run(*arguments):
        return _run([sys.executable, *arguments], tmp_path)

    return run
