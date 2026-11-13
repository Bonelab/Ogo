"""Installed-command checks shared by the pilot characterization tests."""

import os
from pathlib import Path
import subprocess
import sys
import sysconfig

import pytest


@pytest.fixture(scope="session")
def installed_ogo_identity(tmp_path_factory):
    """Do not let an inherited installation silently supply baseline results."""
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    env.pop("PYTHONHOME", None)
    result = subprocess.run(
        [sys.executable, "-c", "import ogo; print(ogo.__file__)"],
        cwd=tmp_path_factory.mktemp("installed-identity"),
        env=env,
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        timeout=60,
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
        env = os.environ.copy()
        env.pop("PYTHONPATH", None)
        env.pop("PYTHONHOME", None)
        return subprocess.run(
            [str(executable), *map(str, arguments)],
            cwd=tmp_path,
            env=env,
            stdin=subprocess.DEVNULL,
            capture_output=True,
            text=True,
            timeout=60,
        )

    return run
