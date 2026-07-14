import ast
import re
from pathlib import Path


def test_setup_py_does_not_bootstrap_setup_dependencies():
    """CI installs build dependencies before installing this package."""
    repo_root = Path(__file__).resolve().parents[1]
    setup_py = repo_root / "setup.py"
    tree = ast.parse(setup_py.read_text(encoding="utf-8"), filename=str(setup_py))

    for node in ast.walk(tree):
        if not isinstance(node, ast.keyword):
            continue
        assert node.arg != "setup_requires", (
            "setup.py should not use setup_requires; that setuptools bootstrap path "
            "is fragile in the Python 3.7 CI environment."
        )


def test_build_system_declares_pbr_dependency():
    """Editable installs need pbr before setup.py metadata is evaluated."""
    repo_root = Path(__file__).resolve().parents[1]
    pyproject_toml = repo_root / "pyproject.toml"

    assert pyproject_toml.exists(), "pyproject.toml should declare build requirements."
    text = pyproject_toml.read_text(encoding="utf-8")
    assert "[build-system]" in text
    assert '"pbr' in text


def test_ci_conda_install_is_non_interactive():
    """Cache misses in CI should not wait for an install confirmation prompt."""
    repo_root = Path(__file__).resolve().parents[1]
    workflow = repo_root / ".github" / "workflows" / "main.yml"
    text = workflow.read_text(encoding="utf-8")

    install_lines = re.findall(r"^\s*(?:conda|mamba) install\b.*$", text, flags=re.MULTILINE)
    assert install_lines, "The CI workflow should install dependencies."
    for line in install_lines:
        assert " -y " in " {0} ".format(line), (
            "CI install commands should include -y so cache misses "
            "cannot block on confirmation."
        )


def test_ci_uses_mamba_for_dependency_resolution():
    """Cold CI environments should use the faster conda-compatible solver."""
    repo_root = Path(__file__).resolve().parents[1]
    workflow = repo_root / ".github" / "workflows" / "main.yml"
    text = workflow.read_text(encoding="utf-8")

    assert "mamba-version:" in text
    assert re.search(r"^\s*mamba install\b", text, flags=re.MULTILINE)
