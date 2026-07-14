import ast
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
