import ast
from pathlib import Path


def test_fe_code_does_not_use_python_310_zip_strict_keyword():
    """The project CI still runs FE tests on Python 3.7."""
    repo_root = Path(__file__).resolve().parents[2]
    for path in repo_root.joinpath("ogo", "fea").glob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            if not isinstance(node.func, ast.Name) or node.func.id != "zip":
                continue
            assert all(keyword.arg != "strict" for keyword in node.keywords), (
                f"{path.relative_to(repo_root)} uses zip(..., strict=...), "
                "which is not available on Python 3.7."
            )
