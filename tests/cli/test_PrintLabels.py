"""Known table content before migrating the installed PrintLabels command."""

import pytest


@pytest.mark.parametrize("mode,header,known_row,last_row", [
    (
        "labels",
        "# ITK-SnAP Label Description File",
        '1 200 68 68 1.0 1 1 "Femur Right"',
        None,
    ),
    (
        "colors16",
        "# IND    R     G     B   DESC  - 16 colors",
        "1 0 0 0 Black",
        "16 0 0 128 Navy",
    ),
    (
        "colors128",
        "# IND    R     G     B   DESC  - 128 colors",
        "1 128 0 0 maroon",
        "128 255 255 255 white",
    ),
])
def test_known_table_entries(run_installed_command, mode, header, known_row, last_row):
    result = run_installed_command("ogoPrintLabels", "--echo", mode)
    assert result.returncode == 0, result.stdout + result.stderr
    assert header in result.stdout
    # Check data independently of column padding and incidental blank lines.
    rows = [" ".join(line.split()) for line in result.stdout.splitlines()
            if line.strip() and not line.lstrip().startswith("#")]
    assert known_row in rows
    if last_row is not None:
        assert rows[-1] == last_row
        assert len(rows) == (16 if mode == "colors16" else 128)
