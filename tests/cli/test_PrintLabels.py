"""Preserved table content and reusable formatting in the PrintLabels pilot."""

import pytest

from ogo.cli.PrintLabels import main, render_labels


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
    result = run_installed_command("ogo", "print-labels", "--echo", mode)
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stderr == ""
    assert header in result.stdout
    # Check data independently of column padding and incidental blank lines.
    rows = [" ".join(line.split()) for line in result.stdout.splitlines()
            if line.strip() and not line.lstrip().startswith("#")]
    assert known_row in rows
    if last_row is not None:
        assert rows[-1] == last_row
        assert len(rows) == (16 if mode == "colors16" else 128)


@pytest.mark.parametrize("mode", ["labels", "colors16", "colors128"])
def test_render_returns_text_without_terminal_output(mode, capsys):
    text = render_labels(mode)
    assert text.startswith("#")
    assert text.endswith("\n")
    assert capsys.readouterr() == ("", "")


@pytest.mark.parametrize("mode", ["lab", "invalid"])
def test_invalid_direct_selection_raises_value_error(mode):
    with pytest.raises(ValueError, match="Unknown label table"):
        render_labels(mode)


def test_main_accepts_explicit_arguments_and_defaults_to_labels(capsys):
    assert main([]) == 0
    output = capsys.readouterr()
    assert output.out == render_labels()
    assert output.err == ""
