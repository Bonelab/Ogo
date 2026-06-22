import csv
import importlib.util
from pathlib import Path


_MODULE_PATH = Path(__file__).resolve().parents[2] / "ogo" / "util" / "faim.py"
_SPEC = importlib.util.spec_from_file_location("faim", str(_MODULE_PATH))
faim = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(faim)


def test_resolve_faim_command_uses_bin_64bit(tmp_path):
    install_root = tmp_path / "Faim 9.0"
    bin_dir = install_root / "bin" / "64bit"
    bin_dir.mkdir(parents=True)
    command = bin_dir / "n88derivedfields"
    command.write_text("", encoding="utf-8")

    assert faim.resolve_faim_bin_dir(install_root=install_root) == bin_dir
    assert faim.resolve_faim_command("n88derivedfields", install_root=install_root) == str(command)


def test_write_results_csv_collects_pistoia_and_loads(tmp_path):
    analysis_file = tmp_path / "model_analysis.txt"
    pistoia_file = tmp_path / "model_pistoia.txt"
    loads_csv = tmp_path / "model_loads.csv"
    pistoia_csv = tmp_path / "model_pistoia.csv"
    results_csv = tmp_path / "model_results.csv"

    analysis_file.write_text("analysis", encoding="utf-8")
    pistoia_file.write_text(
        "\n".join(
            [
                "Critical volume (%): 2.000",
                "Critical EES: 7.0000E-03",
                "EES at vol_crit: 3.5000E-03",
                "Failure load (RF * factor): 1.0000E+02 2.0000E+02 3.0000E+02",
                "Axial stiffness: 1.0000E+03 2.0000E+03 3.0000E+03",
            ]
        ),
        encoding="utf-8",
    )
    loads_csv.write_text("fz_ns1\n123.4\n", encoding="utf-8")
    pistoia_csv.write_text("pis_fz_fail,pis_stiffz\n345.6,789.0\n", encoding="utf-8")

    faim.write_results_csv(
        output_file=results_csv,
        model_file=tmp_path / "model.n88model",
        analysis_file=analysis_file,
        pistoia_file=pistoia_file,
        loads_csv=loads_csv,
        pistoia_csv=pistoia_csv,
        analysis_var="fz_ns1",
        pistoia_vars=["pis_fz_fail", "pis_stiffz"],
        failure_axis="z",
        applied_displacement="-1.0",
        target_displacement=0.2,
    )

    with results_csv.open(newline="") as handle:
        row = next(csv.DictReader(handle))

    assert row["reaction_force_N"] == "123.4"
    assert row["load_at_0p2_percent_N"] == "24.680000000000003"
    assert row["pistoia_failure_load_N"] == "345.6"
    assert row["stiffness_N_per_mm"] == "789.0"
    assert row["critical_volume_pct"] == "2.0"
    assert row["rescale_factor_to_target"] == "0.2"


def test_run_faim_pipeline_dry_run_prints_full_workflow(tmp_path, capsys):
    model = tmp_path / "model.n88model"
    model.write_text("", encoding="utf-8")

    faim.run_faim_pipeline(
        model_file=model,
        analysis_var="fz_ns1",
        pistoia_vars=["pis_fz_fail", "pis_stiffz"],
        failure_axis="z",
        dry_run=True,
        threads=2,
        compress=False,
    )

    out = capsys.readouterr().out
    assert "n88modelinfo --solutions" in out
    assert "faim --engine=mt --threads=2" in out
    assert "n88derivedfields" in out
    assert "n88postfaim --output_file" in out
    assert "n88pistoia" in out
    assert "n88tabulate -H -d , -V fz_ns1" in out
    assert "n88tabulate -H -d , -V pis_fz_fail,pis_stiffz" in out
