import csv
import importlib.util
from pathlib import Path

import pytest


_MODULE_PATH = Path(__file__).resolve().parents[2] / "ogo" / "util" / "faim.py"
_SPEC = importlib.util.spec_from_file_location("faim", str(_MODULE_PATH))
faim = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(faim)


def write_minimal_n88_model(path, coords, node_sets):
    netCDF4 = pytest.importorskip("netCDF4")

    with netCDF4.Dataset(str(path), "w") as root:
        parts = root.createGroup("Parts")
        part = parts.createGroup("Part1")
        part.createDimension("NumberOfNodes", len(coords))
        part.createDimension("NumberOfCoordinates", 3)
        coordinates = part.createVariable(
            "NodeCoordinates",
            "f8",
            ("NumberOfNodes", "NumberOfCoordinates"),
        )
        coordinates[:] = coords

        sets = root.createGroup("Sets")
        node_sets_group = sets.createGroup("NodeSets")
        for name, node_numbers in node_sets.items():
            node_set = node_sets_group.createGroup(name)
            node_set.createDimension("NumberOfNodes", len(node_numbers))
            variable = node_set.createVariable("NodeNumber", "i4", ("NumberOfNodes",))
            variable[:] = node_numbers


def write_minimal_n88_model_with_constraints(path, coords, constraints):
    netCDF4 = pytest.importorskip("netCDF4")

    with netCDF4.Dataset(str(path), "w") as root:
        parts = root.createGroup("Parts")
        part = parts.createGroup("Part1")
        part.createDimension("NumberOfNodes", len(coords))
        part.createDimension("NumberOfCoordinates", 3)
        coordinates = part.createVariable(
            "NodeCoordinates",
            "f8",
            ("NumberOfNodes", "NumberOfCoordinates"),
        )
        coordinates[:] = coords

        constraints_group = root.createGroup("Constraints")
        for name, node_numbers in constraints.items():
            constraint = constraints_group.createGroup(name)
            constraint.createDimension("NumberOfValues", len(node_numbers))
            node_number = constraint.createVariable("NodeNumber", "i4", ("NumberOfValues",))
            node_number[:] = node_numbers
            sense = constraint.createVariable("Sense", "i4", ("NumberOfValues",))
            sense[:] = [3] * len(node_numbers)
            value = constraint.createVariable("Value", "f8", ("NumberOfValues",))
            value[:] = [0.0] * len(node_numbers)


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
    assert row["load_at_target_displacement_N"] == "24.680000000000003"
    assert "load_at_0p68_percent_N" not in row
    assert row["pistoia_failure_load_N"] == "345.6"
    assert row["stiffness_N_per_mm"] == "123.4"
    assert row["critical_volume_pct"] == "2.0"
    assert row["rescale_factor_to_target"] == "0.2"


def test_write_results_csv_reports_generic_target_without_0p68_alias(tmp_path):
    analysis_file = tmp_path / "model_analysis.txt"
    pistoia_file = tmp_path / "model_pistoia.txt"
    loads_csv = tmp_path / "model_loads.csv"
    pistoia_csv = tmp_path / "model_pistoia.csv"
    results_csv = tmp_path / "model_results.csv"

    analysis_file.write_text("analysis", encoding="utf-8")
    pistoia_file.write_text("", encoding="utf-8")
    loads_csv.write_text("fy_ns1\n100.0\n", encoding="utf-8")
    pistoia_csv.write_text("pis_fy_fail,pis_stiffy\n50.0,25.0\n", encoding="utf-8")

    faim.write_results_csv(
        output_file=results_csv,
        model_file=tmp_path / "model.n88model",
        analysis_file=analysis_file,
        pistoia_file=pistoia_file,
        loads_csv=loads_csv,
        pistoia_csv=pistoia_csv,
        analysis_var="fy_ns1",
        pistoia_vars=["pis_fy_fail", "pis_stiffy"],
        failure_axis="y",
        applied_displacement="1.0",
        target_displacement=4.0,
    )

    with results_csv.open(newline="") as handle:
        row = next(csv.DictReader(handle))

    assert row["load_at_target_displacement_N"] == "400.0"
    assert row["load_at_target_displacement_kN"] == "0.4"
    assert "load_at_0p68_percent_N" not in row


def test_write_results_csv_femur_profile_reports_4_percent_and_stiffness(tmp_path):
    analysis_file = tmp_path / "model_analysis.txt"
    model_file = tmp_path / "model.n88model"
    pistoia_file = tmp_path / "model_pistoia.txt"
    loads_csv = tmp_path / "model_loads.csv"
    pistoia_csv = tmp_path / "model_pistoia.csv"
    results_csv = tmp_path / "model_results.csv"

    analysis_file.write_text("analysis", encoding="utf-8")
    pistoia_file.write_text("", encoding="utf-8")
    loads_csv.write_text("fy_ns1\n100.0\n", encoding="utf-8")
    pistoia_csv.write_text("", encoding="utf-8")
    write_minimal_n88_model(
        model_file,
        coords=[
            [0.0, 0.0, 0.0],
            [0.0, 100.0, 0.0],
        ],
        node_sets={
            "Femoral_Head_PMMA_Nodes": [1],
            "Greater_Trochanter_PMMA_Nodes": [2],
        },
    )

    faim.write_results_csv(
        output_file=results_csv,
        model_file=model_file,
        analysis_file=analysis_file,
        pistoia_file=pistoia_file,
        loads_csv=loads_csv,
        pistoia_csv=pistoia_csv,
        analysis_var="fy_ns1",
        pistoia_vars=[],
        failure_axis="y",
        applied_displacement="1.0",
        target_displacement=4.0,
        report_profile="femur",
    )

    with results_csv.open(newline="") as handle:
        row = next(csv.DictReader(handle))

    assert list(row) == [
        "model_file",
        "analysis_file",
        "analysis_var",
        "applied_displacement",
        "reaction_force_N",
        "stiffness_N_per_mm",
        "characteristic_length_mm",
    ]
    assert row["analysis_var"] == "fy_ns1"
    assert row["applied_displacement"] == "1.0"
    assert row["reaction_force_N"] == "100.0"
    assert row["characteristic_length_mm"] == "100.0"
    assert row["stiffness_N_per_mm"] == "100.0"


def test_write_results_csv_spine_profile_reports_0p68_percent(tmp_path):
    analysis_file = tmp_path / "model_analysis.txt"
    model_file = tmp_path / "model.n88model"
    pistoia_file = tmp_path / "model_pistoia.txt"
    loads_csv = tmp_path / "model_loads.csv"
    pistoia_csv = tmp_path / "model_pistoia.csv"
    results_csv = tmp_path / "model_results.csv"

    analysis_file.write_text("analysis", encoding="utf-8")
    pistoia_file.write_text("", encoding="utf-8")
    loads_csv.write_text("fz_ns1\n200.0\n", encoding="utf-8")
    pistoia_csv.write_text("", encoding="utf-8")
    write_minimal_n88_model(
        model_file,
        coords=[
            [0.0, 0.0, 100.0],
            [0.0, 0.0, 0.0],
        ],
        node_sets={
            "body_top": [1],
            "body_bottom": [2],
        },
    )

    faim.write_results_csv(
        output_file=results_csv,
        model_file=model_file,
        analysis_file=analysis_file,
        pistoia_file=pistoia_file,
        loads_csv=loads_csv,
        pistoia_csv=pistoia_csv,
        analysis_var="fz_ns1",
        pistoia_vars=[],
        failure_axis="z",
        applied_displacement="-2.0",
        target_displacement=0.68,
        report_profile="spine",
    )

    with results_csv.open(newline="") as handle:
        row = next(csv.DictReader(handle))

    assert list(row) == [
        "model_file",
        "analysis_file",
        "analysis_var",
        "applied_displacement",
        "reaction_force_N",
        "stiffness_N_per_mm",
        "characteristic_length_mm",
    ]
    assert row["analysis_var"] == "fz_ns1"
    assert row["applied_displacement"] == "-2.0"
    assert row["reaction_force_N"] == "200.0"
    assert row["characteristic_length_mm"] == "100.0"
    assert row["stiffness_N_per_mm"] == "100.0"


def test_write_results_csv_spine_profile_converts_percent_from_geometry(tmp_path):
    analysis_file = tmp_path / "model_analysis.txt"
    model_file = tmp_path / "model.n88model"
    pistoia_file = tmp_path / "model_pistoia.txt"
    loads_csv = tmp_path / "model_loads.csv"
    pistoia_csv = tmp_path / "model_pistoia.csv"
    results_csv = tmp_path / "model_results.csv"

    analysis_file.write_text("analysis", encoding="utf-8")
    pistoia_file.write_text("", encoding="utf-8")
    loads_csv.write_text("fz_ns1\n200.0\n", encoding="utf-8")
    pistoia_csv.write_text("", encoding="utf-8")
    write_minimal_n88_model(
        model_file,
        coords=[
            [0.0, 0.0, 50.0],
            [0.0, 0.0, 0.0],
        ],
        node_sets={
            "body_top": [1],
            "body_bottom": [2],
        },
    )

    faim.write_results_csv(
        output_file=results_csv,
        model_file=model_file,
        analysis_file=analysis_file,
        pistoia_file=pistoia_file,
        loads_csv=loads_csv,
        pistoia_csv=pistoia_csv,
        analysis_var="fz_ns1",
        pistoia_vars=[],
        failure_axis="z",
        applied_displacement="-2.0",
        target_displacement=0.68,
        report_profile="spine",
    )

    with results_csv.open(newline="") as handle:
        row = next(csv.DictReader(handle))

    assert row["characteristic_length_mm"] == "50.0"
    assert row["reaction_force_N"] == "200.0"
    assert row["stiffness_N_per_mm"] == "100.0"


def test_write_results_csv_spine_profile_excludes_pistoia_fields(tmp_path):
    analysis_file = tmp_path / "model_analysis.txt"
    model_file = tmp_path / "model.n88model"
    pistoia_file = tmp_path / "model_pistoia.txt"
    loads_csv = tmp_path / "model_loads.csv"
    pistoia_csv = tmp_path / "model_pistoia.csv"
    results_csv = tmp_path / "model_results.csv"

    analysis_file.write_text("analysis", encoding="utf-8")
    pistoia_file.write_text("Critical volume (%): 2.000\n", encoding="utf-8")
    loads_csv.write_text("fz_ns1\n200.0\n", encoding="utf-8")
    pistoia_csv.write_text("pis_fz_fail,pis_stiffz\n-300.0,120.0\n", encoding="utf-8")
    write_minimal_n88_model(
        model_file,
        coords=[
            [0.0, 0.0, 100.0],
            [0.0, 0.0, 0.0],
        ],
        node_sets={
            "body_top": [1],
            "body_bottom": [2],
        },
    )

    faim.write_results_csv(
        output_file=results_csv,
        model_file=model_file,
        analysis_file=analysis_file,
        pistoia_file=pistoia_file,
        loads_csv=loads_csv,
        pistoia_csv=pistoia_csv,
        analysis_var="fz_ns1",
        pistoia_vars=["pis_fz_fail", "pis_stiffz"],
        failure_axis="z",
        applied_displacement="-2.0",
        target_displacement=0.68,
        report_profile="spine",
    )

    with results_csv.open(newline="") as handle:
        row = next(csv.DictReader(handle))

    assert "pistoia_failure_load_N" not in row
    assert "pistoia_failure_load_kN" not in row
    assert "pistoia_stiffness_N_per_mm" not in row


def test_write_results_csv_spine_profile_can_use_constraint_groups(tmp_path):
    analysis_file = tmp_path / "model_analysis.txt"
    model_file = tmp_path / "model.n88model"
    pistoia_file = tmp_path / "model_pistoia.txt"
    loads_csv = tmp_path / "model_loads.csv"
    pistoia_csv = tmp_path / "model_pistoia.csv"
    results_csv = tmp_path / "model_results.csv"

    analysis_file.write_text("analysis", encoding="utf-8")
    pistoia_file.write_text("", encoding="utf-8")
    loads_csv.write_text("fz_ns1\n200.0\n", encoding="utf-8")
    pistoia_csv.write_text("", encoding="utf-8")
    write_minimal_n88_model_with_constraints(
        model_file,
        coords=[
            [0.0, 0.0, 100.0],
            [0.0, 0.0, 0.0],
        ],
        constraints={
            "top_displacement": [1],
            "bottom_fixed_z": [2],
        },
    )

    faim.write_results_csv(
        output_file=results_csv,
        model_file=model_file,
        analysis_file=analysis_file,
        pistoia_file=pistoia_file,
        loads_csv=loads_csv,
        pistoia_csv=pistoia_csv,
        analysis_var="fz_ns1",
        pistoia_vars=[],
        failure_axis="z",
        applied_displacement="-2.0",
        target_displacement=0.68,
        report_profile="spine",
    )

    with results_csv.open(newline="") as handle:
        row = next(csv.DictReader(handle))

    assert row["characteristic_length_mm"] == "100.0"
    assert row["reaction_force_N"] == "200.0"


def test_set_prescribed_displacement_from_percent_updates_constraints(tmp_path):
    model_file = tmp_path / "model.n88model"
    write_minimal_n88_model_with_constraints(
        model_file,
        coords=[
            [0.0, 0.0, 100.0],
            [0.0, 0.0, 0.0],
        ],
        constraints={
            "top_displacement": [1],
            "bottom_fixed_z": [2],
            "convergence_set": [1],
        },
    )

    update = faim.set_prescribed_displacement_from_percent(
        model_file,
        report_profile="spine",
        failure_axis="z",
        target_displacement_percent=0.68,
        displacement_sign=-1.0,
    )

    assert update["characteristic_length_mm"] == 100.0
    assert update["target_displacement_mm"] == 0.68
    assert update["applied_displacement_mm"] == -0.68
    assert faim.read_prescribed_displacement(model_file, "spine") == -0.68


def test_set_prescribed_displacement_from_percent_updates_femur_constraints(tmp_path):
    model_file = tmp_path / "model.n88model"
    write_minimal_n88_model_with_constraints(
        model_file,
        coords=[
            [0.0, 100.0, 0.0],
            [0.0, 0.0, 0.0],
        ],
        constraints={
            "top_displacement": [1],
            "bottom_fixed_y_PMMA": [2],
            "convergence_set": [1],
        },
    )

    update = faim.set_prescribed_displacement_from_percent(
        model_file,
        report_profile="femur",
        failure_axis="y",
        target_displacement_percent=4.0,
        displacement_sign=1.0,
    )

    assert update["characteristic_length_mm"] == 100.0
    assert update["target_displacement_mm"] == 4.0
    assert update["applied_displacement_mm"] == 4.0
    assert faim.read_prescribed_displacement(model_file, "femur") == 4.0


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


def test_run_faim_pipeline_can_skip_pistoia(tmp_path, capsys):
    model = tmp_path / "model.n88model"
    model.write_text("", encoding="utf-8")

    faim.run_faim_pipeline(
        model_file=model,
        analysis_var="fy_ns1",
        pistoia_vars=["pis_fy_fail", "pis_stiffy"],
        failure_axis="y",
        run_pistoia=False,
        dry_run=True,
        compress=False,
    )

    out = capsys.readouterr().out
    assert "n88pistoia" not in out
    assert "pis_fy_fail,pis_stiffy" not in out
    assert "n88tabulate -H -d , -V fy_ns1" in out
