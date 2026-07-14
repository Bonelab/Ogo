import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest


_MODULE_PATH = Path(__file__).resolve().parents[2] / "ogo" / "cli" / "GenerateFEM.py"
_SPEC = importlib.util.spec_from_file_location("GenerateFEM", str(_MODULE_PATH))
GenerateFEM = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(GenerateFEM)


BENCHMARK_LINEAR_ARGS = [
    "--fe_displacement",
    "-0.2",
    "--elastic_E_func",
    "kopperdahl_trab_E",
    "--cort_elastic_E_func",
    "kopperdahl_trab_E",
]

BENCHMARK_NONLINEAR_ARGS = [
    "--fe_displacement",
    "-2.0",
    "--pmma_yield_compression",
    "70.0",
    "--pmma_yield_tension",
    "70.0",
    "--elastic_E_func",
    "kopperdahl_trab_E",
    "--yield_comp_func",
    "kopperdahl_trab_yc",
    "--yield_tens_func",
    "kopperdahl_trab_yc",
    "--cort_elastic_E_func",
    "kopperdahl_trab_E",
    "--cort_yield_comp_func",
    "kopperdahl_trab_yc",
    "--cort_yield_tens_func",
    "kopperdahl_trab_yc",
]

SPINE_DEFAULT_QC_ARGS = ["--quality_control", "False"]


@pytest.fixture
def solve_calls(monkeypatch):
    calls = []
    monkeypatch.setattr(
        GenerateFEM,
        "solve_model",
        lambda model_path, args, model_type, generator_argv: calls.append(
            (Path(model_path), model_type, list(generator_argv))
        ),
    )
    return calls


def test_spine_runs_each_requested_vertebra(monkeypatch, solve_calls):
    calls = []
    monkeypatch.setattr(GenerateFEM, "run_spine_command", lambda argv: calls.append(list(argv)))

    GenerateFEM.main(
        [
            "spine",
            "density.nii.gz",
            "spine_mask.nii.gz",
            "--output_path",
            "models",
            "--vertebra",
            "L1:2:1",
            "--vertebra",
            "T12:4:3",
            "--iso_resolution",
            "1.2",
        ]
    )

    assert calls == [
        [
            "density.nii.gz",
            "spine_mask.nii.gz",
            "--mask_threshold",
            "2",
            "--process_mask_threshold",
            "1",
            "--appendix",
            "L1",
            "--output_path",
            "models",
            *BENCHMARK_LINEAR_ARGS,
            "--iso_resolution",
            "1.2",
            *SPINE_DEFAULT_QC_ARGS,
        ],
        [
            "density.nii.gz",
            "spine_mask.nii.gz",
            "--mask_threshold",
            "4",
            "--process_mask_threshold",
            "3",
            "--appendix",
            "T12",
            "--output_path",
            "models",
            *BENCHMARK_LINEAR_ARGS,
            "--iso_resolution",
            "1.2",
            *SPINE_DEFAULT_QC_ARGS,
        ],
    ]
    assert solve_calls == [
        (Path("models/density_L1.n88model"), "spine", calls[0]),
        (Path("models/density_T12.n88model"), "spine", calls[1]),
    ]


def test_hip_defaults_to_both_sides(monkeypatch, solve_calls):
    calls = []
    monkeypatch.setattr(GenerateFEM, "run_femur_command", lambda argv: calls.append(list(argv)))

    GenerateFEM.main(["hip", "density.nii.gz", "hip_mask.nii.gz", "--output_path", "models"])

    assert calls == [
        ["density.nii.gz", "hip_mask.nii.gz", "--femur_side", "1", "--output_path", "models"],
        ["density.nii.gz", "hip_mask.nii.gz", "--femur_side", "2", "--output_path", "models"],
    ]
    assert solve_calls == [
        (Path("models/density_LF.n88model"), "hip", calls[0]),
        (Path("models/density_RF.n88model"), "hip", calls[1]),
    ]


def test_hip_can_run_one_side(monkeypatch, solve_calls):
    calls = []
    monkeypatch.setattr(GenerateFEM, "run_femur_command", lambda argv: calls.append(list(argv)))

    GenerateFEM.main(["hip", "density.nii.gz", "hip_mask.nii.gz", "--side", "right"])

    assert calls == [["density.nii.gz", "hip_mask.nii.gz", "--femur_side", "2"]]
    assert solve_calls == [(Path("density_RF.n88model"), "hip", calls[0])]


def test_dry_run_prints_commands_without_running(monkeypatch, capsys):
    calls = []
    monkeypatch.setattr(GenerateFEM, "run_spine_command", lambda argv: calls.append(list(argv)))

    GenerateFEM.main(
        [
            "spine",
            "density.nii.gz",
            "spine_mask.nii.gz",
            "--vertebra",
            "L2,6,5",
            "--dry-run",
            "--no-solve",
            "--preset",
            "none",
        ]
    )

    assert calls == []
    assert (
        capsys.readouterr().out.strip()
        == "OgoSpineCompressionFe density.nii.gz spine_mask.nii.gz "
        "--mask_threshold 6 --process_mask_threshold 5 --appendix L2 "
        "--quality_control False"
    )


def test_spine_preset_none_keeps_command_minimal(monkeypatch, solve_calls):
    calls = []
    monkeypatch.setattr(GenerateFEM, "run_spine_command", lambda argv: calls.append(list(argv)))

    GenerateFEM.main(
        [
            "spine",
            "density.nii.gz",
            "spine_mask.nii.gz",
            "--vertebra",
            "L3:8:7",
            "--preset",
            "none",
        ]
    )

    assert calls == [
        [
            "density.nii.gz",
            "spine_mask.nii.gz",
            "--mask_threshold",
            "8",
            "--process_mask_threshold",
            "7",
            "--appendix",
            "L3",
            *SPINE_DEFAULT_QC_ARGS,
        ]
    ]


def test_spine_nonlinear_benchmark_preset(monkeypatch, solve_calls):
    calls = []
    monkeypatch.setattr(GenerateFEM, "run_spine_command", lambda argv: calls.append(list(argv)))

    GenerateFEM.main(
        [
            "spine",
            "density.nii.gz",
            "spine_mask.nii.gz",
            "--vertebra",
            "T10:3:2",
            "--preset",
            "benchmark-nonlinear",
        ]
    )

    assert calls == [
        [
            "density.nii.gz",
            "spine_mask.nii.gz",
            "--mask_threshold",
            "3",
            "--process_mask_threshold",
            "2",
            "--appendix",
            "T10",
            *BENCHMARK_NONLINEAR_ARGS,
            *SPINE_DEFAULT_QC_ARGS,
        ]
    ]


def test_spine_explicit_options_follow_preset_for_overrides(monkeypatch, solve_calls):
    calls = []
    monkeypatch.setattr(GenerateFEM, "run_spine_command", lambda argv: calls.append(list(argv)))

    GenerateFEM.main(
        [
            "spine",
            "density.nii.gz",
            "spine_mask.nii.gz",
            "--vertebra",
            "L4:2:1",
            "--fe_displacement",
            "-0.25",
        ]
    )

    explicit_index = len(calls[0]) - len(SPINE_DEFAULT_QC_ARGS) - 2
    assert calls[0][explicit_index:explicit_index + 2] == ["--fe_displacement", "-0.25"]
    assert calls[0].count("--fe_displacement") == 2


def _metadata_args(tmp_path, no_solve=True):
    return type(
        "Args",
        (),
        {
            "dry_run": False,
            "calibrated_image": tmp_path / "density.nii.gz",
            "bone_mask": tmp_path / "mask.nii.gz",
            "output_path": tmp_path,
            "skip_bc_audit": False,
            "bc_audit_flat_tolerance": 1.0e-4,
            "debug": False,
            "no_solve": no_solve,
            "target_displacement": None,
            "exclude": 5000,
            "critical_volume": 2.0,
            "critical_strain": 0.007,
            "no_compress": False,
        },
    )()


def _solve_args(**overrides):
    data = {
        "threads": 4,
        "faim_env": None,
        "conda_executable": "conda",
        "faim_install_root": None,
        "faim_bin_dir": None,
        "faim_license_dir": None,
        "faim_command": None,
        "n88modelinfo_command": None,
        "n88derivedfields_command": None,
        "n88postfaim_command": None,
        "n88pistoia_command": None,
        "n88tabulate_command": None,
        "n88copymodel_command": None,
        "critical_volume": 2.0,
        "critical_strain": 0.007,
        "exclude": 5000,
        "no_compress": False,
        "require_pistoia": False,
        "dry_run": False,
        "target_displacement": None,
    }
    data.update(overrides)
    return SimpleNamespace(**data)


def test_solve_model_sets_femur_displacement_from_percent(monkeypatch, tmp_path):
    import ogo.util.faim as faim_module

    calls = []
    monkeypatch.setattr(faim_module, "run_faim_pipeline", lambda **kwargs: calls.append(kwargs))

    GenerateFEM.solve_model(
        tmp_path / "density_LT_FEMUR_SF.n88model",
        _solve_args(),
        "hip",
        ["density.nii.gz", "mask.nii.gz", "--fe_displacement", "-4.0"],
    )

    assert calls[0]["report_profile"] == "femur"
    assert calls[0]["target_displacement"] == 4.0
    assert calls[0]["solve_displacement_percent"] == 4.0


def test_spine_modeling_metadata_records_materials_and_bcs(tmp_path):
    model = tmp_path / "density_vertebra_2_L1.n88model"
    model.write_text("dummy")
    argv = [
        "density.nii.gz",
        "mask.nii.gz",
        "--mask_threshold",
        "2",
        "--process_mask_threshold",
        "3",
        "--appendix",
        "L1",
        "--elastic_E_func",
        "kopperdahl_trab_E",
        "--cort_elastic_E_func",
        "kopperdahl_trab_E",
        "--fe_displacement",
        "-0.2",
        "--pmma_intrusion",
        "3",
    ]

    path = GenerateFEM.write_modeling_metadata(model, "spine", argv, _metadata_args(tmp_path))
    data = json.loads(path.read_text())

    assert data["target"] == {"vertebra": "L1", "body_label": 2, "process_label": 3}
    assert data["alignment"]["registration_backend"] == "vtk"
    assert data["geometry"]["model_coordinates"] == "preprocessed_image_physical_space"
    assert data["materials"]["trabecular"]["elastic_E_func"] == "kopperdahl_trab_E"
    assert data["materials"]["cortical"]["material_id_range"] == [129, 256]
    assert data["materials"]["pmma"]["material_id"] == 5000
    assert data["boundary_conditions"]["fixture_geometry"]["superior_cap"]["pmma_intrusion_mm"] == 3.0
    assert data["boundary_conditions"]["fixture_geometry"]["inferior_cap"]["pmma_intrusion_mm"] == 3.0
    assert data["boundary_conditions"]["constraints"][0]["node_set"] == "body_top"
    assert "initial_generator_value_mm" not in data["boundary_conditions"]["constraints"][0]
    assert data["boundary_conditions"]["constraints"][0]["target_displacement_percent"] == 0.68
    assert data["boundary_conditions"]["constraints"][0]["value_source"].startswith("target_displacement_percent")
    assert data["boundary_conditions"]["constraints"][1]["node_set"] == "body_bottom"
    assert data["boundary_conditions"]["constraints"][1]["axes"] == ["x", "y", "z"]
    assert data["boundary_conditions"]["constraints"][1]["value_mm"] == 0.0
    assert data["solve_and_reporting"]["target_displacement_percent"] == 0.68
    assert data["solve_and_reporting"]["run_pistoia"] is False


def test_femur_modeling_metadata_records_materials_shaft_and_bcs(tmp_path):
    model = tmp_path / "density_LT_FEMUR_SF.n88model"
    model.write_text("dummy")
    argv = [
        "density.nii.gz",
        "mask.nii.gz",
        "--femur_side",
        "1",
        "--pmma_thick",
        "6",
        "--pmma_intrusion",
        "6",
        "--femur_lesser_trochanter_distal_offset",
        "50",
        "--compartment_mask",
        "trab_cort.nii.gz",
    ]

    path = GenerateFEM.write_modeling_metadata(model, "hip", argv, _metadata_args(tmp_path))
    data = json.loads(path.read_text())

    assert data["target"]["side"] == "left"
    assert data["geometry"]["model_coordinates"] == "preprocessed_image_physical_space"
    assert data["segmentation"]["compartment_labels"] == {"cortical": 1, "trabecular": 2}
    assert data["shaft_standardization"]["lesser_trochanter_distal_offset_mm"] == 50.0
    assert data["materials"]["include_cortical_region"] is True
    assert data["materials"]["cortical"]["material_id_range"] == [129, 256]
    assert data["boundary_conditions"]["fixture_geometry"]["femoral_head"]["relative_to"] == "model_bbox"
    assert data["boundary_conditions"]["fixture_geometry"]["femoral_head"]["shape"] == "rectangle fixture cap"
    assert data["boundary_conditions"]["fixture_geometry"]["greater_trochanter"]["relative_to"] == "model_bbox"
    assert data["boundary_conditions"]["fixture_geometry"]["greater_trochanter"]["shape"] == "rectangle fixture cap"
    assert data["boundary_conditions"]["fixture_geometry"]["femoral_head"]["pmma_intrusion_mm"] == 6.0
    assert data["boundary_conditions"]["constraints"][0]["node_set"] == "Femoral_Head_PMMA_Nodes"
    assert "initial_generator_value_mm" not in data["boundary_conditions"]["constraints"][0]
    assert data["boundary_conditions"]["constraints"][0]["target_displacement_percent"] == 4.0
    assert data["boundary_conditions"]["constraints"][0]["value_source"].startswith("target_displacement_percent")
    assert data["solve_and_reporting"]["target_displacement_percent"] == 4.0
