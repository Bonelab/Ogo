import importlib.util
from pathlib import Path

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
        ],
    ]
    assert solve_calls == [
        (Path("models/density_vertebra_2_L1.n88model"), "spine", calls[0]),
        (Path("models/density_vertebra_4_T12.n88model"), "spine", calls[1]),
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
        (Path("models/density_LT_FEMUR_SF.n88model"), "hip", calls[0]),
        (Path("models/density_RT_FEMUR_SF.n88model"), "hip", calls[1]),
    ]


def test_hip_can_run_one_side(monkeypatch, solve_calls):
    calls = []
    monkeypatch.setattr(GenerateFEM, "run_femur_command", lambda argv: calls.append(list(argv)))

    GenerateFEM.main(["hip", "density.nii.gz", "hip_mask.nii.gz", "--side", "right"])

    assert calls == [["density.nii.gz", "hip_mask.nii.gz", "--femur_side", "2"]]
    assert solve_calls == [(Path("density_RT_FEMUR_SF.n88model"), "hip", calls[0])]


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
        "--mask_threshold 6 --process_mask_threshold 5 --appendix L2"
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

    assert calls[0][-2:] == ["--fe_displacement", "-0.25"]
    assert calls[0].count("--fe_displacement") == 2


def test_legacy_model_type_dispatch_still_works(monkeypatch):
    calls = []
    monkeypatch.setattr(GenerateFEM, "run_spine_command", lambda argv: calls.append(list(argv)))

    GenerateFEM.main(
        [
            "--model_type",
            "vertebra",
            "density.nii.gz",
            "mask.nii.gz",
            "--mask_threshold",
            "2",
            "--process_mask_threshold",
            "1",
        ]
    )

    assert calls == [
        ["density.nii.gz", "mask.nii.gz", "--mask_threshold", "2", "--process_mask_threshold", "1"]
    ]
