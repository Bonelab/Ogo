import csv
import json
from pathlib import Path
import shutil

import pytest

from ogo.fea.spine import build_benchmark_sample_model, find_spinefe_benchmark_dir


def _node_set_count(model, name):
    node_set = model.GetNodeSet(name)
    if node_set is None:
        raise AssertionError(f"Model is missing node set {name!r}")
    return node_set.GetNumberOfTuples()


def _model_count(model, method_names):
    for method_name in method_names:
        method = getattr(model, method_name, None)
        if method is not None:
            return method()
    raise AssertionError(f"Model does not expose any of {method_names!r}")


def _sample_paths():
    benchmark_dir = find_spinefe_benchmark_dir(Path(__file__).resolve().parents[2])
    if benchmark_dir is None:
        pytest.skip("spineFE-benchmark checkout not found")

    image = benchmark_dir / "data" / "nii_files" / "sub-verse004_0000_vertebra_21_im.nii.gz"
    mask = benchmark_dir / "data" / "nii_files" / "sub-verse004_0000_vertebra_21_seg.nii.gz"
    if not image.exists() or not mask.exists():
        pytest.skip("spineFE-benchmark sample NIfTI files not found")
    return image, mask


@pytest.mark.integration
@pytest.mark.parametrize("case_name, nonlinear", [("linear", False), ("nonlinear", True)])
def test_spinefe_benchmark_sample_model_generation_matches_golden_metadata(
    tmp_path, case_name, nonlinear
):
    pytest.importorskip("vtk")
    pytest.importorskip("vtkbone")
    pytest.importorskip("SimpleITK")

    image, mask = _sample_paths()
    expected_path = Path(__file__).parent / "fixtures" / "spinefe_sub_verse004_vertebra21_metadata.json"
    expected = json.loads(expected_path.read_text(encoding="utf-8"))[case_name]

    model = build_benchmark_sample_model(
        image,
        mask,
        tmp_path / f"sub-verse004_0000_vertebra_21_{case_name}.n88model",
        nonlinear=nonlinear,
    )

    assert _model_count(model, ["GetNumberOfCells", "GetNumberOfElements"]) == expected["mesh"][
        "elements"
    ]
    assert _model_count(model, ["GetNumberOfPoints", "GetNumberOfNodes"]) == expected["mesh"][
        "nodes"
    ]
    for node_set_name, expected_count in expected["node_sets"].items():
        assert _node_set_count(model, node_set_name) == expected_count


@pytest.mark.integration
def test_spinefe_benchmark_solve_outputs_match_golden_csv(tmp_path):
    pytest.importorskip("vtk")
    pytest.importorskip("vtkbone")
    pytest.importorskip("SimpleITK")

    required_tools = [
        "n88modelinfo",
        "n88derivedfields",
        "n88postfaim",
        "n88pistoia",
        "n88tabulate",
        "n88copymodel",
    ]
    if not (shutil.which("faim") or shutil.which("n88solver_slt")):
        pytest.skip("FAIM/N88 solver command not found")
    missing_tools = [tool for tool in required_tools if shutil.which(tool) is None]
    if missing_tools:
        pytest.skip("Incomplete FAIM/N88 postprocessing toolchain: " + ", ".join(missing_tools))

    golden_csv = Path(__file__).parent / "fixtures" / "spinefe_sub_verse004_vertebra21_solve.csv"
    if not golden_csv.exists():
        pytest.skip("Solve-level golden CSV has not been generated on a full FAIM/N88 install")

    from ogo.util.faim import run_faim_pipeline

    image, mask = _sample_paths()
    model_path = tmp_path / "sub-verse004_0000_vertebra_21_linear.n88model"
    build_benchmark_sample_model(image, mask, model_path, nonlinear=False)
    outputs = run_faim_pipeline(
        model_file=model_path,
        analysis_var="fz_ns1",
        pistoia_vars=["pis_fz_fail", "pis_stiffz"],
        failure_axis="z",
        applied_displacement=-0.2,
        target_displacement=0.68,
        report_profile="spine",
    )

    with golden_csv.open(newline="") as handle:
        expected = next(csv.DictReader(handle))
    with outputs["results_csv"].open(newline="") as handle:
        actual = next(csv.DictReader(handle))

    keys = ["reaction_force_N", "stiffness_N_per_mm"]
    if any(key not in expected for key in keys):
        pytest.skip("Solve-level golden CSV uses the old result schema")
    for key in keys:
        assert float(actual[key]) == pytest.approx(float(expected[key]), rel=1e-4)
