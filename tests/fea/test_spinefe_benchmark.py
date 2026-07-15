import csv
import ast
import inspect
import json
from pathlib import Path
import shutil
import textwrap

import pytest

from ogo.fea.model import filter_node_set_to_dominant_coordinate_plane
from ogo.fea.spine import build_benchmark_sample_model, find_spinefe_benchmark_dir


def test_filter_node_set_to_dominant_coordinate_plane_removes_small_rim():
    vtk = pytest.importorskip("vtk")

    class FakeModel:
        def __init__(self):
            self.points = [
                (0.0, 0.0, 2.0),
                (1.0, 0.0, 2.0),
                (0.0, 1.0, 2.0),
                (1.0, 1.0, 2.0),
                (0.0, 0.0, 1.0),
            ]
            self.node_sets = {"cap": vtk.vtkIdTypeArray()}
            self.node_sets["cap"].SetName("cap")
            for point_id in range(len(self.points)):
                self.node_sets["cap"].InsertNextValue(point_id)

        def GetNodeSet(self, name):
            return self.node_sets.get(name)

        def AddNodeSet(self, node_ids):
            self.node_sets[node_ids.GetName()] = node_ids

        def GetPoint(self, point_id):
            return self.points[point_id]

    model = FakeModel()
    filter_node_set_to_dominant_coordinate_plane(model, "cap", axis="z")

    filtered = model.GetNodeSet("cap")
    assert filtered.GetNumberOfTuples() == 4
    assert [filtered.GetValue(i) for i in range(filtered.GetNumberOfTuples())] == [0, 1, 2, 3]


def test_spine_qc_filename_parser_handles_compact_generated_names():
    from ogo.fea.spine import parse_filename

    metadata = parse_filename("/tmp/density_L4_BCcheck.csv")

    assert metadata == {
        "ID": "density",
        "TREATMENT": "L4",
        "LOCATION": "BCcheck",
        "NUMBER": "",
        "filename": "density_L4_BCcheck",
    }


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


def test_merge_vtk_images_can_preserve_existing_material_voxels():
    np = pytest.importorskip("numpy")
    vtk = pytest.importorskip("vtk")
    from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

    from ogo.fea.spine import merge_vtk_images

    material = np.zeros((3, 3, 3), dtype=np.float32)
    material[1, 1, 1] = 42
    disk = np.zeros_like(material)
    disk[1, 1, 1] = 1
    disk[2, 1, 1] = 1

    def vtk_image(array):
        image = vtk.vtkImageData()
        image.SetDimensions(array.shape)
        image.SetOrigin(0, 0, 0)
        image.SetSpacing(1, 1, 1)
        image.GetPointData().SetScalars(numpy_to_vtk(array.ravel(order="F"), deep=True))
        return image

    merged = merge_vtk_images(
        [vtk_image(material), vtk_image(disk)],
        [None, 5000],
        overwrite_existing=False,
    )

    out = vtk_to_numpy(merged.GetPointData().GetScalars()).reshape(material.shape, order="F")
    assert out[1, 1, 1] == 42
    assert out[2, 1, 1] == 5000


def test_spine_workflow_uses_projected_contact_disks():
    pytest.importorskip("vtk")

    from ogo.fea import spine

    source = textwrap.dedent(inspect.getsource(spine.process_vertebra))
    tree = ast.parse(source)
    called_names = {
        node.func.id
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }

    assert "generate_projected_material_disk_vtk" in called_names
    assert "generate_bone_cap_vtk" not in called_names


def test_spine_registration_scaling_uses_voxel_surface_points():
    pytest.importorskip("vtk")

    from ogo.fea import spine

    source = textwrap.dedent(inspect.getsource(spine.get_icp_with_scaling))

    assert "surface_points_from_vtk_mask" in source
    assert "point_cloud_axis_lengths" in source
    assert "perform_marching_cubes(body)" not in source


def test_spine_registration_uses_preprocessed_isotropic_body():
    pytest.importorskip("vtk")

    from ogo.fea import spine

    source = textwrap.dedent(inspect.getsource(spine.process_vertebra))

    assert "SPINE_PREPROCESSING_CROP_MARGIN_MM" in source
    assert "resample_vtk_image_to_spacing" in source
    assert "registration_body = threshold(" in source
    assert "get_icp_with_scaling(\n        registration_body," in source
    assert source.index("resample_vtk_image_to_spacing") < source.index(
        "registration_body = threshold("
    )
    assert source.index("registration_body = threshold(") < source.index(
        "get_icp_with_scaling("
    )


def test_spine_workflow_resamples_to_explicit_reference_grid():
    pytest.importorskip("vtk")

    from ogo.fea import spine

    source = textwrap.dedent(inspect.getsource(spine.process_vertebra))

    assert "output_grid_for_point_transform" in source
    assert "resample_vtk_image_with_point_transform" in source
    assert "transform_resample(" not in source


def test_spine_contact_planes_use_generated_model_bbox():
    pytest.importorskip("vtk")

    from ogo.fea import spine

    source = textwrap.dedent(inspect.getsource(spine.process_vertebra))

    assert "bounds_with_reference_extent" not in source
    assert "body_bounds = transformed_body_bounds" in source


def test_spine_projected_disks_pad_to_contact_plane_bounds():
    pytest.importorskip("vtk")

    from ogo.fea import spine

    source = textwrap.dedent(inspect.getsource(spine.process_vertebra))

    assert "projected_material_disk_required_bounds" in source
    assert "pad_vtk_images_to_physical_bounds" in source


def test_spine_boundary_conditions_fix_inferior_cap_in_all_directions():
    vtkbone = pytest.importorskip("vtkbone")

    from ogo.fea.model import apply_spine_boundary_conditions

    class RecordingModel:
        def __init__(self):
            self.calls = []

        def ApplyBoundaryCondition(self, *args):
            self.calls.append(args)

    model = RecordingModel()
    apply_spine_boundary_conditions(model, fe_displacement=-0.2924)

    assert model.calls == [
        ("body_top", vtkbone.vtkboneConstraint.SENSE_Z, -0.2924, "top_displacement"),
        ("body_bottom", vtkbone.vtkboneConstraint.SENSE_X, 0, "bottom_fixed_x"),
        ("body_bottom", vtkbone.vtkboneConstraint.SENSE_Y, 0, "bottom_fixed_y"),
        ("body_bottom", vtkbone.vtkboneConstraint.SENSE_Z, 0, "bottom_fixed_z"),
    ]


def test_spine_crop_to_bounding_box_updates_origin_and_extent():
    np = pytest.importorskip("numpy")
    vtk = pytest.importorskip("vtk")
    from vtk.util.numpy_support import numpy_to_vtk

    from ogo.fea.spine import crop_to_bounding_box

    values = np.zeros((6, 7, 8), dtype=np.uint8)
    values[2:5, 3:6, 4:7] = 1
    image = vtk.vtkImageData()
    image.SetDimensions(values.shape)
    image.SetOrigin(10.0, 20.0, 30.0)
    image.SetSpacing(2.0, 3.0, 4.0)
    image.GetPointData().SetScalars(numpy_to_vtk(values.ravel(order="F"), deep=True))

    cropped, bounds = crop_to_bounding_box(image)
    output = cropped.GetOutput()

    assert bounds == [2, 4, 3, 5, 4, 6]
    assert output.GetExtent() == (0, 2, 0, 2, 0, 2)
    assert output.GetOrigin() == pytest.approx((14.0, 29.0, 46.0))


def test_spine_compression_read_preserves_reader_get_output_contract(tmp_path):
    np = pytest.importorskip("numpy")
    sitk = pytest.importorskip("SimpleITK")
    pytest.importorskip("vtk")

    from ogo.fea.spine import read

    image_path = tmp_path / "image.nii.gz"
    image = sitk.GetImageFromArray(np.ones((2, 3, 4), dtype=np.uint8))
    sitk.WriteImage(image, str(image_path))

    reader = read(str(image_path))

    assert reader.GetOutput().GetDimensions() == (4, 3, 2)


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
