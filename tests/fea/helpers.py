"""Fixtures and builders for FEA regression tests."""

from pathlib import Path

import numpy as np
import pytest
import vtk
import vtkbone
from vtk.util.numpy_support import numpy_to_vtkIdTypeArray, vtk_to_numpy

from ogo.fea.spine import (
    DEFAULT_SPINE_POISSONS_RATIO,
    BENCHMARK_PMMA_MATERIAL_ID,
    DEFAULT_SPINE_PMMA_E_MPA,
    DEFAULT_SPINE_PMMA_POISSONS_RATIO,
    BENCHMARK_LINEAR_FE_DISPLACEMENT_MM,
    BENCHMARK_LINEAR_TARGET_DISPLACEMENT_PERCENT,
    BENCHMARK_NONLINEAR_FE_DISPLACEMENT_MM,
    BENCHMARK_NONLINEAR_TARGET_DISPLACEMENT_PERCENT,
    read,
    convert_image_to_material,
    merge_vtk_images,
)
from tests.fea.test_gt_disk_crop import vtk_image


@pytest.fixture
def hip_model():
    import ogo.util.Helper as helper

    data = np.zeros((4, 8, 12), dtype=np.uint16)
    data[1:3, 2:6, 1:11] = 1
    data[1:3, 1:2, 7:10] = 5000
    data[1:3, 6:7, 8:11] = 5000
    grid = helper.Image2Mesh(vtk_image(data))
    model = vtkbone.vtkboneFiniteElementModel()
    model.ShallowCopy(grid)
    xyz = vtk_to_numpy(model.GetPoints().GetData())
    for name, active in (("Femoral_Head_PMMA_Nodes", xyz[:, 1] == xyz[:, 1].max()),
                         ("Greater_Trochanter_PMMA_Nodes", xyz[:, 1] == xyz[:, 1].min()),
                         ("Distal_Femur_Nodes", xyz[:, 2] == xyz[:, 2].min())):
        ids = numpy_to_vtkIdTypeArray(np.flatnonzero(active), deep=True)
        ids.SetName(name)
        model.AddNodeSet(ids)
    return model


def benchmark_linear_params():
    """Return the linear spineFE-benchmark model parameters."""
    return {
        "poissons_ratio": DEFAULT_SPINE_POISSONS_RATIO,
        "pmma_mat_id": BENCHMARK_PMMA_MATERIAL_ID,
        "pmma_E": DEFAULT_SPINE_PMMA_E_MPA,
        "pmma_v": DEFAULT_SPINE_PMMA_POISSONS_RATIO,
        "top_node_set_id": 6,
        "bottom_node_set_id": 5,
        "top_direction": (0, 0, 1),
        "bottom_direction": (0, 0, -1),
        "fe_displacement": BENCHMARK_LINEAR_FE_DISPLACEMENT_MM,
        "target_displacement_percent": BENCHMARK_LINEAR_TARGET_DISPLACEMENT_PERCENT,
        "pmma_yield_compression": None,
        "pmma_yield_tension": None,
        "cort_poissons_ratio": None,
        "elastic_E_func": "kopperdahl_trab_E",
        "yield_comp_func": None,
        "yield_tens_func": None,
        "cort_elastic_E_func": "kopperdahl_trab_E",
        "cort_yield_comp_func": None,
        "cort_yield_tens_func": None,
    }


def benchmark_nonlinear_params():
    """Return the nonlinear spineFE-benchmark model parameters."""
    params = benchmark_linear_params()
    params.update(
        {
            "fe_displacement": BENCHMARK_NONLINEAR_FE_DISPLACEMENT_MM,
            "target_displacement_percent": BENCHMARK_NONLINEAR_TARGET_DISPLACEMENT_PERCENT,
            "pmma_yield_compression": 70.0,
            "pmma_yield_tension": 70.0,
            "yield_comp_func": "kopperdahl_trab_yc",
            "yield_tens_func": "kopperdahl_trab_yc",
            "cort_yield_comp_func": "kopperdahl_trab_yc",
            "cort_yield_tens_func": "kopperdahl_trab_yc",
        }
    )
    return params


def label_mask_from_vtk(vtk_image, condition_fn):
    """Create a binary VTK mask using a NumPy condition over voxel labels."""
    import numpy as np
    import vtk
    from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

    arr = vtk_to_numpy(vtk_image.GetPointData().GetScalars()).reshape(
        vtk_image.GetDimensions(), order="F"
    )
    out_arr = condition_fn(arr).astype(np.uint8)
    out_vtk = vtk_image.NewInstance()
    out_vtk.DeepCopy(vtk_image)
    out_vtk.GetPointData().SetScalars(
        numpy_to_vtk(out_arr.ravel(order="F"), deep=True, array_type=vtk.VTK_UNSIGNED_CHAR)
    )
    return out_vtk


def prepare_benchmark_images(input_image_path, input_mask_path, n_bins=128):
    """Reproduce the spineFE-benchmark notebook image and mask preparation."""
    import numpy as np

    input_image = read(str(input_image_path)).GetOutput()
    input_mask_with_disk = read(str(input_mask_path)).GetOutput()

    cortical_mask = label_mask_from_vtk(input_mask_with_disk, lambda x: np.isin(x, [2, 4]))
    disk_mask = label_mask_from_vtk(input_mask_with_disk, lambda x: np.isin(x, [5, 6]))
    vertebra_mask = label_mask_from_vtk(
        input_mask_with_disk,
        lambda x: np.where(np.isin(x, [1, 2]), 1, np.where(np.isin(x, [3, 4]), 2, 0)),
    )

    binned_image, bin_centers = convert_image_to_material(
        input_image, vertebra_mask, n_bins=n_bins, cort_mask=cortical_mask
    )
    image_with_disk = merge_vtk_images(
        [binned_image, disk_mask],
        [None, 300],
        overwrite_existing=False,
    )
    return image_with_disk, input_mask_with_disk, bin_centers


def build_benchmark_sample_model(input_image_path, input_mask_path, output_model_path, nonlinear=False):
    """Build the public spineFE-benchmark sample model and write it to disk."""
    from ogo.fea.model import create_microfe_model, write_model
    from ogo.util.faim import set_prescribed_displacement_from_percent

    image_with_disk, mask_with_disk, bin_centers = prepare_benchmark_images(
        input_image_path, input_mask_path
    )
    params = benchmark_nonlinear_params() if nonlinear else benchmark_linear_params()
    target_displacement_percent = params.pop("target_displacement_percent")
    model = create_microfe_model(image_with_disk, mask_with_disk, bin_centers, **params)
    write_model(model, output_model_path)
    set_prescribed_displacement_from_percent(
        output_model_path,
        report_profile="spine",
        failure_axis="z",
        target_displacement_percent=target_displacement_percent,
        displacement_sign=-1.0 if params["fe_displacement"] < 0 else 1.0,
    )
    return model


def find_spinefe_benchmark_dir(start_dir=None, env_var="SPINEFE_BENCHMARK_DIR"):
    """Locate the optional public spineFE-benchmark checkout."""
    import os

    env_path = os.environ.get(env_var)
    if env_path:
        candidate = Path(env_path).expanduser()
        if candidate.exists():
            return candidate

    base = Path(start_dir or Path.cwd()).resolve()
    for parent in [base] + list(base.parents):
        candidate = parent / "spineFE-benchmark"
        if candidate.exists():
            return candidate
    return None
