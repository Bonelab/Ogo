"""Spine-specific FE helpers and spineFE benchmark presets."""

from pathlib import Path


def benchmark_linear_params():
    """Return the linear spineFE-benchmark model parameters."""
    return {
        "poissons_ratio": 0.3,
        "pmma_mat_id": 300,
        "pmma_E": 2500,
        "pmma_v": 0.3,
        "top_node_set_id": 6,
        "bottom_node_set_id": 5,
        "top_direction": (0, 0, 1),
        "bottom_direction": (0, 0, -1),
        "fe_displacement": -0.2,
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
            "fe_displacement": -2.0,
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

    from ogo.cli.ref.SpineCompressionFe import convert_image_to_material, merge_vtk_images, read

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
    image_with_disk = merge_vtk_images([binned_image, disk_mask], [None, 300])
    return image_with_disk, input_mask_with_disk, bin_centers


def build_benchmark_sample_model(input_image_path, input_mask_path, output_model_path, nonlinear=False):
    """Build the public spineFE-benchmark sample model and write it to disk."""
    from ogo.fea.model import create_microfe_model, write_model

    image_with_disk, mask_with_disk, bin_centers = prepare_benchmark_images(
        input_image_path, input_mask_path
    )
    params = benchmark_nonlinear_params() if nonlinear else benchmark_linear_params()
    model = create_microfe_model(image_with_disk, mask_with_disk, bin_centers, **params)
    write_model(model, output_model_path)
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

