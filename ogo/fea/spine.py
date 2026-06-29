"""Spine-specific defaults, helpers, and spineFE benchmark presets.

Default spine workflow:
1. Read the calibrated density image and labelled vertebra mask. The high-level
   wrapper uses ``--vertebra LEVEL:BODY_LABEL:PROCESS_LABEL`` so body and
   posterior-process labels are explicit and traceable.
2. Crop the vertebral body, posterior process, calibrated image, and combined
   vertebra mask to the same bounding box.
3. ICP-align the vertebral body to the bundled L4 compression reference. By
   default the reference is PCA-scaled within 0.8,0.8,0.75 to 1.2,1.2,1.3 before
   rigid ICP; users can override the scale/reference explicitly.
4. Apply the ICP transform and set the final output spacing to
   1.0 x 1.0 x 1.0 mm in one shared VTK reslice helper. Image data use cubic
   interpolation; body/process labels use nearest-neighbor interpolation.
5. Smooth the transformed body/process masks with one binary close/open pass
   only when at least one input spacing dimension is coarser than 2 mm. Then
   generate fixed-thickness anatomy PMMA caps on the superior and inferior body
   surfaces. The maintained cap thickness is 10 mm and the intrusion depth is
   6 mm.
6. Build materials with the same shared bone/PMMA material-table helper used by
   the femur workflow. The spine convention is trabecular material IDs 1..128
   and cortical IDs 129..256; PMMA is a separate material ID.
7. Apply axial compression boundary conditions: prescribed displacement at the
   superior PMMA cap toward the inferior cap and fixed displacement at the
   inferior cap. The wrapper solves and reports at 0.68 percent strain by
   default.
"""

from pathlib import Path


DEFAULT_SPINE_ISO_RESOLUTION_MM = 1.0
DEFAULT_SPINE_PMMA_THICKNESS_MM = 10
DEFAULT_SPINE_PMMA_INTRUSION_MM = 6
DEFAULT_SPINE_PMMA_MATERIAL_ID = 5000
DEFAULT_SPINE_POISSONS_RATIO = 0.3
DEFAULT_SPINE_PMMA_E_MPA = 2500
DEFAULT_SPINE_PMMA_POISSONS_RATIO = 0.3
DEFAULT_SPINE_FE_DISPLACEMENT_MM = -1.0
DEFAULT_SPINE_TARGET_DISPLACEMENT_PERCENT = 0.68
DEFAULT_SPINE_MASK_SMOOTHING_SPACING_THRESHOLD_MM = 2.0
DEFAULT_SPINE_TOP_NODE_SET_ID = 4
DEFAULT_SPINE_BOTTOM_NODE_SET_ID = 3
DEFAULT_SPINE_REGISTRATION_SCALE = None
DEFAULT_SPINE_REGISTRATION_MIN_SCALE = "0.8,0.8,0.75"
DEFAULT_SPINE_REGISTRATION_MAX_SCALE = "1.2,1.2,1.3"
DEFAULT_SPINE_REGISTRATION_BACKEND = "vtk"
DEFAULT_SPINE_REFERENCE_FILENAME = "L4_BODY_SPINE_COMPRESSION_REF.vtk"

SPINE_ALIGNMENT_METHOD = "scaled ICP to reference vertebral body"

BENCHMARK_LINEAR_FE_DISPLACEMENT_MM = -0.2
BENCHMARK_NONLINEAR_FE_DISPLACEMENT_MM = -2.0
BENCHMARK_PMMA_MATERIAL_ID = 300


def default_spine_reference_path():
    """Return the bundled reference body used for spine ICP alignment."""
    return Path(__file__).resolve().parents[1] / "dat" / DEFAULT_SPINE_REFERENCE_FILENAME


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

    from ogo.fea.spine_compression import convert_image_to_material, merge_vtk_images, read

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
