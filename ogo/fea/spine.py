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
DEFAULT_SPINE_LABEL_SMOOTHING_SIGMA_MM = 1.0
DEFAULT_SPINE_LABEL_SMOOTHING_THRESHOLD = 0.5
SPINE_PREPROCESSING_CROP_MARGIN_MM = 8.0
DEFAULT_SPINE_MASK_SMOOTHING_SPACING_THRESHOLD_MM = 2.0
DEFAULT_SPINE_TOP_NODE_SET_ID = 4
DEFAULT_SPINE_BOTTOM_NODE_SET_ID = 3
DEFAULT_SPINE_REGISTRATION_SCALE = None
DEFAULT_SPINE_REGISTRATION_MIN_SCALE = "0.8,0.8,0.75"
DEFAULT_SPINE_REGISTRATION_MAX_SCALE = "1.2,1.2,1.3"
DEFAULT_SPINE_REGISTRATION_BACKEND = "vtk"
DEFAULT_SPINE_REFERENCE_FILENAME = "L4_BODY_SPINE_COMPRESSION_REF.vtk"
SPINE_CONTACT_SIZE_FRACTION = (1.6, 1.6)
SPINE_SUPERIOR_CONTACT_CENTER_FRACTION = (0.5, 0.5, 1.05)
SPINE_INFERIOR_CONTACT_CENTER_FRACTION = (0.5, 0.5, -0.05)

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


# -----------------------------------------------------------------------------
# Spine compression workflow builder
# -----------------------------------------------------------------------------

import os
import argparse
import os
from glob import glob
import ogo.fea.material_laws

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import SimpleITK as sitk
from matplotlib import pyplot as plt
from scipy.ndimage import (
    binary_dilation,
    binary_erosion,
    find_objects,
    label,
)
from skimage.exposure import rescale_intensity
import vtk
import vtkbone
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

from ogo.fea.alignment import (
    estimate_rigid_icp,
    invert_point_transform,
    output_grid_for_point_transform,
    point_cloud_axis_lengths,
    polydata_points,
    resample_vtk_image_with_point_transform,
    sample_points,
    surface_points_from_vtk_mask,
)
from ogo.fea.boundary import (
    bbox_relative_contact_plane,
    foreground_voxel_center_bounds,
    generate_projected_material_disk_vtk,
    pad_vtk_images_to_physical_bounds,
    projected_material_disk_required_bounds,
    resample_vtk_image_to_spacing,
    should_smooth_resampled_mask,
    smooth_binary_mask_vtk,
    smooth_label_mask_vtk,
)
from ogo.fea.model import (
    apply_spine_boundary_conditions as apply_boundary_conditions,
    create_microfe_model,
    find_and_add_visible_nodes,
    log_fe_arguments,
)
import ogo.util.Helper as ogo
from ogo.util.echo_arguments import echo_arguments
from scipy.ndimage import gaussian_filter


vtk.vtkObject.GlobalWarningDisplayOff()


def _as_bool(value):
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"1", "true", "t", "yes", "y", "on"}:
        return True
    if text in {"0", "false", "f", "no", "n", "off"}:
        return False
    raise argparse.ArgumentTypeError(f"Invalid boolean value: {value}")


###################################################################### HELPERS (Python)
## General helper functions:
###################################################################### HELPERS (Python)

def remove_extension(filename):
    while True:
        filename, ext = os.path.splitext(filename)
        if not ext:
            break
    return filename

def print_matrix(matrix):
    for i in range(4):
        # Create a string for the entire row
        row_message = " ".join(f"{matrix.GetElement(i, j):0.4f}" for j in range(4))
        # Send the row as a single message
        ogo.message(row_message)


def parse_scale_triplet(value, name="scale"):
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return np.array([float(value)] * 3, dtype=float)
    parts = [part.strip() for part in str(value).replace("x", ",").split(",") if part.strip()]
    try:
        values = [float(part) for part in parts]
    except ValueError as exc:
        raise ValueError(f"{name} must be a float or three comma-separated floats.") from exc
    if len(values) == 1:
        values = values * 3
    if len(values) != 3:
        raise ValueError(f"{name} must be a float or three comma-separated floats.")
    return np.array(values, dtype=float)


###################################################################### HELPERS (VTK)
## VTK Helper functions:
###################################################################### HELPERS (VTK)


class _NiftiImageHandle:
    """Reader-like wrapper for NIfTI images loaded through Ogo's orientation path."""

    def __init__(self, image):
        self._image = image

    def GetOutput(self):
        return self._image


def read(input_mask):
    return _NiftiImageHandle(ogo.readNii(input_mask))

def check_vertebra_presence(mask_reader, vertebra):
    if vertebra not in vtk_to_numpy(mask_reader.GetOutput().GetPointData().GetScalars()).ravel():
        raise ValueError(f'Mask does not contain label {vertebra} for "{vertebra}" vertebra.')

def threshold(mask_output, label):
    threshold = vtk.vtkImageThreshold()
    threshold.SetInputData(mask_output)
    threshold.ThresholdBetween(float(label), float(label))
    threshold.ReplaceInOn()
    threshold.SetInValue(1)
    threshold.ReplaceOutOn()
    threshold.SetOutValue(0)
    threshold.SetOutputScalarTypeToUnsignedChar()
    threshold.Update()
    return threshold

def label_mask(mask, label_value):
    # Use vtkImageMathematics to change 1 values to a specific label value

    caster = vtk.vtkImageCast()
    caster.SetInputData(mask)
    caster.SetOutputScalarTypeToUnsignedChar()
    caster.ClampOverflowOn()  # Optional, clamps values that are out of range of the target data type
    caster.Update()

    math = vtk.vtkImageMathematics()
    math.SetInputData(caster.GetOutput())
    math.SetOperationToMultiplyByK()
    math.SetConstantK(label_value)
    math.Update()
    return math

def combine_mask(mask1, mask2):
    # Combine two masks using vtkImageLogic with a logical OR operation
    logic = vtk.vtkImageLogic()
    logic.SetInput1Data(mask1)
    logic.SetInput2Data(mask2)
    logic.SetOperationToOr()
    logic.Update()
    return logic

def add_masks(mask1, mask2):
    # Combine two masks by adding their values
    math = vtk.vtkImageMathematics()
    math.SetInput1Data(mask1)
    math.SetInput2Data(mask2)
    math.SetOperationToAdd()
    math.Update()

    # Get the result of addition
    added_mask = math.GetOutput()

    # Apply thresholding to ensure the max value is 2
    threshold = vtk.vtkImageThreshold()
    threshold.SetInputData(added_mask)
    threshold.ThresholdByUpper(2)  # Set upper threshold limit to 2
    threshold.SetReplaceIn(True)
    threshold.SetInValue(2)  # Values above the upper threshold will be set to 2
    threshold.SetOutputScalarTypeToUnsignedChar()
    threshold.Update()

    return threshold

def crop_to_bounding_box(image_data, bb=None):
    if bb is None:

        extents = image_data.GetExtent()
        scalars = vtk_to_numpy(image_data.GetPointData().GetScalars())
        reshaped_scalars = scalars.reshape((extents[5]-extents[4]+1, extents[3]-extents[2]+1, extents[1]-extents[0]+1))

        nz = np.nonzero(reshaped_scalars)

        if not nz[0].size:
            raise ValueError("No non-zero voxels found after thresholding. Check earlier steps.")

        min_x, max_x = np.min(nz[2])+extents[0], np.max(nz[2])+extents[0]
        min_y, max_y = np.min(nz[1])+extents[2], np.max(nz[1])+extents[2]
        min_z, max_z = np.min(nz[0])+extents[4], np.max(nz[0])+extents[4]
        bb = [min_x, max_x, min_y, max_y, min_z, max_z]
    else:
        min_x, max_x, min_y, max_y, min_z, max_z = bb


    ogo.message(f"Cropping to extents: X({min_x} to {max_x}), Y({min_y} to {max_y}), Z({min_z} to {max_z})")

    crop = vtk.vtkImageClip()
    crop.SetInputData(image_data)
    crop.SetOutputWholeExtent(min_x, max_x, min_y, max_y, min_z, max_z)
    crop.ClipDataOn()
    crop.Update()

    if crop.GetOutput().GetNumberOfPoints() == 0:
        raise ValueError("Cropping resulted in an empty image. Extents might be incorrect.")

    output = crop.GetOutput()
    spacing = output.GetSpacing()
    origin = image_data.GetOrigin()
    output.SetOrigin(
        float(origin[0]) + float(min_x) * float(spacing[0]),
        float(origin[1]) + float(min_y) * float(spacing[1]),
        float(origin[2]) + float(min_z) * float(spacing[2]),
    )
    dims = output.GetDimensions()
    output.SetExtent(0, dims[0] - 1, 0, dims[1] - 1, 0, dims[2] - 1)

    return crop, bb

def pad_vtk_image(vtk_image, axis='x', pmma_thick=5, pad_value=0):
    """
    Pad a VTK image along the specified axis with a given pmma_thick while preserving VTK metadata.

    Parameters:
    - vtk_image: The input VTK image to pad.
    - axis: The axis along which to pad ('x', 'y', or 'z').
    - pmma_thick: Number of voxels to pad on each side along the specified axis.
    - pad_value: The value to pad with (default: 0).

    Returns:
    - padded_vtk_image: The padded VTK image with updated metadata.
    """
    # Map the axis to the corresponding dimension index
    axis_map = {'x': 0, 'y': 1, 'z': 2}
    if axis not in axis_map:
        raise ValueError("Axis must be one of 'x', 'y', or 'z'.")
    extrusion_axis = axis_map[axis]

    # Convert VTK image data to a NumPy array
    input_array = vtk_to_numpy(vtk_image.GetPointData().GetScalars()).reshape(
        vtk_image.GetDimensions(), order='F'
    )

    # Define the padding configuration
    pad_width = [(0, 0)] * 3  # Default no padding
    pad_width[extrusion_axis] = (pmma_thick, pmma_thick)  # Apply padding along the target axis

    # Apply padding
    padded_array = np.pad(input_array, pad_width=pad_width, mode='constant', constant_values=pad_value)

    # Update VTK metadata
    padded_origin = list(vtk_image.GetOrigin())
    padded_dimensions = list(padded_array.shape)
    spacing = vtk_image.GetSpacing()

    # Adjust the origin for padding in the negative direction
    padded_origin[extrusion_axis] -= pmma_thick * spacing[extrusion_axis]

    # Convert the padded array back to VTK format
    padded_vtk_image = vtk_image.NewInstance()
    padded_vtk_image.SetDimensions(padded_dimensions)
    padded_vtk_image.SetOrigin(padded_origin)
    padded_vtk_image.SetSpacing(spacing)
    padded_vtk_image.GetPointData().SetScalars(numpy_to_vtk(padded_array.ravel(order='F')))

    return padded_vtk_image


def merge_vtk_images(image_list, label_list=None, overwrite_existing=True):
    """
    Merges an arbitrary list of vtkImageData objects into a single vtkImageData.
    Optionally assigns specific label values for each image.

    Parameters:
    - image_list: List of vtkImageData objects to merge.
    - label_list: Optional list of label values corresponding to each image.
                  If None or an entry in the list is None, uses the original values in the image.
    - overwrite_existing: If False, labeled images only write into zero-valued voxels already
                          present in the merged image.

    Returns:
    - merged_image: vtkImageData containing the merged result.
    """
    if not image_list:
        raise ValueError("The image list is empty.")

    if label_list is not None and len(label_list) != len(image_list):
        raise ValueError("The label list must be the same length as the image list.")

    # Get metadata from the first image
    reference_image = image_list[0]
    dims = reference_image.GetDimensions()
    spacing = reference_image.GetSpacing()
    origin = reference_image.GetOrigin()

    # Initialize the merged array as zeros with the same dimensions as the reference image
    merged_array = np.zeros(dims, dtype=np.float32)

    # Merge each image into the merged array
    for idx, image in enumerate(image_list):
        image_array = vtk_to_numpy(image.GetPointData().GetScalars()).reshape(dims, order="F")
        if label_list is not None and label_list[idx] is not None:
            label_value = label_list[idx]
            # Replace non-zero values in the image with the label value
            write_mask = image_array > 0
            if not overwrite_existing:
                write_mask &= merged_array == 0
            merged_array[write_mask] = label_value
        else:
            # Add the original values of the image where label is None
            merged_array += image_array

    # Convert the merged array back to VTK format
    merged_vtk_array = numpy_to_vtk(merged_array.ravel(order="F"), deep=True)

    # Create a new vtkImageData to hold the merged result
    merged_image = vtk.vtkImageData()
    merged_image.SetDimensions(dims)
    merged_image.SetSpacing(spacing)
    merged_image.SetOrigin(origin)
    merged_image.GetPointData().SetScalars(merged_vtk_array)

    return merged_image

###################################################################### REGISTRATION
## functions related to ICP registration
###################################################################### REGISTRATION

def _expand_bounding_box_by_margin(bb, image_data, margin_mm):
    spacing = image_data.GetSpacing()
    extent = image_data.GetExtent()
    margins = [
        int(np.ceil(float(margin_mm) / max(float(spacing[0]), 1.0e-12))),
        int(np.ceil(float(margin_mm) / max(float(spacing[1]), 1.0e-12))),
        int(np.ceil(float(margin_mm) / max(float(spacing[2]), 1.0e-12))),
    ]
    return [
        max(int(extent[0]), int(bb[0]) - margins[0]),
        min(int(extent[1]), int(bb[1]) + margins[0]),
        max(int(extent[2]), int(bb[2]) - margins[1]),
        min(int(extent[3]), int(bb[3]) + margins[1]),
        max(int(extent[4]), int(bb[4]) - margins[2]),
        min(int(extent[5]), int(bb[5]) + margins[2]),
    ]


def crop_and_transform(fullvertebra, body_output, process_output, image, label_output=None, *, margin_mm=0.0):

    _, bb = crop_to_bounding_box(fullvertebra.GetOutput())
    if float(margin_mm) > 0.0:
        bb = _expand_bounding_box_by_margin(bb, fullvertebra.GetOutput(), margin_mm)
    isolated_vertebra, _ = crop_to_bounding_box(body_output, bb)
    isolated_process, _ = crop_to_bounding_box(process_output, bb)
    isolated_image, _ = crop_to_bounding_box(image, bb)
    if label_output is None:
        return isolated_vertebra, isolated_process, isolated_image
    isolated_label, _ = crop_to_bounding_box(label_output, bb)
    return isolated_vertebra, isolated_process, isolated_image, isolated_label

def perform_marching_cubes(body_output):
    mcubes = vtk.vtkImageMarchingCubes()
    mcubes.SetInputConnection(body_output.GetOutputPort())
    mcubes.SetValue(1, 1.0)
    mcubes.Update()
    if mcubes.GetOutput().GetNumberOfPoints() == 0:
        raise ValueError("No data was generated by the Marching Cubes algorithm.")
    return mcubes


def get_icp_with_scaling(
    body,
    reference_path,
    scale_factors=None,
    min_scale=(0.8, 0.8, 0.75),
    max_scale=(1.2, 1.2, 1.3),
):

    sample_surface_points = surface_points_from_vtk_mask(
        body.GetOutput(),
        max_points=8000,
        sample_mode="stride",
    )
    reference_polydata = ogo.readPolyData(reference_path)
    reference_points = sample_points(
        polydata_points(reference_polydata),
        max_points=8000,
        mode="linspace",
    )

    sample_lengths = point_cloud_axis_lengths(sample_surface_points)
    ref_lengths = point_cloud_axis_lengths(reference_points)
    if scale_factors is None:
        scale_factors = sample_lengths / ref_lengths
        min_scale = parse_scale_triplet(min_scale, "registration_min_scale")
        max_scale = parse_scale_triplet(max_scale, "registration_max_scale")
        scale_factors = np.minimum(scale_factors, max_scale)
        scale_factors = np.maximum(scale_factors, min_scale)
        scale_source = "voxel-surface PCA axis lengths"
    else:
        scale_factors = parse_scale_triplet(scale_factors, "registration_scale")
        scale_source = "manual override"

    ogo.message(f"PCA axis lengths (sample): {np.round(sample_lengths, 2)}")
    ogo.message(f"PCA axis lengths (reference): {np.round(ref_lengths, 2)}")
    ogo.message(f"Scale factors applied to reference ({scale_source}): {np.round(scale_factors, 3)}")

    transform = estimate_rigid_icp(
        moving_points=reference_points * scale_factors,
        fixed_points=sample_surface_points,
        iterations=50,
        tolerance=1.0e-4,
        start_by_matching_centroids_only=False,
        convergence="delta",
        distance_mode="mean",
    )

    matrix = vtk.vtkMatrix4x4()
    matrix.Identity()
    rotation = transform["rotation"]
    translation = transform["translation"]
    for row in range(3):
        for col in range(3):
            matrix.SetElement(row, col, float(rotation[row, col]))
        matrix.SetElement(row, 3, float(translation[row]))

    vtk_transform = vtk.vtkTransform()
    vtk_transform.SetMatrix(matrix)
    ogo.message(
        "ICP (with scaled reference) iterations=%d mean_distance=%0.4f"
        % (transform["iterations"], transform["mean_distance"])
    )
    ogo.message("ICP reference-to-sample Matrix:")
    print_matrix(vtk_transform.GetMatrix())

    return vtk_transform

def get_icp(body, reference_path):

    mcubes = perform_marching_cubes(body)
    mcubes_output = mcubes.GetOutput()
    reference_bone = ogo.readPolyData(reference_path)

    icp = vtk.vtkIterativeClosestPointTransform()
    icp.SetTarget(mcubes_output)
    icp.SetSource(reference_bone)

    icp.StartByMatchingCentroidsOn()
    icp.GetLandmarkTransform().SetModeToRigidBody()
    icp.SetMeanDistanceModeToRMS()
    icp.SetMaximumMeanDistance(0.05)
    icp.CheckMeanDistanceOn()
    icp.SetMaximumNumberOfLandmarks(250)
    icp.SetMaximumNumberOfIterations(75)
    icp.Update()

    ogo.message("ICP Matrix:")
    print_matrix(icp.GetMatrix())
    return icp

def transform_resample(image, matrix, iso_resolution, interpolation='cubic'):
    output = ogo.transformResample(image, matrix, iso_resolution, interpolation=interpolation)
    if output.GetNumberOfPoints() == 0:
        ogo.message("Reslice output contains no points. Check the input and transformation.")
    return output


def _matrix_rotation_translation(matrix):
    rotation = np.asarray(
        [[matrix.GetElement(row, col) for col in range(3)] for row in range(3)],
        dtype=float,
    )
    translation = np.asarray([matrix.GetElement(row, 3) for row in range(3)], dtype=float)
    return rotation, translation

###################################################################### DISK GENERATION
## Functions related to disk generation (image processing)
###################################################################### DISK GENERATION

def identify_boundary_surface(input_vtk_image, top_value, bottom_value):
    """Transforms a VTK image to highlight the topmost and bottommost components.

    Parameters:
        input_vtk_image (vtk.vtkImageData): The input binary mask as a VTK image.
        top_value (int): The value to assign to the top mask.
        bottom_value (int): The value to assign to the bottom mask.

    Returns:
        vtk.vtkImageData: The modified VTK image with the top and bottom masks highlighted.
    """
    # Convert VTK image data to NumPy array
    input_array = vtk_to_numpy(input_vtk_image.GetPointData().GetScalars()).reshape(input_vtk_image.GetDimensions(), order='F')
    input_array = np.swapaxes(input_array, 0, 2)  # Swap axes for processing along correct dimension

    binary_mask = np.copy(input_array) == 2

    # Apply boundary conditions from the original method
    eroded_mask = binary_erosion(binary_mask, iterations=1)
    #dilated_mask = binary_dilation(binary_mask, iterations=1)
    outline_mask = (eroded_mask == 0) & (binary_mask == 1)

    # Highlight the outline
    struct_element_erosion = np.array([[[0, 1, 0], [1, 1, 1], [0, 1, 0]]])
    highlighted_mask = binary_erosion(outline_mask, structure=struct_element_erosion)

    struct_element_dilation = np.zeros((3, 3, 3))
    struct_element_dilation[1, 1, :] = 1
    highlighted_mask = binary_dilation(highlighted_mask, iterations=2, structure=struct_element_dilation)
    highlighted_mask = binary_dilation(highlighted_mask, iterations=1)
    highlighted_mask[binary_mask == 0] = 0

    labeled_mask, num_features = label(highlighted_mask)
    sizes = np.bincount(labeled_mask.ravel())[1:]
    sorted_indices = np.argsort(sizes)[::-1]
    object_slices = find_objects(labeled_mask)

    top_positions = [(i + 1, object_slices[i][0].start) for i in sorted_indices[:2]]
    top_positions.sort(key=lambda x: x[1])
    top_component, bottom_component = top_positions

    top_mask = binary_dilation(labeled_mask == top_component[0], iterations=2) & outline_mask
    bottom_mask = binary_dilation(labeled_mask == bottom_component[0], iterations=2) & outline_mask

    # Assign values and convert back to VTK
    result_array = input_array.copy()
    result_array[top_mask] = top_value
    result_array[bottom_mask] = bottom_value
    result_array[input_array==0]=0
    result_array = np.swapaxes(result_array, 0, 2)  # Swap axes for processing along correct dimension

    vtk_data_array = numpy_to_vtk(result_array.ravel(order='F'), deep=True)
    output_vtk_image = vtk.vtkImageData()
    output_vtk_image.DeepCopy(input_vtk_image)
    output_vtk_image.GetPointData().SetScalars(vtk_data_array)

    return output_vtk_image

###################################################################### MICRO-FE (FAIM)
## Functions related to micro-FE model
###################################################################### MICRO-FE (FAIM)

def generate_cortical_mask(image_vtk, mask_vtk, threshold=0.2, min_th=1, max_th=3):
    """
    Generate a cortical mask from the input VTK image and mask using a threshold.
    Ensures the cortical mask is at least as thick as the erosion shell.

    Parameters:
    - image_vtk: The input image data (VTK Image Data).
    - mask_vtk: The input mask data (VTK Image Data).
    - threshold: The threshold value for binarization (default: 0.2).
    - min_th: Minimum thickness constraint (default: 1 voxel).
    - max_th: Maximum thickness constraint (default: 3 voxels).

    Returns:
    - cortical_mask_vtk: The binary mask for the cortical region as a VTK Image Data.
    """

    # Convert VTK images to numpy arrays
    image_np = vtk_to_numpy(image_vtk.GetPointData().GetScalars()).reshape(image_vtk.GetDimensions(), order='F').astype(np.float32)
    #image_np = gaussian_filter(image_np, sigma=1)  # Adjust sigma as needed

    mask_np = vtk_to_numpy(mask_vtk.GetPointData().GetScalars()).reshape(mask_vtk.GetDimensions(), order='F') > 0

    # Apply thresholding to the image (convert to g/ccm by dividing by 1000)
    binary_image = (image_np / 1000) > threshold

    # Compute cortical mask as intersection of thresholded image and mask
    cortical_mask = binary_image & mask_np

    ogo.message(f'Average Image intensity for cort mask {np.mean(image_np[cortical_mask])}')
    ogo.message(f'Average Image intensity for trab mask {np.mean(image_np[~cortical_mask & mask_np])}')


    # Compute the shell (minimum thickness enforcement)
    if min_th>0:
        min_erosion = binary_erosion(mask_np, iterations=min_th, border_value=0)
        min_shell = mask_np & ~min_erosion  # Shell is the difference between the mask and its eroded version
    else:
        min_shell = np.zeros_like(mask_np)

    # Apply maximum thickness constraint
    if max_th>0:
        max_erosion = binary_erosion(mask_np, iterations=max_th, border_value=0)
        max_shell = mask_np & ~max_erosion  # Shell is the difference between the mask and its eroded version
    else:
        max_shell = np.ones_like(mask_np)

    # Ensure cortical mask is at least the shell
    cortical_mask = (cortical_mask | min_shell) & max_shell  # Ensure at least the shell is included

    # Convert the binary mask back to VTK format
    cortical_mask_vtk = vtk.vtkImageData()
    cortical_mask_vtk.DeepCopy(mask_vtk)  # Copy spatial properties

    # Convert numpy array back to VTK format
    mask_vtk_array = numpy_to_vtk(cortical_mask.flatten(order='F'), deep=True, array_type=vtk.VTK_UNSIGNED_CHAR)
    cortical_mask_vtk.GetPointData().SetScalars(mask_vtk_array)

    ogo.message(f'Average Image intensity for cort mask {np.mean(image_np[cortical_mask])}')
    ogo.message(f'Average Image intensity for trab mask {np.mean(image_np[~cortical_mask & mask_np])}')


    return cortical_mask_vtk


def convert_image_to_material(image, mask, n_bins=128, cort_mask=None):
    change = ogo.prepareFiniteElementImage(image)
    mask_change = ogo.prepareFiniteElementImage(mask)
    cort_mask_change = ogo.prepareFiniteElementImage(cort_mask) if cort_mask is not None else None
    connected_image = ogo.imageConnectivity(change)
    thr_image = ogo.bmd_preprocess(connected_image, -31)
    #this was done previously - but if someone wants this they should just do it in the calib functions
    #ash_image = ogo.bmd_K2hpo4ToAsh(thr_image)
    cast_image = ogo.cast2short(thr_image)
    bone_image = ogo.applyMaskByArray(cast_image, mask_change)
    binned_image, bin_centers = ogo.density2materialID(bone_image, n_bins=n_bins, cort_mask=cort_mask_change)

    return binned_image, bin_centers

def resolve_func(func_or_name, module):
    """Resolve a material-law function by object or configured name."""
    from ogo.fea.materials import resolve_material_func

    return resolve_material_func(func_or_name, module)

###################################################################### QUALITY CONTROL
## Functions to check image and boundary conditions
###################################################################### QUALITY CONTROL

## Functions to check the results
def check_image_values(image):
    array = vtk.util.numpy_support.vtk_to_numpy(image.GetPointData().GetScalars())
    unique_values = np.unique(array)
    ogo.message("Unique values in the image:", unique_values)

def calculate_features(image, label):
    label_map = sitk.BinaryThreshold(image, lowerThreshold=label, upperThreshold=label, insideValue=1, outsideValue=0)
    stats = sitk.LabelShapeStatisticsImageFilter()
    stats.Execute(label_map)
    volume = stats.GetPhysicalSize(1)
    centroid = stats.GetCentroid(1)
    principal_axes = stats.GetPrincipalAxes(1)
    extent = stats.GetBoundingBox(1)
    return {'volume': volume, 'centroid': centroid, 'principal_axes': principal_axes, 'extent': extent}

def check_values(features, checks):
    results = {}
    pass_all = True
    for check, condition in checks.items():
        result = eval(condition, {"features": features})
        results[check] = result
        if not result:
            pass_all = False
    results['pass'] = pass_all
    return results

def parse_filename(filepath):
    base_name = os.path.splitext(os.path.basename(filepath))[0]
    parts = base_name.split('_')
    return {
        'ID': parts[0] if len(parts) > 0 else '',
        'TREATMENT': parts[1] if len(parts) > 1 else '',
        'LOCATION': parts[2] if len(parts) > 2 else '',
        'NUMBER': parts[4] if len(parts) > 4 else '',
        'filename': '_'.join(parts[:-3]) if len(parts) > 3 else base_name,
    }

def visualize_slice(image, filepath):
    # Get the dimensions of the vtkImageData object
    dimensions = image.GetDimensions()  # (x, y, z)
    mid_x = dimensions[0] // 2  # Middle slice along the X-axis

    # Convert vtkImageData to a NumPy array
    vtk_array = vtk_to_numpy(image.GetPointData().GetScalars())
    numpy_array = vtk_array.reshape(dimensions[::-1])  # Reverse dimensions to match NumPy order (z, y, x)

    # Extract the middle X slice
    slice_img = numpy_array[:, :, mid_x]

    # Clip values at the 90th percentile
    percentile_90 = np.percentile(slice_img, 90)
    clipped_slice = np.clip(slice_img, 0, percentile_90)  # Cap values at the 90th percentile

    # Normalize the clipped slice to [0, 1] range
    normalized_slice = (clipped_slice - clipped_slice.min()) / (clipped_slice.max() - clipped_slice.min() + 1e-8)

    # Plot the normalized slice
    plt.figure(figsize=(4.0, 4.0))
    plt.imshow(normalized_slice, cmap='viridis', origin="lower")
    plt.title('Middle X', fontsize=9, pad=2)
    plt.xlabel('Y', fontsize=8)
    plt.ylabel('Z', fontsize=8)
    plt.tick_params(axis='both', which='major', labelsize=8, length=2.5, pad=1.5)

    # Save the output image
    plt.tight_layout(pad=0.2)
    plt.savefig(filepath, dpi=300)
    plt.close()
    ogo.message(f"Slice image saved to {filepath}")

def check_image(vtkimage, output_filename=None):
    try:
        # Assuming 'image' is a vtkImageData object
        numpy_array = vtk_to_numpy(vtkimage.GetPointData().GetScalars())
        numpy_array = numpy_array.reshape(vtkimage.GetDimensions(), order='F')
        numpy_array = np.transpose(numpy_array, (2, 1, 0))  # Adjust depending on your specific data
        image = sitk.GetImageFromArray(numpy_array)
        image.SetSpacing(vtkimage.GetSpacing())
        image.SetOrigin(vtkimage.GetOrigin())

        features = {label: calculate_features(image, label) for label in range(1, 5)}

        checks = {
            "VCHECK_PROCESS": "features[1]['volume'] > 4000",
            "VCHECK_BODY": "features[2]['volume'] > 3000",
            "VCHECK_BCBOT": "features[3]['volume'] > 200",
            "VCHECK_BCTOP": "features[4]['volume'] > 200",
            "DCHECK_BC": "features[4]['centroid'][2] > features[2]['centroid'][2] > features[3]['centroid'][2]",
            "EXTENT_Z_BCBOT": "features[3]['extent'][5] < 20",
            "EXTENT_Z_BCTOP": "features[4]['extent'][5] < 20"
        }

        results = check_values(features, checks)
        metadata = parse_filename(output_filename)
        data = {**metadata, **{f'{label}_{key}': val for label, feats in features.items() for key, val in feats.items()}, **results}

        df = pd.DataFrame([data])

        if output_filename is not None:
            check_ending = f"{data['pass']}.csv"
            output_filename = output_filename.replace('.csv', check_ending)
            ogo.message(f"Saving QC results to {output_filename}")
            df.to_csv(output_filename, index=False)
        else:
            ogo.message(df)

    except Exception as e:
        ogo.message(e)

        data = {
            "VCHECK_BODY": False,
            "VCHECK_PROCESS": False,
            "VCHECK_BCBOT": False,
            "VCHECK_BCTOP": False,
            "DCHECK_BC": False,
            "EXTENT_Z_BCBOT": False,
            "EXTENT_Z_BCTOP": False,
            "pass": False
        }

        df = pd.DataFrame([data])
        if output_filename is not None:
            check_ending = f"{data['pass']}.csv"
            output_filename = output_filename.replace('.csv', check_ending)
            ogo.message(f"Saving QC results to {output_filename}")
            df.to_csv(output_filename, index=False)
        else:
            ogo.message(df)

    return data['pass']

def resample_to_match(target_image, source_image, interpolation='nearest'):
    resampler = vtk.vtkImageResample()
    resampler.SetInputData(source_image)

    spacing = target_image.GetSpacing()
    resampler.SetAxisOutputSpacing(0, spacing[0])
    resampler.SetAxisOutputSpacing(1, spacing[1])
    resampler.SetAxisOutputSpacing(2, spacing[2])

    if interpolation == 'nearest':
        resampler.SetInterpolationModeToNearestNeighbor()
    elif interpolation == 'linear':
        resampler.SetInterpolationModeToLinear()
    elif interpolation == 'cubic':
        resampler.SetInterpolationModeToCubic()

    resampler.Update()
    return resampler

def export_nifti_outputs(
    image_vtk,
    body_labeled_vtk,
    process_labeled_vtk,
    cortical_mask_vtk,
    inferior_disk_vtk,
    superior_disk_vtk,
    output_path,
):
    """
    Export resampled grayscale image and labeled segmentation for FE input.

    Labels:
    1 = trabecular body
    2 = cortical body
    3 = trabecular process
    4 = cortical process
    5 = inferior disk
    6 = superior disk
    """
    import SimpleITK as sitk
    import numpy as np
    from vtk.util.numpy_support import vtk_to_numpy
    import os

    def vtk_to_sitk(vtk_image):
        dims = vtk_image.GetDimensions()
        np_array = vtk_to_numpy(vtk_image.GetPointData().GetScalars()).reshape(dims[::-1])  # (z, y, x)
        sitk_image = sitk.GetImageFromArray(np_array)
        sitk_image.SetSpacing(vtk_image.GetSpacing())
        sitk_image.SetOrigin(vtk_image.GetOrigin())
        return sitk_image

    def get_mask_np(vtk_image):
        dims = vtk_image.GetDimensions()
        return vtk_to_numpy(vtk_image.GetPointData().GetScalars()).reshape(dims, order="F")

    ogo.message("Exporting NIfTI files for grayscale image and labeled segmentation...")

    dims = image_vtk.GetDimensions()
    combined_array = np.zeros(dims, dtype=np.uint8)

    # Convert to NumPy
    body_mask = get_mask_np(body_labeled_vtk) > 0
    process_mask = get_mask_np(process_labeled_vtk) > 0
    cortical_mask = get_mask_np(cortical_mask_vtk) > 0
    inferior_mask = get_mask_np(inferior_disk_vtk) > 0
    superior_mask = get_mask_np(superior_disk_vtk) > 0
    im = get_mask_np(image_vtk)


    # Assign labels
    combined_array[(body_mask) & (~cortical_mask)] = 1
    combined_array[(body_mask) & (cortical_mask)] = 2
    combined_array[(process_mask) & (~cortical_mask)] = 3
    combined_array[(process_mask) & (cortical_mask)] = 4
    combined_array[inferior_mask] = 5
    combined_array[superior_mask] = 6
    im[(combined_array>4) | (combined_array==0)] = 0

    # Convert segmentation to SimpleITK
    seg_sitk = sitk.GetImageFromArray(np.transpose(combined_array, (2, 1, 0)))
    seg_sitk.SetSpacing(image_vtk.GetSpacing())
    seg_sitk.SetOrigin(image_vtk.GetOrigin())

    im_sitk = sitk.GetImageFromArray(np.transpose(im, (2, 1, 0)))
    im_sitk.SetSpacing(image_vtk.GetSpacing())
    im_sitk.SetOrigin(image_vtk.GetOrigin())

    outbase = os.path.splitext(output_path)[0]
    im_out = outbase + "_im.nii.gz"
    seg_out = outbase + "_seg.nii.gz"

    sitk.WriteImage(im_sitk, im_out)
    sitk.WriteImage(seg_sitk, seg_out)

    ogo.message(f"Export complete: {im_out}, {seg_out}")

###################################################################### VERTEBRA PIPELINE
# Main pipeline to process a vertebra FE
###################################################################### VERTEBRA PIPELINE


def process_vertebra(input_mask, input_image, n88model_output_path, body_label, process_label, reference_path, **kwargs):

    pmma_mat_id = kwargs.get("pmma_mat_id", 5000)
    iso_resolution = kwargs.get("iso_resolution", DEFAULT_SPINE_ISO_RESOLUTION_MM)
    pmma_thick = kwargs.get("pmma_thick", DEFAULT_SPINE_PMMA_THICKNESS_MM)
    pmma_intrusion = kwargs.get("pmma_intrusion", DEFAULT_SPINE_PMMA_INTRUSION_MM)
    top_node_set_id = kwargs.get("top_node_set_id", DEFAULT_SPINE_TOP_NODE_SET_ID)
    bottom_node_set_id = kwargs.get("bottom_node_set_id", DEFAULT_SPINE_BOTTOM_NODE_SET_ID)
    quality_control = kwargs.get("quality_control", True)
    registration_scale = kwargs.get("registration_scale", DEFAULT_SPINE_REGISTRATION_SCALE)
    registration_min_scale = kwargs.get("registration_min_scale", DEFAULT_SPINE_REGISTRATION_MIN_SCALE)
    registration_max_scale = kwargs.get("registration_max_scale", DEFAULT_SPINE_REGISTRATION_MAX_SCALE)
    mask_smoothing_spacing_threshold = kwargs.get(
        "mask_smoothing_spacing_threshold",
        DEFAULT_SPINE_MASK_SMOOTHING_SPACING_THRESHOLD_MM,
    )
    label_smoothing_sigma_mm = float(
        kwargs.get("label_smoothing_sigma_mm", DEFAULT_SPINE_LABEL_SMOOTHING_SIGMA_MM)
    )
    label_smoothing_threshold = float(
        kwargs.get("label_smoothing_threshold", DEFAULT_SPINE_LABEL_SMOOTHING_THRESHOLD)
    )

    #Read Image and Mask
    ogo.message(f"Reading image...: {input_image}")
    image_reader = read(input_image)
    ogo.message(f"Reading mask...: {input_mask}")
    mask_reader = resample_to_match(image_reader.GetOutput(), read(input_mask).GetOutput())
    input_spacing = image_reader.GetOutput().GetSpacing()

    # Get and Check labels
    ogo.message(f"Checking if labels {body_label} (body) and {process_label} (process) are present...")
    check_vertebra_presence(mask_reader, body_label)
    check_vertebra_presence(mask_reader, process_label)

    # Threshold images to extract body, process and full vertebra
    ogo.message(f"Thresholding...")
    body = threshold(mask_reader.GetOutput(), body_label)
    process = threshold(mask_reader.GetOutput(), process_label)
    fullvertebra = combine_mask(body.GetOutput(), process.GetOutput())

    # Crop everything to same BB
    ogo.message(f"Cropping to common bounding box...")
    (
        isolated_vertebra,
        isolated_process,
        isolated_image,
        isolated_labels,
    ) = crop_and_transform(
        fullvertebra,
        body.GetOutput(),
        process.GetOutput(),
        image_reader.GetOutput(),
        mask_reader.GetOutput(),
        margin_mm=SPINE_PREPROCESSING_CROP_MARGIN_MM,
    )

    ogo.message("Resampling cropped spine inputs to isotropic spacing before ICP...")
    preprocessed_image = resample_vtk_image_to_spacing(
        isolated_image.GetOutput(),
        iso_resolution,
        interpolation="bspline",
    )
    preprocessed_labels = resample_vtk_image_to_spacing(
        isolated_labels.GetOutput(),
        iso_resolution,
        interpolation="nearest",
    )
    if label_smoothing_sigma_mm > 0:
        ogo.message(
            "smoothing isotropic label mask "
            f"(sigma={label_smoothing_sigma_mm} mm, threshold={label_smoothing_threshold})..."
        )
        preprocessed_labels = smooth_label_mask_vtk(
            preprocessed_labels,
            sigma_mm=label_smoothing_sigma_mm,
            threshold=label_smoothing_threshold,
        )
    registration_body = threshold(preprocessed_labels, body_label)
    isolated_vertebra = registration_body
    isolated_process = threshold(preprocessed_labels, process_label)

    # Marching Cubes and Registration
    ogo.message(f"Starting ICP registration to reference...: {reference_path} ")
    icp = get_icp_with_scaling(
        registration_body,
        reference_path,
        scale_factors=registration_scale,
        min_scale=registration_min_scale,
        max_scale=registration_max_scale,
    )

    # Transform Images and resample at the same time (less interpolation)
    reference_to_sample_rotation, reference_to_sample_translation = _matrix_rotation_translation(
        icp.GetMatrix()
    )
    sample_to_reference_rotation, sample_to_reference_translation = invert_point_transform(
        reference_to_sample_rotation,
        reference_to_sample_translation,
    )
    output_spacing = (float(iso_resolution), float(iso_resolution), float(iso_resolution))
    isolated_model_mask = combine_mask(
        isolated_vertebra.GetOutput(),
        isolated_process.GetOutput(),
    )
    model_surface_points = surface_points_from_vtk_mask(
        isolated_model_mask.GetOutput(),
        max_points=None,
    )
    output_origin, output_size = output_grid_for_point_transform(
        model_surface_points,
        rotation=sample_to_reference_rotation,
        translation=sample_to_reference_translation,
        spacing=output_spacing,
        margin_voxels=int(kwargs.get("registration_margin_voxels", 4)),
    )
    transformed_vertebra = resample_vtk_image_with_point_transform(
        isolated_vertebra.GetOutput(),
        rotation=sample_to_reference_rotation,
        translation=sample_to_reference_translation,
        output_origin=output_origin,
        output_size=output_size,
        output_spacing=output_spacing,
        interpolation="nearest",
    )
    transformed_process = resample_vtk_image_with_point_transform(
        isolated_process.GetOutput(),
        rotation=sample_to_reference_rotation,
        translation=sample_to_reference_translation,
        output_origin=output_origin,
        output_size=output_size,
        output_spacing=output_spacing,
        interpolation="nearest",
    )
    transformed_image = resample_vtk_image_with_point_transform(
        preprocessed_image,
        rotation=sample_to_reference_rotation,
        translation=sample_to_reference_translation,
        output_origin=output_origin,
        output_size=output_size,
        output_spacing=output_spacing,
        interpolation="bspline",
    )

    if should_smooth_resampled_mask(input_spacing, mask_smoothing_spacing_threshold):
        ogo.message(
            "smoothing resampled body/process masks because input spacing "
            f"{input_spacing} exceeds {mask_smoothing_spacing_threshold} mm..."
        )
        transformed_vertebra = smooth_binary_mask_vtk(transformed_vertebra, close_iter=1, open_iter=1)
        transformed_process = smooth_binary_mask_vtk(transformed_process, close_iter=1, open_iter=1)
    else:
        ogo.message(
            "skipping body/process mask smoothing because input spacing "
            f"{input_spacing} is <= {mask_smoothing_spacing_threshold} mm in all dimensions..."
        )

    ogo.message(f"relabelling mask and identifying boundary surfaces...")
    # Assuming mask1_data and mask2_data are your initial binary masks
    mask1_labeled = label_mask(transformed_vertebra, 2)  # Apply label 1 to the first mask
    mask2_labeled = label_mask(transformed_process, 1)  # Apply label 2 to the second mask

    # Combine masks here (otherwise it does weird interpolations)
    combined_masks = add_masks(mask1_labeled.GetOutput(), mask2_labeled.GetOutput())  # Combine the labeled masks

    # This creates an image identifying the top and bottom surfaces of the vertebra
    boundary_masks = identify_boundary_surface(combined_masks.GetOutput(), bottom_node_set_id, top_node_set_id)

    ogo.message("creating cortical mask....")
    #1 - 3 alterantive (to create a continous cortical shell - changed to 0 5 )

    cort_mask = generate_cortical_mask(transformed_image, combined_masks.GetOutput(), threshold=1, min_th=2, max_th=5)

    #writer = vtk.vtkXMLImageDataWriter()
    #writer.SetFileName('/home/matthias.walle/work/fem/CT_FE_TEMPLATE/MODELS/10001_QCT_vertebra_20_cortmask.vti')
    #writer.SetInputData(cort_mask)  # For VTK 8+
    #writer.Write()

    ogo.message(f"converting image to material ID...")
    # This converts the raw image to a material ID mapped image
    n_bins = 128
    bone_image, bin_centers = convert_image_to_material(transformed_image, combined_masks.GetOutput(), n_bins=n_bins, cort_mask = cort_mask)

    ogo.message(f"padding images...")
    # Pad images before we add the disks (otherwise they may not have right pmma_thick)
    padded_image = pad_vtk_image(bone_image, axis='x', pmma_thick=pmma_thick, pad_value=0)
    padded_transformed_image = pad_vtk_image(transformed_image, axis='x', pmma_thick=pmma_thick, pad_value=0)
    padded_mask = pad_vtk_image(boundary_masks, axis='x', pmma_thick=pmma_thick, pad_value=0)
    padded_cort_mask = pad_vtk_image(cort_mask, axis='x', pmma_thick=pmma_thick, pad_value=0)
    padded_mask1 = pad_vtk_image(mask1_labeled.GetOutput(), axis='x', pmma_thick=pmma_thick, pad_value=0)
    padded_mask2 = pad_vtk_image(mask2_labeled.GetOutput(), axis='x', pmma_thick=pmma_thick, pad_value=0)

    ogo.message(
        "generating fixed-thickness anatomy body-cap disks "
        f"(pmma_thick = {pmma_thick} mm total, pmma_intrusion = {pmma_intrusion} mm)..."
    )
    transformed_body_bounds = foreground_voxel_center_bounds(padded_mask1)
    body_bounds = transformed_body_bounds
    inferior_plane = bbox_relative_contact_plane(
        body_bounds,
        center_fraction=SPINE_INFERIOR_CONTACT_CENTER_FRACTION,
        size_fraction=SPINE_CONTACT_SIZE_FRACTION,
        projection_axis="z",
        shape="anatomy",
    )
    superior_plane = bbox_relative_contact_plane(
        body_bounds,
        center_fraction=SPINE_SUPERIOR_CONTACT_CENTER_FRACTION,
        size_fraction=SPINE_CONTACT_SIZE_FRACTION,
        projection_axis="z",
        shape="anatomy",
    )
    ogo.message(
        "Inferior PMMA Plane: center=%s normal=%s u=%s v=%s size=%s"
        % (
            inferior_plane["center"],
            inferior_plane["normal"],
            inferior_plane["u_axis"],
            inferior_plane["v_axis"],
            inferior_plane["size"],
        )
    )
    ogo.message(
        "Superior PMMA Plane: center=%s normal=%s u=%s v=%s size=%s"
        % (
            superior_plane["center"],
            superior_plane["normal"],
            superior_plane["u_axis"],
            superior_plane["v_axis"],
            superior_plane["size"],
        )
    )
    required_bounds = [
        projected_material_disk_required_bounds(
            padded_mask1,
            center=plane["center"],
            normal=plane["normal"],
            u_axis=plane["u_axis"],
            v_axis=plane["v_axis"],
            size=plane["size"],
            shape=plane["shape"],
            thickness=pmma_thick,
            intrusion=pmma_intrusion,
        )
        for plane in (inferior_plane, superior_plane)
    ]
    required_bounds = [bounds for bounds in required_bounds if bounds is not None]
    if required_bounds:
        bounds_array = np.asarray(required_bounds, dtype=float)
        desired_bounds = (
            float(bounds_array[:, 0].min()),
            float(bounds_array[:, 1].max()),
            float(bounds_array[:, 2].min()),
            float(bounds_array[:, 3].max()),
            float(bounds_array[:, 4].min()),
            float(bounds_array[:, 5].max()),
        )
        (
            padded_image,
            padded_transformed_image,
            padded_mask,
            padded_cort_mask,
            padded_mask1,
            padded_mask2,
        ), contact_padding = pad_vtk_images_to_physical_bounds(
            [
                padded_image,
                padded_transformed_image,
                padded_mask,
                padded_cort_mask,
                padded_mask1,
                padded_mask2,
            ],
            desired_bounds=desired_bounds,
            constants=[0, 0, 0, 0, 0, 0],
        )
        if any(contact_padding["lower"]) or any(contact_padding["upper"]):
            ogo.message(
                "expanded image canvas for projected disks: lower=%s upper=%s"
                % (contact_padding["lower"], contact_padding["upper"])
            )
    inferior_disk = generate_projected_material_disk_vtk(
        padded_image,
        surface_vtk_image=padded_mask1,
        exclusion_vtk_image=padded_mask1,
        center=inferior_plane["center"],
        normal=inferior_plane["normal"],
        u_axis=inferior_plane["u_axis"],
        v_axis=inferior_plane["v_axis"],
        size=inferior_plane["size"],
        shape=inferior_plane["shape"],
        thickness=pmma_thick,
        intrusion=pmma_intrusion,
        anatomy_constrained=True,
        output_value=1,
    )
    superior_disk = generate_projected_material_disk_vtk(
        padded_image,
        surface_vtk_image=padded_mask1,
        exclusion_vtk_image=padded_mask1,
        center=superior_plane["center"],
        normal=superior_plane["normal"],
        u_axis=superior_plane["u_axis"],
        v_axis=superior_plane["v_axis"],
        size=superior_plane["size"],
        shape=superior_plane["shape"],
        thickness=pmma_thick,
        intrusion=pmma_intrusion,
        anatomy_constrained=True,
        output_value=1,
    )

    ogo.message(f"merging disks with image...")
    # Now we create the image with the disks and the boundary image as well
    image_with_pads = merge_vtk_images(
        [padded_image, superior_disk, inferior_disk],
        [None, pmma_mat_id, pmma_mat_id],
        overwrite_existing=False,
    )

    boundary_masks_with_pads = merge_vtk_images([
        padded_mask, inferior_disk, superior_disk], [None, bottom_node_set_id, top_node_set_id])

    # This is the final model that we will use for the FEA.
    # For the future --> pull out material table to have option without disks
    ogo.message("generating n88model file...")
    model = create_microfe_model(
        image_with_pads,
        boundary_masks_with_pads,
        bin_centers,
        top_boundary_mask_image=superior_disk,
        bottom_boundary_mask_image=inferior_disk,
        top_boundary_mask_label=1,
        bottom_boundary_mask_label=1,
        **kwargs,
    )

    if kwargs.get("export_nifti", False):
        ogo.message("generating nifti files...")

        export_nifti_outputs(
            padded_transformed_image,
            padded_mask1,
            padded_mask2,
            padded_cort_mask,
            inferior_disk,
            superior_disk,
            n88model_output_path
        )

    # Quality control and output. If quality control is activated output only provided if it passes
    if quality_control:
            image_path = n88model_output_path.replace(".n88model", ".png")
            dataframe_path = n88model_output_path.replace(".n88model", "_BCcheck.csv")

            ogo.message(f"Starting QC...")
            test = check_image(boundary_masks, dataframe_path)
            visualize_slice(image_with_pads, image_path)

    ogo.message(f"Writing n88model file: {n88model_output_path}")
    writer = vtkbone.vtkboneN88ModelWriter()
    writer.SetInputData(model)
    writer.SetFileName(n88model_output_path)
    writer.Update()



###################################################################### MAIN
## Main function containing file I/O via argparse
###################################################################### MAIN


def main():
    description = '''
    This script sets up the L4 vertebral compression FE model from the
        density (K2HPO4) calibrated image. This script sets up the model for a L4 vertebra
        (including arch and pedicles). The analysis resamples the image to isotropic voxels,
        transforms the image, applies the bone mask and bins the data. It then creates the FE
        model for solving using FAIM (v8.1, Numerics Solutions Ltd, Calgary, Canada - Steven
        Boyd).

        Input: Calibrated K2HPO4 Image (*.nii), Bone Mask (*_MASK.nii)

        Optional Parameters:
        1) Mask Threshold
        2) Isotropic resample voxel size
        3) Power-law exponent
        4) Power-law coefficient
        5) Bone Poissons ratio
        6) PMMA Elastic Modulus
        7) PMMA Poissons ratio
        8) PMMA pmma_thick
        9) PMMA material ID
        10) FE displacement

        Output: N88 Model (*.n88model)
    '''

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoFEA-spine-builder",
        description=description
    )

    parser.add_argument("calibrated_image", help="*_K2HPO4.nii image file")
    parser.add_argument("bone_mask", help="*_MASK.nii mask image of bone")

    parser.add_argument("--mask_threshold", type=int, required=True,
                        help="Label value for the vertebral body in the bone mask.")
    parser.add_argument("--process_mask_threshold", type=int, required=True,
                        help="Label value for the vertebral process in the bone mask.")
    parser.add_argument("--output_path", type=str, default=None,
                        help="Set output path for the N88 model file. (default: same as input image)")
    parser.add_argument("--quality_control", type=_as_bool, default=True,
                        help="Set quality control flag (visualise output and add volume checks). (default: %(default)s)")
    parser.add_argument("--iso_resolution", type=float, default=DEFAULT_SPINE_ISO_RESOLUTION_MM,
                        help="Set the isotropic voxel size [in mm]. (default: %(default)s [mm])")
    parser.add_argument("--mask_smoothing_spacing_threshold", type=float, default=DEFAULT_SPINE_MASK_SMOOTHING_SPACING_THRESHOLD_MM,
                        help="Smooth resampled body/process masks only when an input spacing dimension exceeds this value. (default: %(default)s [mm])")
    #parser.add_argument("--elastic_exponent", type=float, default=2.29,
    #                    help="Sets the exponent (b) for power law: E=A(den)^b. (default: %(default)s)")
    #parser.add_argument("--elastic_Emax", type=float, default=10500,
    #                    help="Sets the coefficient (A) Elastic Modulus value for the power law: E=A(den)^b. (default: %(default)s [MPa])")
    parser.add_argument("--poissons_ratio", type=float, default=DEFAULT_SPINE_POISSONS_RATIO,
                        help="Sets the Poisson's ratio for the material(s) in the FE model. (default: %(default)s)")
    parser.add_argument("--pmma_E", type=float, default=DEFAULT_SPINE_PMMA_E_MPA,
                        help="Sets the Elastic Modulus for PMMA caps in the FE model. (default: %(default)s [MPa])")
    parser.add_argument("--pmma_v", type=float, default=DEFAULT_SPINE_PMMA_POISSONS_RATIO,
                        help="Sets the Poisson's ratio for the PMMA material(s) in the FE model. (default: %(default)s)")
    parser.add_argument("--pmma_thick", type=int, default=DEFAULT_SPINE_PMMA_THICKNESS_MM,
                        help="Sets the total PMMA cap thickness in the FE model. (default: %(default)s [mm])")
    parser.add_argument("--pmma_intrusion", type=int, default=DEFAULT_SPINE_PMMA_INTRUSION_MM,
                        help="Sets how far anatomy can occupy the fixed PMMA cap thickness. (default: %(default)s [mm])")
    parser.add_argument("--pmma_mat_id", type=int, default=DEFAULT_SPINE_PMMA_MATERIAL_ID,
                        help="Sets the material ID for the PMMA blocks. (default: %(default)s)")
    parser.add_argument("--fe_displacement", type=float, default=DEFAULT_SPINE_FE_DISPLACEMENT_MM,
                        help="Sets the applied displacement in [mm] to the FE model. (default: %(default)s [mm])")
    parser.add_argument("--reference_path", type=str, required=False, default=None,
                        help="Path to the reference vtk file for ICP registration. (default: None)")
    parser.add_argument("--registration_scale", type=str, default=DEFAULT_SPINE_REGISTRATION_SCALE,
                        help="Manual affine scale applied to the reference body before ICP. Use a scalar or sx,sy,sz. If omitted, PCA-based scaling is used.")
    parser.add_argument("--registration_min_scale", type=str, default=DEFAULT_SPINE_REGISTRATION_MIN_SCALE,
                        help="Minimum sx,sy,sz clamp for automatic PCA registration scaling. (default: %(default)s)")
    parser.add_argument("--registration_max_scale", type=str, default=DEFAULT_SPINE_REGISTRATION_MAX_SCALE,
                        help="Maximum sx,sy,sz clamp for automatic PCA registration scaling. (default: %(default)s)")
    parser.add_argument("--registration_backend", choices=("vtk",), default=DEFAULT_SPINE_REGISTRATION_BACKEND,
                        help="Spine reference alignment backend. (default: %(default)s)")
    parser.add_argument("--top_node_set_id", type=int, default=DEFAULT_SPINE_TOP_NODE_SET_ID,
                        help="ID for the top node set. (default: %(default)s)")
    parser.add_argument("--bottom_node_set_id", type=int, default=DEFAULT_SPINE_BOTTOM_NODE_SET_ID,
                        help="ID for the bottom node set. (default: %(default)s)")
    parser.add_argument("--pmma_yield_compression", type=float, default=None,
                        help="Sets the yield strength in compression for PMMA material in the FE model. (default: %(default)s [MPa])")
    parser.add_argument("--pmma_yield_tension", type=float, default=None,
                        help="Sets the yield strength in tension for PMMA material in the FE model. (default: %(default)s [MPa])")
#    parser.add_argument("--bone_yield_compression", type=float, default=None,
#                        help="Sets the yield strength in compression for bone material in the FE model. (default: %(default)s [MPa])")
#    parser.add_argument("--bone_yield_tension", type=float, default=None,
#                        help="Sets the yield strength in tension for bone material in the FE model. (default: %(default)s [MPa])")
#    parser.add_argument("--cort_elastic_Emax", type=float,
#                        help="Sets the maximum elastic modulus for cortical bone in the FE model. (default: %(default)s)")
#    parser.add_argument("--cort_elastic_exponent", type=float,
#                        help="Sets the elastic exponent for cortical bone in the FE model.")
    parser.add_argument("--cort_poissons_ratio", type=float,
                        help="Sets the Poisson's ratio for cortical bone in the FE model.")
#    parser.add_argument("--cort_yield_compression", type=float,
#                        help="Sets the yield strength in compression for cortical bone in the FE model.")
#    parser.add_argument("--cort_yield_tension", type=float,
#                        help="Sets the yield strength in tension for cortical bone in the FE model.")
    parser.add_argument("--elastic_E_func", type=str, default="default_E",
        help="Function name for trabecular bone Young’s modulus. (default: default_E)")
    parser.add_argument("--yield_comp_func", type=str,
        help="Function name for trabecular compression yield. (default: default_yc)")
    parser.add_argument("--yield_tens_func", type=str,
        help="Function name for trabecular tension yield. (default: default_yt)")

    parser.add_argument("--cort_elastic_E_func", type=str,
        help="Function name for cortical bone Young’s modulus. (default: default_E)")
    parser.add_argument("--cort_yield_comp_func", type=str,
        help="Function name for cortical compression yield. (default: default_yc)")
    parser.add_argument("--cort_yield_tens_func", type=str,
        help="Function name for cortical tension yield. (default: default_yt)")
    parser.add_argument("--appendix", default=None, type=str,
        help="Extra appendix for output file")
    parser.add_argument("--export_nifti", action='store_true',
        help="If set, exports resampled grayscale image and labeled segmentation.")

    args = parser.parse_args()

    # Print arguments
    ogo.message(echo_arguments("ogoFEA-spine-builder", vars(args)))

    # Prepare the output file
    basename = remove_extension(os.path.basename(args.calibrated_image))

    if args.output_path is None:
        output_dir = os.path.dirname(args.calibrated_image)
    else:
        output_dir = args.output_path

    output_dir = os.path.abspath(output_dir)

    if args.appendix is None:
        output_file = os.path.join(output_dir, f"{basename}_{args.mask_threshold}.n88model")
    else:
        output_file = os.path.join(output_dir, f"{basename}_{args.appendix}.n88model")

    ogo.message(f'N88Model File path: {output_file}')
    # Set default reference path if not provided
    reference_path = args.reference_path
    if reference_path is None:
        reference_path = str(default_spine_reference_path())

    # Extract all kwargs
    kwargs = vars(args).copy()  # Convert parsed arguments to a dictionary
    kwargs.pop("calibrated_image")  # Remove positional argument
    kwargs.pop("bone_mask")  # Remove positional argument
    kwargs.pop("mask_threshold")  # Remove required argument
    kwargs.pop("process_mask_threshold")  # Remove required argument
    kwargs.pop("reference_path")  # Ensure updated reference path

    # Run the vertebra processing
    process_vertebra(
        args.bone_mask,
        args.calibrated_image,
        output_file,
        args.mask_threshold,
        args.process_mask_threshold,
        reference_path,
        **kwargs  # Pass remaining arguments as kwargs
    )

###################################################################### MAIN
if __name__ == '__main__':
    main()
