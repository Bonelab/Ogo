#####
# SpineCompressionFE.py
#
# This script replaces the L4 vertebral compression FE model. It sets up a FAIM micro-FE model 
# from the density (K2HPO4) calibrated image. 
# This script sets up the model for a L4 vertebra (including arch and
# pedicles). The analysis resamples the image to isotropic voxels, transforms
# the image, applies the bone mask and bins the data. It then creates the FE model for
# solving using FAIM (>v8.0, Numerics Solutions Ltd, Calgary, Canada - Steven  Boyd).
#
# Changes to the original script:
# - Updated to include pedicles and body int he FE model
# - Requires the input of the body and process labels (does not identify body - use nnUnet-Model
# - Generates Jelly Bean disks that fit the surface of the vertebra perfectly 

# To-Do: 
# - Include option for no-disks (direct boundary conditions)
# - Include option for square and round disks

#####
#
# Matthias Walle
# Postdoctoral Associate
# University of Calgary
# Jan 24, 2025x
# Modifed to Py3: March 25, 2020
#####


import os
import argparse
import os
from glob import glob
import ogo.cli.ref.material_laws 

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import SimpleITK as sitk
from matplotlib import pyplot as plt
from scipy.ndimage import (
    binary_dilation,
    binary_erosion,
    binary_fill_holes,
    find_objects,
    label,
)
from skimage.exposure import rescale_intensity
import vtk
import vtkbone
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

import ogo.util.Helper as ogo
from ogo.util.echo_arguments import echo_arguments
from scipy.ndimage import gaussian_filter


vtk.vtkObject.GlobalWarningDisplayOff()

###################################################################### HELPERS (Python)
## General helper functions: 
###################################################################### HELPERS (Python)

def get_bounding_box(binary_mask):
    """
    Get the bounding box of a binary mask.
    Returns a tuple of slices defining the bounding box.
    """
    coords = np.array(np.where(binary_mask))
    min_coords = coords.min(axis=1)
    max_coords = coords.max(axis=1) + 1
    return tuple(slice(min_coords[i], max_coords[i]) for i in range(binary_mask.ndim))

def getLargestCC(segmentation):
    labels, _ = label(segmentation)
    assert( labels.max() != 0 ) # assume at least 1 CC
    largestCC = labels == np.argmax(np.bincount(labels.flat)[1:])+1
    return largestCC

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


###################################################################### HELPERS (VTK)
## VTK Helper functions:
###################################################################### HELPERS (VTK)


def read(input_mask):
    mask_reader = vtk.vtkNIFTIImageReader()
    mask_reader.SetFileName(input_mask)
    mask_reader.Update()
    return mask_reader

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


def merge_vtk_images(image_list, label_list=None):
    """
    Merges an arbitrary list of vtkImageData objects into a single vtkImageData.
    Optionally assigns specific label values for each image.

    Parameters:
    - image_list: List of vtkImageData objects to merge.
    - label_list: Optional list of label values corresponding to each image.
                  If None or an entry in the list is None, uses the original values in the image.

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
            merged_array[image_array > 0] = label_value
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

def crop_and_transform(fullvertebra, body_output, process_output, image):
  
    _, bb = crop_to_bounding_box(fullvertebra.GetOutput())
    isolated_vertebra, _ = crop_to_bounding_box(body_output, bb)
    isolated_process, _ = crop_to_bounding_box(process_output, bb)
    isolated_image, _ = crop_to_bounding_box(image, bb)
    
    
    return isolated_vertebra, isolated_process, isolated_image

def perform_marching_cubes(body_output):
    mcubes = vtk.vtkImageMarchingCubes()
    mcubes.SetInputConnection(body_output.GetOutputPort())
    mcubes.SetValue(1, 1.0)
    mcubes.Update()
    if mcubes.GetOutput().GetNumberOfPoints() == 0:
        raise ValueError("No data was generated by the Marching Cubes algorithm.")
    return mcubes

def get_icp_with_scaling(body, reference_path):

    def get_principal_axes_lengths(polydata):
        points = vtk_to_numpy(polydata.GetPoints().GetData())
        centered = points - np.mean(points, axis=0)
        cov = np.cov(centered.T)
        eigvals, eigvecs = np.linalg.eigh(cov)
        lengths = np.sqrt(eigvals) * 2  # Approx. diameter along each PCA axis
        return lengths

    # Get sample surface from marching cubes
    mcubes = perform_marching_cubes(body)
    mcubes_output = mcubes.GetOutput()

    # Load reference polydata
    reference_reader = vtk.vtkPolyDataReader()
    reference_reader.SetFileName(reference_path)
    reference_reader.Update()
    reference_polydata = reference_reader.GetOutput()

    # Compute PCA-based axis lengths
    sample_lengths = get_principal_axes_lengths(mcubes_output)
    ref_lengths = get_principal_axes_lengths(reference_polydata)
    scale_factors = sample_lengths / ref_lengths

    # Limit maximum scale factors
    min_scale = np.array([0.8, 0.8, 0.75])
    max_scale = np.array([1.2, 1.2, 1.3])
    scale_factors = np.minimum(scale_factors, max_scale)
    scale_factors = np.maximum(scale_factors, min_scale)

    ogo.message(f"PCA axis lengths (sample): {np.round(sample_lengths, 2)}")
    ogo.message(f"PCA axis lengths (reference): {np.round(ref_lengths, 2)}")
    ogo.message(f"Scale factors applied to reference: {np.round(scale_factors, 3)}")

    # Scale the reference to match the sample
    scale_transform = vtk.vtkTransform()
    scale_transform.Scale(*scale_factors)

    transform_filter = vtk.vtkTransformPolyDataFilter()
    transform_filter.SetInputData(reference_polydata)
    transform_filter.SetTransform(scale_transform)
    transform_filter.Update()
    scaled_reference = transform_filter.GetOutput()

    # ICP with scaled reference
    icp = vtk.vtkIterativeClosestPointTransform()
    icp.SetTarget(mcubes_output)
    icp.SetSource(scaled_reference)
    icp.StartByMatchingCentroidsOn()
    icp.GetLandmarkTransform().SetModeToRigidBody()
    icp.SetMeanDistanceModeToRMS()
    icp.SetMaximumMeanDistance(0.05)
    icp.CheckMeanDistanceOn()
    icp.SetMaximumNumberOfLandmarks(250)
    icp.SetMaximumNumberOfIterations(75)
    icp.Update()

    ogo.message("ICP (with scaled reference) Matrix:")
    print_matrix(icp.GetMatrix())

    return icp

def get_icp(body, reference_path):

    mcubes = perform_marching_cubes(body)
    mcubes_output = mcubes.GetOutput()
    reference_bone_reader = vtk.vtkPolyDataReader()
    reference_bone_reader.SetFileName(reference_path)
    reference_bone_reader.Update()
    
    icp = vtk.vtkIterativeClosestPointTransform()
    icp.SetTarget(mcubes_output)
    icp.SetSource(reference_bone_reader.GetOutput())
    
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

def transform_resample(image, matrix, iso_resolution): 
    transform = vtk.vtkTransform()
    transform.SetMatrix(matrix)
    transform.Update()
    
    reslice = vtk.vtkImageReslice()
    reslice.SetInputData(image)
    reslice.SetInterpolationModeToCubic()
    reslice.SetResliceTransform(transform)
    reslice.AutoCropOutputOn()
    reslice.Update()
    
    # After reslice
    if reslice.GetOutput().GetNumberOfPoints() == 0:
        ogo.message("Reslice output contains no points. Check the input and transformation.")
    
    image_resample = vtk.vtkImageResample()
    image_resample.SetInputData(reslice.GetOutput())
    image_resample.SetInterpolationModeToCubic()
    image_resample.SetDimensionality(3)
    image_resample.SetAxisOutputSpacing(0, iso_resolution)
    image_resample.SetAxisOutputSpacing(1, iso_resolution)
    image_resample.SetAxisOutputSpacing(2, iso_resolution)
    image_resample.Update()
    
    # After reslice
    if image_resample.GetOutput().GetNumberOfPoints() == 0:
        ogo.message("Resample output contains no points. Check the input and transformation.")
        
    return image_resample

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
    eroded_mask = binary_erosion(binary_mask, iterations=2)
    #dilated_mask = binary_dilation(binary_mask, iterations=1)
    outline_mask = (eroded_mask == 0) & (binary_mask == 1)

    # Highlight the outline
    struct_element_erosion = np.array([[[0, 1, 0], [1, 1, 1], [0, 1, 0]]])
    highlighted_mask = binary_erosion(outline_mask, structure=struct_element_erosion, iterations=2)

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



def extrude_disk(boundary_labelled_vtk_image, value, pmma_thick=5, direction='up', axis='x'):
    """
    Extend the surface along the specified axis by projecting the surface in the
    remaining two axes and filling the bounding box for each slice along the chosen axis.
    Values in the original image are masked out.

    Parameters:
    - boundary_labelled_vtk_image: VTK image with labeled boundary surface.
    - value: Integer value representing the labeled surface to extrude.
    - pmma_thick: Number of slices to extrude.
    - direction: Direction of extrusion ('up' or 'down').
    - axis: Axis along which to extrude ('x', 'y', 'z').

    Returns:
    - extruded_vtk_image: VTK image with the extruded and bounded surface.
    """
    # Convert VTK image data to NumPy array
    input_array = vtk_to_numpy(boundary_labelled_vtk_image.GetPointData().GetScalars()).reshape(
        boundary_labelled_vtk_image.GetDimensions(), order='F'
    )
    input_array = np.swapaxes(input_array, 0, 2)  # Swap axes for correct orientation

    # Create a binary mask for the labeled surface
    binary_mask = input_array == value

    # Determine bounding box of the binary mask
    bb = get_bounding_box(binary_mask)

    # Adjust the bounding box based on direction and pmma_thick
    axis_map = {'x': 0, 'y': 1, 'z': 2}
    extrusion_axis = axis_map[axis]
    bb_list = list(bb)  # Convert tuple to list for modification

    if direction == 'up':
        bb_list[extrusion_axis] = slice(
            bb[extrusion_axis].start, 
            min(bb[extrusion_axis].stop + pmma_thick, binary_mask.shape[extrusion_axis])
        )
    elif direction == 'down':
        bb_list[extrusion_axis] = slice(
            max(bb[extrusion_axis].start - pmma_thick, 0), 
            bb[extrusion_axis].stop
        )
    else:
        raise ValueError("Direction must be 'up' or 'down'.")

    # Apply the adjusted bounding box
    bb = tuple(bb_list)
    ogo.message(f"Adjusted bounding box: {bb}")

    # Crop and project the mask
    cropped_mask = binary_mask[bb]
    cross_section = cropped_mask.any(axis=extrusion_axis)

    # Fill each slice in the cropped mask with the projection
    for idx in range(cropped_mask.shape[extrusion_axis]):
        if extrusion_axis == 0:  # Extrusion along x
            cropped_mask[idx, :, :] = cross_section
        elif extrusion_axis == 1:  # Extrusion along y
            cropped_mask[:, idx, :] = cross_section
        elif extrusion_axis == 2:  # Extrusion along z
            cropped_mask[:, :, idx] = cross_section

    # fill holes
    pad_width = [(1, 1) if i == extrusion_axis else (0, 0) for i in range(3)]

    # Pad the binary mask with True (1) to simulate closed boundaries
    padded_mask = np.pad(cropped_mask, pad_width=pad_width, mode='constant', constant_values=1)
    # Fill holes in the padded mask
    filled_padded = binary_fill_holes(padded_mask)
    # Remove the padding to return to original size
    cropped_mask = filled_padded[tuple(slice(1, -1) if i == extrusion_axis else slice(None) for i in range(3))]
    
    # Create the output array and mask out original values
    output_array = np.zeros_like(input_array, dtype=np.uint16)
    output_array[bb] = cropped_mask
    output_array[input_array != 0] = False
    largestcc = getLargestCC(output_array)
    # Convert back to VTK format
    final_array = np.swapaxes(largestcc, 0, 2)

    # Adjust VTK metadata
    origin = list(boundary_labelled_vtk_image.GetOrigin())
    spacing = boundary_labelled_vtk_image.GetSpacing()
    dimensions = list(final_array.shape)

    vtk_array = numpy_to_vtk(final_array.ravel(order='F').astype(np.uint16), deep=True)


    # Convert back to VTK format
    extruded_vtk_image = boundary_labelled_vtk_image.NewInstance()
    extruded_vtk_image.SetDimensions(dimensions)
    extruded_vtk_image.SetOrigin(origin)
    extruded_vtk_image.SetSpacing(spacing)
    extruded_vtk_image.GetPointData().SetScalars(vtk_array)

    return extruded_vtk_image

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
    change = ogo.changeInfo(image)
    connected_image = ogo.imageConnectivity(change)
    thr_image = ogo.bmd_preprocess(connected_image, -31)
    #this was done previously - but if someone wants this they should just do it in the calib functions
    #ash_image = ogo.bmd_K2hpo4ToAsh(thr_image)
    cast_image = ogo.cast2short(thr_image)
    bone_image = ogo.applyMask(cast_image, mask)
    binned_image, bin_centers = ogo.density2materialID(bone_image, n_bins=n_bins, cort_mask=cort_mask)

    return binned_image, bin_centers

def find_and_add_visible_nodes(model, bc_geometry, normal_vector, bone_material_id, node_set_name):
    visibleNodesIds = vtk.vtkIdTypeArray()
    vtkbone.vtkboneNodeSetsByGeometry.FindNodesOnVisibleSurface(
        visibleNodesIds, bc_geometry, normal_vector, bone_material_id)
    visibleNodesIds.SetName(node_set_name)
    model.AddNodeSet(visibleNodesIds)
    ogo.message(f"Found {visibleNodesIds.GetNumberOfTuples()} visible nodes for {node_set_name}.")
    return model 

def resolve_func(func_or_name, module):
    if callable(func_or_name):
        return func_or_name
    elif isinstance(func_or_name, str):
        return getattr(module, func_or_name)
    else:
        return None

def apply_boundary_conditions(model,**kwargs):

    fe_displacement = kwargs.get("fe_displacement", -1.0)

    ogo.message('Applying boundary conditions...')
    model.ApplyBoundaryCondition("body_top", vtkbone.vtkboneConstraint.SENSE_Z, fe_displacement, "top_displacement")
    model.ApplyBoundaryCondition("body_bottom", vtkbone.vtkboneConstraint.SENSE_Z, 0, "bottom_fixed_z")
    return model


def log_fe_arguments(**kwargs):
    import inspect

    ogo.message("========= FE MODEL ARGUMENTS =========")
    for key, value in kwargs.items():
        if inspect.isfunction(value):
            ogo.message(f"{key}: function -> {value.__name__}")
        else:
            ogo.message(f"{key}: {value}")
    ogo.message("======================================")

def create_microfe_model(
    image_with_pads,
    boundary_masks_with_pads,
    bin_centers,
    **kwargs
):
    import ogo.cli.ref.material_laws as material_laws

    # === Load defaults and dynamic functions === #
    n_bins = 128
    poissons_ratio = kwargs.get("poissons_ratio", 0.3)
    pmma_mat_id = kwargs.get("pmma_mat_id", 5000)
    pmma_E = kwargs.get("pmma_E", 2500)
    pmma_v = kwargs.get("pmma_v", 0.3)
    top_displacement = kwargs.get("top_displacement", "top_displacement")
    top_direction = kwargs.get("top_direction", (0, 0, 1))
    bottom_direction = kwargs.get("bottom_direction", (0, 0, -1))
    top_node_set_id = kwargs.get("top_node_set_id", 4)
    bottom_node_set_id = kwargs.get("bottom_node_set_id", 3)
    top_node_set_name = kwargs.get("top_node_set_name", "body_top")
    bottom_node_set_name = kwargs.get("bottom_node_set_name", "body_bottom")
    pmma_yield_compression = kwargs.get("pmma_yield_compression", None)
    pmma_yield_tension = kwargs.get("pmma_yield_tension", None)

    # === Load material law functions from kwargs === #
    elastic_E_func =  resolve_func(kwargs.get("elastic_E_func"), material_laws) or default_E
    yield_comp_func = resolve_func(kwargs.get("yield_comp_func"), material_laws)
    yield_tens_func = resolve_func(kwargs.get("yield_tens_func"), material_laws)

    cort_elastic_E_func = resolve_func(kwargs.get("cort_elastic_E_func"), material_laws) or elastic_E_func
    cort_yield_comp_func = resolve_func(kwargs.get("cort_yield_comp_func"), material_laws) or yield_comp_func
    cort_yield_tens_func = resolve_func(kwargs.get("cort_yield_tens_func"), material_laws) or yield_tens_func
    cort_poissons_ratio = poissons_ratio if kwargs.get('cort_poissons_ratio') is None else kwargs['cort_poissons_ratio']

    log_fe_arguments(
        n_bins=n_bins,
        poissons_ratio=poissons_ratio,
        pmma_mat_id=pmma_mat_id,
        pmma_E=pmma_E,
        pmma_v=pmma_v,
        top_displacement=top_displacement,
        top_direction=top_direction,
        bottom_direction=bottom_direction,
        top_node_set_id=top_node_set_id,
        bottom_node_set_id=bottom_node_set_id,
        top_node_set_name=top_node_set_name,
        bottom_node_set_name=bottom_node_set_name,
        pmma_yield_compression=pmma_yield_compression,
        pmma_yield_tension=pmma_yield_tension,
        elastic_E_func=elastic_E_func,
        yield_comp_func=yield_comp_func,
        yield_tens_func=yield_tens_func,
        cort_elastic_E_func=cort_elastic_E_func,
        cort_yield_comp_func=cort_yield_comp_func,
        cort_yield_tens_func=cort_yield_tens_func,
        cort_poissons_ratio=cort_poissons_ratio
    )

    # === Mesh images === #
    ogo.message(f"Casting to Short Integer datatype...")
    image_pads_short = ogo.cast2short(image_with_pads)

    ogo.message(f"Filtering connected components...")
    conn = ogo.imageConnectivity(image_pads_short)
    conn_bc = ogo.imageConnectivity(boundary_masks_with_pads)

    ogo.message(f"Meshing...")
    mesh = ogo.Image2Mesh(conn)
    temp_bc_mesh = ogo.Image2Mesh(conn_bc)

    # === Build material table === #
    ogo.message("Setting up the Finite Element Material Table...")
    material_table = vtkbone.vtkboneMaterialTable()

    # Trabecular material
    material_table = ogo.add_bone_material(
        material_table, bin_centers,
        elastic_E_func=elastic_E_func,
        mu=poissons_ratio,
        yield_comp_func=yield_comp_func,
        yield_tens_func=yield_tens_func,
        bin_range=(0, n_bins),
        material_name='TrabBone'
    )

    # Cortical material
    material_table = ogo.add_bone_material(
        material_table, bin_centers,
        elastic_E_func=cort_elastic_E_func,
        mu=cort_poissons_ratio,
        yield_comp_func=cort_yield_comp_func,
        yield_tens_func=cort_yield_tens_func,
        bin_range=(n_bins, 2 * n_bins),
        material_name='CortBone'
    )

    # PMMA caps
    material_table = ogo.add_pmma_material(
        material_table,
        pmma_mat_id,
        pmma_E,
        pmma_v,
        pmma_yield_tension=pmma_yield_tension,
        pmma_yield_compression=pmma_yield_compression
    )

    # === Build and annotate model === #
    ogo.message("Constructing the Finite Element Model...")
    model = ogo.applyTestBase(mesh, material_table)
    model.ComputeBounds()

    ogo.message(f"Identifying boundary nodes...")
    model = find_and_add_visible_nodes(model, temp_bc_mesh, top_direction, top_node_set_id, top_node_set_name)
    model = find_and_add_visible_nodes(model, temp_bc_mesh, bottom_direction, bottom_node_set_id, bottom_node_set_name)
    model = apply_boundary_conditions(model, **kwargs)

    ogo.message(f"Setting convergence criteria...")
    model.ConvergenceSetFromConstraint(top_displacement)

    # === Postprocessing sets === #
    ogo.message('Postprocessing...')
    info = model.GetInformation()
    pp_node_sets_key = vtkbone.vtkboneSolverParameters.POST_PROCESSING_NODE_SETS()
    pp_elem_sets_key = vtkbone.vtkboneSolverParameters.POST_PROCESSING_ELEMENT_SETS()
    for setname in [top_node_set_name, bottom_node_set_name]:
        pp_node_sets_key.Append(info, setname)
        elementSet = model.GetAssociatedElementsFromNodeSet(setname)
        model.AddElementSet(elementSet)
        pp_elem_sets_key.Append(info, setname)

    return model


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
    base_name = os.path.basename(filepath)
    parts = base_name.split('_')
    return {'ID': parts[0], 'TREATMENT': parts[1], 'LOCATION': parts[2], 'NUMBER': parts[4], 'filename':'_'.join(parts[:-3])}

def visualize_qc_figure(image_vtk, label_vtk, filepath, title=None, body_label=2, offset=10):
    """
    Create a QC figure with:
    - 3 orthogonal views (sagittal, coronal, axial)
    - slices centered on the vertebral body centroid
    - grayscale image background
    - colored transparent label overlay
    - grid + crosshairs
    - 3 stacked 3D views with different azimuths, keeping axial direction as z-axis
    - color legend using actual swatches

    Parameters
    ----------
    image_vtk : vtk.vtkImageData
        Grayscale image.
    label_vtk : vtk.vtkImageData
        Integer label image for overlay / 3D rendering.
    filepath : str
        Output PNG path.
    title : str, optional
        Figure title.
    body_label : int, optional
        Label used to define vertebral body center. Default = 2.
    offset : int, optional
        Slice offsets around the body center in voxels. Default = 10.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from vtk.util.numpy_support import vtk_to_numpy
    from skimage.measure import marching_cubes

    def vtk_to_numpy_zyx(vtk_img):
        dims = vtk_img.GetDimensions()              # (x, y, z)
        arr = vtk_to_numpy(vtk_img.GetPointData().GetScalars())
        arr = arr.reshape(dims, order="F")          # (x, y, z)
        arr = np.transpose(arr, (2, 1, 0))          # -> (z, y, x)
        return arr

    def robust_window(img, lower=1, upper=99):
        vals = img[np.isfinite(img)]
        if vals.size == 0:
            return img, 0.0, 1.0
        vmin, vmax = np.percentile(vals, [lower, upper])
        if vmax <= vmin:
            vmax = vmin + 1.0
        img_clip = np.clip(img, vmin, vmax)
        return img_clip, vmin, vmax

    def extract_slice(vol, axis, idx):
        if axis == 0:      # axial: z slice
            return vol[idx, :, :]
        elif axis == 1:    # coronal: y slice
            return vol[:, idx, :]
        elif axis == 2:    # sagittal: x slice
            return vol[:, :, idx]
        raise ValueError("axis must be 0, 1, or 2")

    def orient_for_display(img2d):
        return np.flipud(img2d)

    def clamp(i, maxv):
        return max(0, min(int(i), maxv - 1))

    def get_body_center(labels_arr, target_label=2):
        coords = np.where(labels_arr == target_label)
        if coords[0].size == 0:
            return (
                labels_arr.shape[0] // 2,
                labels_arr.shape[1] // 2,
                labels_arr.shape[2] // 2,
            )
        zc = int(np.mean(coords[0]))
        yc = int(np.mean(coords[1]))
        xc = int(np.mean(coords[2]))
        return zc, yc, xc

    def get_bbox(labels_arr, pad=6):
        coords = np.where(labels_arr > 0)
        if coords[0].size == 0:
            zdim, ydim, xdim = labels_arr.shape
            return (0, zdim, 0, ydim, 0, xdim)

        zmin = max(0, coords[0].min() - pad)
        zmax = min(labels_arr.shape[0], coords[0].max() + pad + 1)
        ymin = max(0, coords[1].min() - pad)
        ymax = min(labels_arr.shape[1], coords[1].max() + pad + 1)
        xmin = max(0, coords[2].min() - pad)
        xmax = min(labels_arr.shape[2], coords[2].max() + pad + 1)
        return zmin, zmax, ymin, ymax, xmin, xmax

    def crop_2d_for_axis(img2d, lab2d, axis_name, bbox):
        zmin, zmax, ymin, ymax, xmin, xmax = bbox

        if axis_name == "Sagittal":
            # shape: (z, y)
            return img2d[zmin:zmax, ymin:ymax], lab2d[zmin:zmax, ymin:ymax]
        elif axis_name == "Coronal":
            # shape: (z, x)
            return img2d[zmin:zmax, xmin:xmax], lab2d[zmin:zmax, xmin:xmax]
        elif axis_name == "Axial":
            # shape: (y, x)
            return img2d[ymin:ymax, xmin:xmax], lab2d[ymin:ymax, xmin:xmax]
        return img2d, lab2d

    def add_surface(ax, mask_xyz, color, alpha):
        if np.any(mask_xyz):
            try:
                verts, faces, _, _ = marching_cubes(mask_xyz.astype(np.uint8), level=0.5)
                mesh = Poly3DCollection(verts[faces], alpha=alpha)
                mesh.set_facecolor(color)
                mesh.set_edgecolor("none")
                ax.add_collection3d(mesh)
            except Exception:
                pass

    # Convert inputs
    image = vtk_to_numpy_zyx(image_vtk).astype(np.float32)
    labels = vtk_to_numpy_zyx(label_vtk).astype(np.int16)

    # Windowing
    image, vmin, vmax = robust_window(image, 1, 99)

    # Dimensions
    zdim, ydim, xdim = image.shape

    # Center slices on vertebral body
    zc, yc, xc = get_body_center(labels, target_label=body_label)

    slice_positions = {
        "Sagittal": (
            2,
            [clamp(xc - offset, xdim), clamp(xc, xdim), clamp(xc + offset, xdim)],
        ),
        "Coronal": (
            1,
            [clamp(yc - offset, ydim), clamp(yc, ydim), clamp(yc + offset, ydim)],
        ),
        "Axial": (
            0,
            [clamp(zc - offset, zdim), clamp(zc, zdim), clamp(zc + offset, zdim)],
        ),
    }

    bbox = get_bbox(labels, pad=6)

    # Overlay colors
    overlay_colors = [
        (0, 0, 0, 0.0),            # 0 background transparent
        (0.90, 0.35, 0.35, 0.55),  # 1 process
        (0.30, 0.60, 0.95, 0.45),  # 2 body
        (0.20, 0.80, 0.30, 0.60),  # 3 inferior BC
        (0.95, 0.75, 0.15, 0.60),  # 4 superior BC
        (0.70, 0.30, 0.90, 0.60),  # 5 inferior disk (optional)
        (0.15, 0.85, 0.85, 0.60),  # 6 superior disk (optional)
    ]
    cmap = ListedColormap(overlay_colors)
    norm = BoundaryNorm(np.arange(-0.5, len(overlay_colors) + 0.5, 1), cmap.N)

    # For legends / 3D surfaces
    label_specs = [
        (1, "Process", (0.90, 0.35, 0.35), 0.45),
        (2, "Body", (0.30, 0.60, 0.95), 0.38),
        (3, "Inferior BC", (0.20, 0.80, 0.30), 0.60),
        (4, "Superior BC", (0.95, 0.75, 0.15), 0.60),
        (5, "Inferior disk", (0.70, 0.30, 0.90), 0.50),
        (6, "Superior disk", (0.15, 0.85, 0.85), 0.50),
    ]

    # Figure layout: 3 rows x 4 cols
    fig = plt.figure(figsize=(20, 14), constrained_layout=True)
    gs = fig.add_gridspec(3, 4, width_ratios=[1, 1, 1, 1.15])

    if title is not None:
        fig.suptitle(title, fontsize=18, fontweight="bold")

    # 2D panels
    row = 0
    for view_name, (axis, idxs) in slice_positions.items():
        for col, idx in enumerate(idxs):
            ax = fig.add_subplot(gs[row, col])

            img2d = extract_slice(image, axis, idx)
            lab2d = extract_slice(labels, axis, idx)

            img2d, lab2d = crop_2d_for_axis(img2d, lab2d, view_name, bbox)

            img2d = orient_for_display(img2d)
            lab2d = orient_for_display(lab2d)

            ax.imshow(img2d, cmap="gray", vmin=vmin, vmax=vmax, interpolation="nearest")
            ax.imshow(lab2d, cmap=cmap, norm=norm, interpolation="nearest")

            h, w = img2d.shape
            ax.axhline(h // 2, color="white", lw=0.8, alpha=0.7)
            ax.axvline(w // 2, color="white", lw=0.8, alpha=0.7)

            ax.set_xticks(np.linspace(0, max(w - 1, 1), 6))
            ax.set_yticks(np.linspace(0, max(h - 1, 1), 6))
            ax.grid(color="white", linestyle=":", linewidth=0.5, alpha=0.35)

            ax.set_title(f"{view_name} @ {idx}", fontsize=12, fontweight="bold")
            ax.set_xlabel("pixels")
            ax.set_ylabel("pixels")

        row += 1

    # 3D panels
    # Convert (z, y, x) -> (x, y, z), so axial stays as vertical z-axis in the rendered plot
    labels_xyz = np.transpose(labels, (2, 1, 0))
    x3, y3, z3 = labels_xyz.shape

    azimuths = [135, 45, -45]
    view_titles = ["3D oblique", "3D lateral", "3D opposite oblique"]

    for i in range(3):
        ax3d = fig.add_subplot(gs[i, 3], projection="3d")
        ax3d.set_title(view_titles[i], fontsize=12, fontweight="bold")

        # translucent overall shell
        add_surface(ax3d, labels_xyz > 0, color=(0.82, 0.84, 0.90), alpha=0.10)

        # label surfaces
        for label_val, _, color, alpha in label_specs:
            add_surface(ax3d, labels_xyz == label_val, color=color, alpha=alpha)

        ax3d.set_xlim(0, x3)
        ax3d.set_ylim(0, y3)
        ax3d.set_zlim(0, z3)
        ax3d.set_box_aspect((x3, y3, z3))

        # Same elevation, only rotate around axial axis
        ax3d.view_init(elev=22, azim=azimuths[i])

        ax3d.set_xlabel("X")
        ax3d.set_ylabel("Y")
        ax3d.set_zlabel("Z")

        # cleaner panes
        try:
            ax3d.xaxis.pane.fill = False
            ax3d.yaxis.pane.fill = False
            ax3d.zaxis.pane.fill = False
        except Exception:
            pass

        # put color legend only on bottom 3D panel
        if i == 2:
            handles = [
                Patch(facecolor=color, edgecolor="none", label=name)
                for _, name, color, _ in label_specs
                if np.any(labels == _[0] if False else labels == _)
            ]
            # safer explicit version
            handles = []
            for label_val, name, color, _alpha in label_specs:
                if np.any(labels == label_val):
                    handles.append(Patch(facecolor=color, edgecolor="none", label=name))

            if handles:
                ax3d.legend(
                    handles=handles,
                    loc="upper left",
                    bbox_to_anchor=(0.0, 0.0),
                    frameon=False,
                    fontsize=9,
                    title="Overlay colors",
                )

    plt.savefig(filepath, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    ogo.message(f"QC figure saved to {filepath}")

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
    plt.figure(figsize=(10, 10))
    plt.imshow(normalized_slice, cmap='viridis', origin="lower")
    plt.title('Middle X Slice (Clipped and Normalized)')
    plt.xlabel('Y axis')
    plt.ylabel('Z axis')

    # Save the output image
    plt.savefig(filepath)
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
    iso_resolution = kwargs.get("iso_resolution", 1.0)
    pmma_thick = kwargs.get("pmma_thick", 5)
    top_node_set_id = kwargs.get("top_node_set_id", 4)
    bottom_node_set_id = kwargs.get("bottom_node_set_id", 3)
    quality_control = kwargs.get("quality_control", True)
   
    #Read Image and Mask
    ogo.message(f"Reading image...: {input_image}")
    image_reader = read(input_image)
    ogo.message(f"Reading mask...: {input_mask}")
    mask_reader = resample_to_match(image_reader.GetOutput(), read(input_mask).GetOutput())

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
    isolated_vertebra, isolated_process, isolated_image = crop_and_transform(fullvertebra, body.GetOutput(), process.GetOutput(), image_reader.GetOutput())

    # Marching Cubes and Registration
    ogo.message(f"Starting ICP registration to refrence...: {reference_path} ")
    icp = get_icp_with_scaling(body, reference_path)
    
    # Transform Images and resample at the same time (less interpolation)
    transformed_vertebra = transform_resample(isolated_vertebra.GetOutput(), icp.GetMatrix(), iso_resolution)
    transformed_process = transform_resample(isolated_process.GetOutput(), icp.GetMatrix(), iso_resolution)
    transformed_image = transform_resample(isolated_image.GetOutput(), icp.GetMatrix(), iso_resolution)

    ogo.message(f"relabelling mask and identifying boundary surfaces...")
    # Assuming mask1_data and mask2_data are your initial binary masks
    mask1_labeled = label_mask(transformed_vertebra.GetOutput(), 2)  # Apply label 1 to the first mask
    mask2_labeled = label_mask(transformed_process.GetOutput(), 1)  # Apply label 2 to the second mask
    
    # Combine masks here (otherwise it does weird interpolations)
    combined_masks = add_masks(mask1_labeled.GetOutput(), mask2_labeled.GetOutput())  # Combine the labeled masks
    
    # This creates an image identifying the top and bottom surfaces of the vertebra
    boundary_masks = identify_boundary_surface(combined_masks.GetOutput(), bottom_node_set_id, top_node_set_id)
    
    ogo.message("creating cortical mask....")
    #1 - 3 alterantive (to create a continous cortical shell - changed to 0 5 )

    cort_mask = generate_cortical_mask(transformed_image.GetOutput(), combined_masks.GetOutput(), threshold=1, min_th=2, max_th=3)
    
    #writer = vtk.vtkXMLImageDataWriter()
    #writer.SetFileName('/home/matthias.walle/work/fem/CT_FE_TEMPLATE/MODELS/10001_QCT_vertebra_20_cortmask.vti')
    #writer.SetInputData(cort_mask)  # For VTK 8+
    #writer.Write()

    ogo.message(f"converting image to material ID...")
    # This converts the raw image to a material ID mapped image
    n_bins = 128
    bone_image, bin_centers = convert_image_to_material(transformed_image.GetOutput(), combined_masks.GetOutput(), n_bins=n_bins, cort_mask = cort_mask)
    
    ogo.message(f"padding images...")
    # Pad images before we add the disks (otherwise they may not have right pmma_thick)
    padded_image = pad_vtk_image(bone_image, axis='x', pmma_thick=pmma_thick, pad_value=0)
    padded_transformed_image = pad_vtk_image(transformed_image.GetOutput(), axis='x', pmma_thick=pmma_thick, pad_value=0)
    padded_mask = pad_vtk_image(boundary_masks, axis='x', pmma_thick=pmma_thick, pad_value=0)
    padded_cort_mask = pad_vtk_image(cort_mask, axis='x', pmma_thick=pmma_thick, pad_value=0)
    padded_mask1 = pad_vtk_image(mask1_labeled.GetOutput(), axis='x', pmma_thick=pmma_thick, pad_value=0)
    padded_mask2 = pad_vtk_image(mask2_labeled.GetOutput(), axis='x', pmma_thick=pmma_thick, pad_value=0)

    ogo.message(f"extruding boundary diks (pmma_thick > {pmma_thick} mm)...")
    # Generate disks for the top and bottom surfaces (this could be easily changed to have no disk)
    inferior_disk = extrude_disk(padded_mask, value=bottom_node_set_id, pmma_thick=pmma_thick, direction='down', axis='x')
    superior_disk = extrude_disk(padded_mask, value=top_node_set_id, pmma_thick=pmma_thick, direction='up', axis='x')

    ogo.message(f"merging disks with image...")
    # Now we create the image with the disks and the boundary image as well
    image_with_pads = merge_vtk_images([
        padded_image, superior_disk, inferior_disk], [None, pmma_mat_id,pmma_mat_id])

    boundary_masks_with_pads = merge_vtk_images([
        padded_mask, inferior_disk, superior_disk], [None, bottom_node_set_id, top_node_set_id])

    # This is the final model that we will use for the FEA. 
    # For the future --> pull out material table to have option without disks
    ogo.message("generating n88model file...")
    model = create_microfe_model(image_with_pads, boundary_masks_with_pads, bin_centers, **kwargs)

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
            visualize_qc_figure(
                image_vtk=padded_transformed_image,
                label_vtk=boundary_masks_with_pads,
                filepath=image_path
               #title=os.path.basename(n88model_output_path)
            )

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
        prog="OgoSpineCompressionFe",
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
    parser.add_argument("--quality_control", type=bool, default=True,  
                        help="Set quality control flag (visualise output and add volume checks). (default: %(default)s)")
    parser.add_argument("--iso_resolution", type=float, default=1.0,
                        help="Set the isotropic voxel size [in mm]. (default: %(default)s [mm])")
    #parser.add_argument("--elastic_exponent", type=float, default=2.29,
    #                    help="Sets the exponent (b) for power law: E=A(den)^b. (default: %(default)s)")
    #parser.add_argument("--elastic_Emax", type=float, default=10500,
    #                    help="Sets the coefficient (A) Elastic Modulus value for the power law: E=A(den)^b. (default: %(default)s [MPa])")
    parser.add_argument("--poissons_ratio", type=float, default=0.3,
                        help="Sets the Poisson's ratio for the material(s) in the FE model. (default: %(default)s)")
    parser.add_argument("--pmma_E", type=float, default=2500,
                        help="Sets the Elastic Modulus for PMMA caps in the FE model. (default: %(default)s [MPa])")
    parser.add_argument("--pmma_v", type=float, default=0.3,
                        help="Sets the Poisson's ratio for the PMMA material(s) in the FE model. (default: %(default)s)")
    parser.add_argument("--pmma_thick", type=int, default=3,
                        help="Sets the minimum pmma_thick for PMMA caps in the FE model. (default: %(default)s [mm])")
    parser.add_argument("--pmma_mat_id", type=int, default=5000,
                        help="Sets the material ID for the PMMA blocks. (default: %(default)s)")
    parser.add_argument("--fe_displacement", type=float, default=-1.0,
                        help="Sets the applied displacement in [mm] to the FE model. (default: %(default)s [mm])")
    parser.add_argument("--reference_path", type=str, required=False, default=None,
                        help="Path to the reference vtk file for ICP registration. (default: None)")
    parser.add_argument("--top_node_set_id", type=int, default=4,
                        help="ID for the top node set. (default: %(default)s)")
    parser.add_argument("--bottom_node_set_id", type=int, default=3,
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
    ogo.message(echo_arguments('OgoSpineCompressionFe', vars(args)))

    # Prepare the output file
    basename = remove_extension(os.path.basename(args.calibrated_image))
    
    if args.output_path is None:
        output_dir = os.path.dirname(args.calibrated_image)
    else:
        output_dir = args.output_path
    
    output_dir = os.path.abspath(output_dir)

    if args.appendix is None:
        output_file = os.path.join(output_dir, f"{basename}_vertebra_{args.mask_threshold}.n88model")
    else:
        output_file = os.path.join(output_dir, f"{basename}_vertebra_{args.mask_threshold}_{args.appendix}.n88model")

    ogo.message(f'N88Model File path: {output_file}')
    # Set default reference path if not provided
    reference_path = args.reference_path
    if reference_path is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        data_dir = os.path.join(script_dir, "../../dat")
        reference_path = os.path.join(data_dir, "L4_BODY_SPINE_COMPRESSION_REF.vtk")

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