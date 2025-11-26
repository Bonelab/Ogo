#!/usr/bin/env python

import SimpleITK as sitk
import numpy as np
import os
import argparse

def segment_air(input_file: str, output_mask_file: str, air_max_threshold: float = -990, air_min_threshold: float = -1010):
    """
    Segments air from a NIfTI image using a simple thresholding method.

    Parameters:
        input_file (str): Path to the input NIfTI file.
        output_mask_file (str): Path to save the output mask file.
        air_max_threshold (float): Intensity value below which voxels are considered air (default: -990).
        air_min_threshold (float): Intensity value above which voxels are consisdered air (default: -1010). 
    """
    
    print(f"Processing {input_file}...")
    image = sitk.ReadImage(input_file)

    image_array = sitk.GetArrayFromImage(image)
    air_mask = np.zeros_like(image_array)
    air_mask[(image_array >= air_min_threshold) & (image_array <= air_max_threshold)] = 1

    air_mask_image = sitk.GetImageFromArray(air_mask)
    component_image = sitk.ConnectedComponent(air_mask_image)
    sorted_component_image = sitk.RelabelComponent(component_image, sortByObjectSize=True)
    
    component1_mask = sorted_component_image == 1
    component2_mask = sorted_component_image == 2

    
    component1_array = sitk.GetArrayFromImage(component1_mask)
    component2_array = sitk.GetArrayFromImage(component2_mask)

    
    mean1 = image_array[component1_array == 1].mean()
    mean2 = image_array[component2_array == 1].mean()

    # Choose the one closer to -1000
    target_mean = -1000
    if abs(mean1 - target_mean) < abs(mean2 - target_mean):
        selected_component_mask = component1_mask
    else:
        selected_component_mask = component2_mask

    
    component1_mask.CopyInformation(image)
    sitk.WriteImage(selected_component_mask, output_mask_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Segment air from a NIfTI image using thresholding.")
    parser.add_argument("input_filename", help="Path to the input NIfTI file.")
    parser.add_argument("output_filename", help="Path to save the output air segmentation mask.")
    parser.add_argument("--air_min_threshold", type=float, default=-1010, help="Minimum intensity threshold for air segmentation (default: -1010).")
    parser.add_argument("--air_max_threshold", type=float, default=-990, help="Maximum intensity threshold for air segmentation (default: -990).")

    args = parser.parse_args()

    # Check if the input file exists
    if not os.path.exists(args.input_filename):
        print(f"Error: Input file does not exist: {args.input}")
    else:
        segment_air(input_file=args.input_filename, 
                    output_mask_file=args.output_filename, 
                    air_min_threshold=args.air_min_threshold, 
                    air_max_threshold=args.air_max_threshold)


