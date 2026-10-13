# /------------------------------------------------------------------------------+
# | 19-JUL-2022                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+

# Imports
import argparse
import os
import SimpleITK as sitk
import math
import numpy as np
from ogo.util.echo_arguments import echo_arguments
import ogo.util.Helper as ogo
import ogo.dat.OgoMasterLabels as lb
from ogo.util.cli_tools import (
    check_overwrite,
    validate_input_file,
    validate_output_file,
    validate_integer_range,
    load_config_from_json,
    save_config_to_json,
    print_image_info
)

# Pixel type value ranges (min, max)
PIXEL_TYPE_RANGES = {
    sitk.sitkInt8:    (-128, 127),
    sitk.sitkUInt8:   (0, 255),
    sitk.sitkInt16:   (-32768, 32767),
    sitk.sitkUInt16:  (0, 65535),
    sitk.sitkInt32:   (-2147483648, 2147483647),
    sitk.sitkUInt32:  (0, 4294967295),
}

def ImageThreshold(input_file, output_file, lower, upper, inside, outside, overwrite, quiet):
    
    # Check if output exists and should overwrite
    if not check_overwrite(output_file, overwrite):
        return 1
    
    # Validate threshold parameters
    if lower > upper:
        ogo.message('[ERROR] Lower threshold cannot be greater than upper threshold.')
        ogo.message(f'        lower={lower}, upper={upper}')
        return 1
    
    if inside == outside:
        if not quiet:
            ogo.message('[WARNING] Setting inside equal to outside makes no sense!')
    
    # Read input
    if not quiet:
        ogo.message(f'Reading input image:')
        ogo.message(f'  {input_file}')
    
    input_image = sitk.ReadImage(input_file)
    
    if not quiet:
        print_image_info(input_file, input_image)
    
    # Checking image type and validate thresholds
    image_type = input_image.GetPixelIDValue()
    if not quiet:
        ogo.message(f'Input image type is {input_image.GetPixelIDTypeAsString()}')
    
    # Get type min/max from constant dictionary
    if image_type not in PIXEL_TYPE_RANGES:
        ogo.message('[ERROR] Unexpected image input type.')
        return 1
    
    type_min, type_max = PIXEL_TYPE_RANGES[image_type]
    
    # Adjust thresholds if they exceed type limits
    lower_adjusted = lower
    upper_adjusted = upper
    
    if lower < type_min:
        lower_adjusted = type_min
        if not quiet:
            ogo.message(f'[WARNING] Lower threshold {lower} < type minimum {type_min}. Adjusted to {type_min}')
    
    if upper > type_max:
        upper_adjusted = type_max
        if not quiet:
            ogo.message(f'[WARNING] Upper threshold {upper} > type maximum {type_max}. Adjusted to {type_max}')
    
    # Get image statistics
    stats = sitk.StatisticsImageFilter()
    stats.Execute(input_image)
    
    # Display threshold information
    if not quiet:
        guard = '!-------------------------------------------------------------------------------'
        print(guard)
        print(f'!> Threshold Parameters:')
        print(f'!>   lower                {lower_adjusted:>12}')
        print(f'!>   upper                {upper_adjusted:>12}')
        print(f'!>   inside               {inside:>12}')
        print(f'!>   outside              {outside:>12}')
        print(guard)
    
    # Apply threshold
    threshold_filter = sitk.BinaryThresholdImageFilter()
    threshold_filter.SetInsideValue(inside)
    threshold_filter.SetOutsideValue(outside)
    threshold_filter.SetLowerThreshold(lower_adjusted)
    threshold_filter.SetUpperThreshold(upper_adjusted)
    output_image = threshold_filter.Execute(input_image)
    
    # Calculate label statistics on the result
    label_stats = sitk.LabelIntensityStatisticsImageFilter()
    label_stats.Execute(output_image, input_image)
    
    n_labels = label_stats.GetNumberOfLabels()
    
    # Display label statistics
    if not quiet:
        if n_labels > 0:
            ogo.message('')
            ogo.message('Label Statistics:')
            ogo.message(f'  {"Label":<20s} {"Volume":<12s} {"N Voxels":<12s} {"Centroid":>19s}')
            
            for label in label_stats.GetLabels():
                # Get label description
                desc = lb.labels_dict.get(label, {}).get('LABEL', 'unknown label')
                
                # Calculate bounding box center
                bbox = label_stats.GetBoundingBox(label)
                center = [
                    bbox[0] + bbox[3] // 2,
                    bbox[1] + bbox[4] // 2,
                    bbox[2] + bbox[5] // 2
                ]
                
                phys_size = label_stats.GetPhysicalSize(label)
                n_pixels = label_stats.GetNumberOfPixels(label)
                
                ogo.message(f'  {desc:<15s}({label:<3d}) '
                           f'{phys_size:>10.1f} mm³ '
                           f'{n_pixels:>10d}   '
                           f'({center[0]:>5d},{center[1]:>5d},{center[2]:>5d})')
        else:
            ogo.message('')
            ogo.message('[INFO] No labels found in output image.')
    
    # Write output
    if not quiet:
        ogo.message(f'Writing output to: {output_file}')
    
    sitk.WriteImage(output_image, output_file)
    
    if not quiet:
        print_image_info(output_file, output_image)
        ogo.message('[DONE]')
    
    return 0


def main():
    # Setup description
    description = '''
Utility to threshold a NIFTI file.

This tool applies binary thresholding to an image, setting pixels within the 
threshold range [lower, upper] to the 'inside' value, and all other pixels to 
the 'outside' value.

The threshold operation uses SimpleITK's BinaryThresholdImageFilter, which 
includes both lower and upper bounds (i.e., lower ≤ pixel ≤ upper).

THRESHOLD BEHAVIOR:
  - Pixels where (lower ≤ value ≤ upper) are set to 'inside' value
  - All other pixels are set to 'outside' value
  
TYPICAL USAGE:
  - Create binary mask: --inside 1 --outside 0
  - Extract specific range: --lower 250 --upper 3000
  - Bone segmentation: --lower 250 --inside 127 --outside 0

JSON CONFIGURATION:
  - To output current arguments as JSON: --json
  - To load arguments from JSON file: --json my_config.json
'''
    epilog = '''
Example calls: 
ogoImageThreshold input.nii.gz output.nii.gz --upper 1500
ogoImageThreshold input.nii.gz output.nii.gz --lower 250 --upper 3000 --inside 1 --outside 0
ogoImageThreshold input.nii.gz output.nii.gz --lower 2500 --upper 5000 --inside 81 --quiet
ogoImageThreshold input.nii.gz output.nii.gz --lower 500 --json > my_config.json
ogoImageThreshold --json my_config.json
'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoImageThreshold",
        description=description,
        epilog=epilog
    )
    
    # JSON configuration argument
    parser.add_argument('--json', nargs='?', const='', metavar='FILE',
                        help='Output arguments as JSON to stdout (no file), or load arguments from JSON file')
    
    parser.add_argument('input_file', nargs='?', type=validate_input_file(['.nii', '.nii.gz']), 
                        help='Input image file (*.nii, *.nii.gz)')
    parser.add_argument('output_file', nargs='?', type=validate_output_file(['.nii', '.nii.gz']), 
                        help='Output image file (*.nii, *.nii.gz)')
    parser.add_argument('--lower','-l', type=int, default=250, metavar='VAL', 
                        help='Lower threshold (default: %(default)s)')
    parser.add_argument('--upper', '-u', type=int, default=3000, metavar='VAL', 
                        help='Upper threshold (default: %(default)s)')
    parser.add_argument('--inside', '-i', type=validate_integer_range(0, 255), default=127, metavar='VAL', 
                        help='Inside value (0-255, default: %(default)s)')
    parser.add_argument('--outside', '-o', type=validate_integer_range(0, 255), default=0, metavar='VAL', 
                        help='Outside value (0-255, default: %(default)s)')
    parser.add_argument('--overwrite', '-ow', action='store_true', 
                        help='Overwrite output without asking')
    parser.add_argument('--quiet', '-q', action='store_true', 
                        help='Suppress informational output to stdout')

    # Parse and display
    args = parser.parse_args()
    
    # Handle --json flag
    if args.json is not None:
        if args.json == '':
            # --json with no argument: output to stdout
            save_config_to_json(args, None)
            return 0
        else:
            # --json file.json: load from file
            if not args.json.endswith('.json'):
                parser.error('JSON file must have .json extension')
            json_config = load_config_from_json(args.json)
            # Update args with JSON values (command line args override JSON)
            for key, value in json_config.items():
                if hasattr(args, key):
                    # Only update if argument was not explicitly provided on command line
                    if getattr(args, key) == parser.get_default(key):
                        setattr(args, key, value)
    
    # Check that required positional arguments are provided
    if args.input_file is None or args.output_file is None:
        parser.error('input_file and output_file are required')
    
    if not args.quiet:
        ogo.message(echo_arguments('ImageThreshold', vars(args)))

    # Run program (exclude json argument)
    run_args = vars(args).copy()
    run_args.pop('json', None)
    result = ImageThreshold(**run_args)
    return result


if __name__ == '__main__':
    main()
