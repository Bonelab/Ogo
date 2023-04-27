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

def ImageThreshold(input_filename, output_filename, lower_threshold, upper_threshold, inside_value, outside_value, overwrite):
    
    # Check if output exists and should overwrite
    if os.path.isfile(output_filename) and not overwrite:
        result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(output_filename))
        if result.lower() not in ['y', 'yes']:
            print('Not overwriting. Exiting...')
            os.sys.exit()
            
    # Read input
    if not os.path.isfile(input_filename):
        os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(input_filename))

    ogo.message('Reading input CT image to be cropped:')
    ogo.message('      \"{}\"'.format(input_filename))
    ct = sitk.ReadImage(input_filename)

    # Checking image type
    ogo.message('Input image type is {}'.format(ct.GetPixelIDTypeAsString()))
    image_type = ct.GetPixelIDValue()
    if image_type is sitk.sitkInt8:
        type_min = -128
        type_max = 127
    elif image_type is sitk.sitkUInt8:
        type_min = 0
        type_max = 255
    elif image_type is sitk.sitkInt16:
        type_min = -32768
        type_max = 32767
    elif image_type is sitk.sitkUInt16:
        type_min = 0
        type_max = 65535
    elif image_type is sitk.sitkInt32:
        type_min = -2147483648
        type_max = 2147483647
    elif image_type is sitk.sitkUInt32:
        type_min = 0
        type_max = 4294967295
    else:
        ogo.message('[ERROR] Unexpected image input type.')
    if lower_threshold < type_min:
        lower_threshold = type_min
        ogo.message('[WARNING] Invalid lower threshold. Changed to {}'.format(lower_threshold))
    if upper_threshold > type_max:
        upper_threshold = type_max
        ogo.message('[WARNING] Invalid upper threshold. Changed to {}'.format(upper_threshold))
    
    # Image statistics
    stats = sitk.StatisticsImageFilter()
    stats.Execute(ct)
    
    if lower_threshold < stats.GetMinimum():
        ogo.message('Lower ')
    # Report information about input image
    dim = ct.GetSize()
    spacing = ct.GetSpacing()
    origin = ct.GetOrigin()
    phys_dim = [x * y for x, y in zip(dim, spacing)]
    position = [math.floor(x / y) for x, y in zip(origin, spacing)]
    
    guard = '!-------------------------------------------------------------------------------'
    print(guard)
    print('!> dim                            {:>8}  {:>8}  {:>8}'.format(*dim))
    print('!> off                            {:>8}  {:>8}  {:>8}'.format('-', '-', '-'))
    print('!> pos                            {:>8}  {:>8}  {:>8}'.format(*position))
    print('!> element size in mm             {:>8.4f}  {:>8.4f}  {:>8.4f}'.format(*spacing))
    print('!> phys dim in mm                 {:>8.4f}  {:>8.4f}  {:>8.4f}'.format(*phys_dim))
    print('!> ')
    print('!> min                            {:>8.0f}'.format(stats.GetMinimum()))
    print('!> max                            {:>8.0f}'.format(stats.GetMaximum()))
    print(guard)
    
    # Threshold
    print(guard)
    print('!> lower_threshold                {:>8}'.format(lower_threshold))
    print('!> upper_threshold                {:>8}'.format(upper_threshold))
    print('!> inside_value                   {:>8}'.format(inside_value))
    print('!> outside_value                  {:>8}'.format(outside_value))
    print(guard)
    
    if lower_threshold > upper_threshold:
        os.sys.exit('[ERROR] Lower threshold cannot be greater than upper threshold.')
    if inside_value == outside_value:
        os.message('[WARNING] Setting inside_value equal to outside_value makes no sense!')
        
    ct_thres = sitk.BinaryThresholdImageFilter()
    ct_thres.SetInsideValue(inside_value)
    ct_thres.SetOutsideValue(outside_value)
    ct_thres.SetLowerThreshold(lower_threshold)
    ct_thres.SetUpperThreshold(upper_threshold)
    ct_out = ct_thres.Execute(ct)
        
    ogo.message('Writing output to file {}'.format(output_filename))
    sitk.WriteImage(ct_out, output_filename)
    
    # Report information about input image
    dim = ct_out.GetSize()
    spacing = ct_out.GetSpacing()
    origin = ct_out.GetOrigin()
    phys_dim = [x * y for x, y in zip(dim, spacing)]
    position = [math.floor(x / y) for x, y in zip(origin, spacing)]
    
    ogo.message('Done ogoImageThreshold!')


def main():
    # Setup description
    description = '''
Utility to threshold a NIFTI file.
'''
    epilog = '''
Example calls: 
ogoImageThreshold --upper_threshold 1500 input.nii.gz output.nii.gz
'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoImageThreshold",
        description=description,
        epilog=epilog
    )
    parser.add_argument('input_filename', help='Input image file (*.nii, *.nii.gz)')
    parser.add_argument('output_filename', help='Output image file (*.nii, *.nii.gz)')
    parser.add_argument('--lower_threshold', type=int, default=-32768, metavar='VAL', help='Lower threshold (default: %(default)s)')
    parser.add_argument('--upper_threshold', type=int, default=32767, metavar='VAL', help='Upper threshold (default: %(default)s)')
    parser.add_argument('--inside_value', type=int, default=0, metavar='VAL', help='Inside value (default: %(default)s)')
    parser.add_argument('--outside_value', type=int, default=127, metavar='VAL', help='Outside value (default: %(default)s)')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite output without asking')

    # Parse and display
    args = parser.parse_args()
    print(echo_arguments('ImageThreshold', vars(args)))

    # Run program
    ImageThreshold(**vars(args))


if __name__ == '__main__':
    main()
