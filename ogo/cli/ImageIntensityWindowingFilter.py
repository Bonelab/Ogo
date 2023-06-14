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

def ImageIntensityWindowingFilter(input_filename, output_filename, output_maximum, output_minimum, window_maximum, window_minimum, overwrite):
    
    # Check if output exists and should overwrite
    if os.path.isfile(output_filename) and not overwrite:
        result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(output_filename))
        if result.lower() not in ['y', 'yes']:
            print('Not overwriting. Exiting...')
            os.sys.exit()
            
    # Read input
    if not os.path.isfile(input_filename):
        os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(input_filename))

    ogo.message('Reading input CT image:')
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
    if output_minimum < type_min:
        output_minimum = type_min
        ogo.message('[WARNING] Invalid output minimum. Value changed to {}'.format(output_minimum))
    if output_maximum > type_max:
        output_maximum = type_max
        ogo.message('[WARNING] Invalid output maximum. Value changed to {}'.format(output_maximum))
    if window_minimum < type_min:
        window_minimum = type_min
        ogo.message('[WARNING] Invalid window minimum. Value changed to {}'.format(window_minimum))
    if window_maximum > type_max:
        window_maximum = type_max
        ogo.message('[WARNING] Invalid window maximum. Value changed to {}'.format(window_maximum))
    
    # Image statistics
    stats = sitk.StatisticsImageFilter()
    stats.Execute(ct)
    
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
    print('!> window_minimum                {:>8}'.format(window_minimum))
    print('!> window_maximum                {:>8}'.format(window_maximum))
    print('!> output_minimum                {:>8}'.format(output_minimum))
    print('!> output_maximum                {:>8}'.format(output_maximum))
    print(guard)
    
    if window_minimum > window_maximum:
        ogo.message('[ERROR] Window minimum cannot be greater than window maximum.')
        os.sys.exit()
    if output_minimum > output_maximum:
        ogo.message('[ERROR] Output minimum cannot be greater than output maximum.')
        os.sys.exit()
        
    ct_scale = sitk.IntensityWindowingImageFilter()
    ct_scale.SetWindowMaximum(window_maximum)
    ct_scale.SetWindowMinimum(window_minimum)
    ct_scale.SetOutputMaximum(output_maximum)
    ct_scale.SetOutputMinimum(output_minimum)
    ct_out = ct_scale.Execute(ct)
    
    # Report information about output image
    stats.Execute(ct_out)
    dim = ct_out.GetSize()
    spacing = ct_out.GetSpacing()
    origin = ct_out.GetOrigin()
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
        
    ogo.message('Writing output to file {}'.format(output_filename))
    sitk.WriteImage(ct_out, output_filename)

def main():
    # Setup description
    description = '''
Applies a linear transformation to the intensity levels of the input Image 
that are inside a user-defined interval. Values below this interval are 
mapped to a constant. Values over the interval are mapped to another constant.
'''
    epilog = '''
Example calls: 
ogoImageIntensityWindowingFilter input.nii.gz output.nii.gz
'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoImageIntensityWindowingFilter",
        description=description,
        epilog=epilog
    )
    parser.add_argument('input_filename', help='Input image file (*.nii, *.nii.gz)')
    parser.add_argument('output_filename', help='Output image file (*.nii, *.nii.gz)')
    parser.add_argument('--output_maximum', type=int, default=1200, metavar='VAL', help='Maximum output image intensity (default: %(default)s)')
    parser.add_argument('--output_minimum', type=int, default=-400, metavar='VAL', help='Minimum output image intensity (default: %(default)s)')
    parser.add_argument('--window_maximum', type=int, default=150, metavar='VAL', help='Maximum input window image intensity (default: %(default)s)')
    parser.add_argument('--window_minimum', type=int, default=-200, metavar='VAL', help='Minimum input window image intensity (default: %(default)s)')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite output without asking')

    # Parse and display
    args = parser.parse_args()
    print(echo_arguments('ImageIntensityWindowingFilter', vars(args)))

    # Run program
    ImageIntensityWindowingFilter(**vars(args))


if __name__ == '__main__':
    main()
