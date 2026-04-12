# /------------------------------------------------------------------------------+
# | 19-JUL-2022                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+

# Imports
import argparse
import os
import numpy as np
from ogo.util.echo_arguments import echo_arguments
import ogo.util.Helper as ogo
import SimpleITK as sitk
from ogo.util.cli_tools import (
    validate_integer_min,
    validate_input_file,
    print_image_info
)


def ImageExam(input_file, header, histogram, bins, quiet):
    # Check valid number of histogram bins
    if bins < 1:
        ogo.message('[ERROR] Number of histogram bins must be greater than zero: bins = {}'.format(bins))
        return 1

    # Read input using SimpleITK
    if not quiet:
        ogo.message('Reading input image ' + input_file)
    
    image = sitk.ReadImage(input_file)
    
    if not quiet:
        print_image_info(input_file, image)

    if header:
        if not quiet:
            ogo.message('Header from NIFTI file (SimpleITK reader)')
        
        for k in image.GetMetaDataKeys(): 
            v = image.GetMetaData(k) 
            if not quiet:
                print('  {:20s} = {}'.format(k, v))
        
    if histogram:
        if not quiet:
            ogo.message('Histogram of NIFTI image data')
            
            # Get numpy array from SimpleITK image
            array = sitk.GetArrayFromImage(image).ravel()
            
            # Compute histogram
            counts, bin_edges = np.histogram(array, bins=bins)
            
            # Display histogram
            guard = '!-------------------------------------------------------------------------------'
            line_format_header = '!>  {:>8s} {:>8s} {:>12s}'
            line_format_data = '!>  {:>8d} {:>8.3f} {:>12d}'
            line_format_summary = '!>  {:>8s} {:>8.3f} {:>12s}'
            
            print(guard)
            print(line_format_header.format('Bin', 'Value', 'Count'))
            print(guard)
            
            # Print histogram bins
            for i in range(bins):
                bin_center = (bin_edges[i] + bin_edges[i + 1]) / 2.0
                print(line_format_data.format(i, bin_center, counts[i]))
            
            print(guard)
            print(line_format_summary.format('Min', np.min(array), ''))
            print(line_format_summary.format('Max', np.max(array), ''))
            print(line_format_summary.format('Mean', np.mean(array), ''))
            print(line_format_summary.format('Std', np.std(array), ''))
            print(guard)
    
    return 0


def main():
    # Setup description
    description = '''
Utility to read a NIFTI file and report characteristics such as image type, 
dimensions, histograms, and header.
'''
    epilog = '''
Example calls: 
ogoImageExam input.nii.gz
ogoImageExam input.nii.gz --histogram
ogoImageExam input.nii.gz -H --histogram -b 256 -q
'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoImageExam",
        description=description,
        epilog=epilog
    )
    parser.add_argument('input_file', type=validate_input_file(['.nii', '.nii.gz']), help='Input image file (*.nii, *.nii.gz)')
    parser.add_argument('--header', '-H', action='store_true', help='Show NIFTI header (default: %(default)s)')
    parser.add_argument('--histogram', action='store_true', help='Show histogram (default: %(default)s)')
    parser.add_argument('--bins', '-b', type=validate_integer_min(1), default=128, metavar='#',
                        help='Number of bins in histogram (default: %(default)s)')
    parser.add_argument('--quiet', '-q', action='store_true', help='Suppress informational output to stdout')

    # Parse and display
    args = parser.parse_args()
    if not args.quiet:
        print(echo_arguments('ImageExam', vars(args)))

    # Run program
    result = ImageExam(**vars(args))
    return result


if __name__ == '__main__':
    main()
