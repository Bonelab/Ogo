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
            
            # Check if data type is integer
            is_integer_type = np.issubdtype(array.dtype, np.integer)
            
            # Compute histogram with integer-aligned bins for integer data
            if is_integer_type:
                # For integer data, align bin edges to integer values
                # Use strict [a, b) intervals by extending beyond max
                data_min = int(np.min(array))
                data_max = int(np.max(array))
                
                # Create integer bin edges from min to max+1
                # The last edge at max+1 ensures that max value falls in [max, max+1)
                # This enforces strict [a, b) for all bins (no special case for last bin)
                bin_edges = np.linspace(data_min, data_max + 1, bins + 1, dtype=float)
                bin_edges = np.round(bin_edges).astype(int)
                
                # Ensure we don't have duplicate edges
                bin_edges = np.unique(bin_edges)
                if len(bin_edges) < bins + 1:
                    # If unique edges are fewer than requested bins, use what we have
                    if not quiet:
                        ogo.message(f'[INFO] Adjusted to {len(bin_edges)-1} bins for integer alignment')
                
                counts, bin_edges = np.histogram(array, bins=bin_edges)
            else:
                # For floating-point data, use default binning
                counts, bin_edges = np.histogram(array, bins=bins)
            
            # Display histogram
            guard = '!-------------------------------------------------------------------------------'
            
            # Print histogram bins differently for integer vs float data
            actual_bins = len(counts)
            
            if is_integer_type:
                # For integer data, show bin range instead of center
                # All bins use [min, max) notation (right edge excluded)
                line_format_header = '!>  {:>18s} {:>8s} {:>12s}'
                line_format_data = '!>  [{:>7.0f}, {:>7.0f}) {:>8d} {:>12d}'
                line_format_summary = '!>  {:>18.0f} {:>8s} {:>12s}'
                
                print(guard)
                print(line_format_header.format('Range [min, max)', 'Bin', 'Count'))
                print(guard)
                
                for i in range(actual_bins):
                    print(line_format_data.format(bin_edges[i], bin_edges[i + 1], i, counts[i]))
            else:
                # For float data, show bin center
                line_format_header = '!>  {:>8s} {:>8s} {:>12s}'
                line_format_data = '!>  {:>8.3f} {:>8d} {:>12d}'
                line_format_summary = '!>  {:>8.3f} {:>8s} {:>12s}'
                
                print(guard)
                print(line_format_header.format('Value', 'Bin', 'Count'))
                print(guard)
                
                for i in range(actual_bins):
                    bin_center = (bin_edges[i] + bin_edges[i + 1]) / 2.0
                    print(line_format_data.format(bin_center, i, counts[i]))
            
            print(guard)
            print(line_format_summary.format(np.min(array), 'Min', ''))
            print(line_format_summary.format(np.max(array), 'Max', ''))
            print(line_format_summary.format(np.mean(array), 'Mean', ''))
            print(line_format_summary.format(np.std(array), 'Std', ''))
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
