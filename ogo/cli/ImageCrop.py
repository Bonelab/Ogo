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

def ImageCrop(input_filename, output_filename, overwrite, voi):
    
    # Check if output exists and should overwrite
    if os.path.isfile(output_filename) and not overwrite:
        result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(output_filename))
        if result.lower() not in ['y', 'yes']:
            print('Not overwriting. Exiting...')
            os.sys.exit()

    # Read input
    if not os.path.isfile(input_filename):
        os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(input_filename))

    ogo.message('Reading input CT image to be calibrated:')
    ogo.message('      \"{}\"'.format(input_filename))
    ct = sitk.ReadImage(input_filename)
    
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
    print(guard)
    
    # Take the subvolume
    if voi[0]<0 or voi[2]<0 or voi[4]<0 or voi[1]>(dim[0]-1) or voi[3]>(dim[1]-1) or voi[5]>(dim[2]-1):
        os.sys.exit('[ERROR] Defined VOI is out of range: {} {} {} {} {} {}'.format(*voi))
    
    ct_out = ct[voi[0]:voi[1], voi[2]:voi[3], voi[4]:voi[5]]

    ogo.message('Writing output to file {}'.format(output_filename))
    sitk.WriteImage(ct_out, output_filename)
    
    ogo.message('Done ogoImageCrop!')


def main():
    # Setup description
    description = '''
Utility to crop a NIFTI file to a subvolume.
'''
    epilog = '''
Example calls: 
ogoImageCrop --voi 0 511 0 511 137 448 input.nii.gz output.nii.gz

python ImageCrop.py \
/Users/skboyd/Documents/projects/CTDXAICI/models/CTDXAICI_0053_V01.nii.gz \
/Users/skboyd/Documents/projects/CTDXAICI/models/CTDXAICI_0053_V01_sub.nii.gz \
--voi 0 511 0 511 137 448 
'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoImageCrop",
        description=description,
        epilog=epilog
    )
    parser.add_argument('input_filename', help='Input image file (*.nii, *.nii.gz)')
    parser.add_argument('output_filename', help='Output image file (*.nii, *.nii.gz)')
    parser.add_argument('--voi', type=int, nargs=6, default=[0,1,0,1,0,1], metavar='0', help='VOI bounds in units pixels (default: %(default)s)')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite output without asking')

    # Parse and display
    args = parser.parse_args()
    print(echo_arguments('ImageCrop', vars(args)))

    # Run program
    ImageCrop(**vars(args))


if __name__ == '__main__':
    main()
