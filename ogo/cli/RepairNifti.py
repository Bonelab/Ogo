# /------------------------------------------------------------------------------+
# | 21-OCT-2022                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+

#Imports
import argparse
import os
import sys
import time
import math
import vtk
import vtkbone
import numpy as np
import SimpleITK as sitk

from ogo.util.echo_arguments import echo_arguments
import ogo.util.Helper as ogo
from vtk.util.numpy_support import vtk_to_numpy

def repair_nifti(input_filename, output_filename, overwrite=False):

    debug = False
    
    # Set output filename to input filename if not already set
    if output_filename is None:
        output_filename = input_filename
        ogo.message('Output is set to input filename:\n         \"{}\"'.format(output_filename))
        
    # Check if output exists and should overwrite
    if os.path.isfile(output_filename) and not overwrite:
        result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(output_filename))
        if result.lower() not in ['y', 'yes']:
            ogo.message('Not overwriting. Exiting...')
            os.sys.exit()
    
    if not (output_filename.lower().endswith('.nii') or output_filename.lower().endswith('.nii.gz')):
        os.sys.exit('[ERROR] Output must be type NIFTI file: \"{}\"'.format(input_filename))
    
    # Read input image
    if not os.path.isfile(input_filename):
        os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(input_filename))

    if not (input_filename.lower().endswith('.nii') or input_filename.lower().endswith('.nii.gz')):
        os.sys.exit('[ERROR] Input must be type NIFTI file: \"{}\"'.format(input_filename))
    
    im = sitk.ReadImage(input_filename)

    if (debug): print(im)
    print('-------------------------------------------------------------------------------')
    print('{:16s} = {:16s}'.format('input_filename',input_filename))
    print('{:16s} = {:16.10f} {:16.10f} {:16.10f}'.format('spacing',im.GetSpacing()[0],im.GetSpacing()[1],im.GetSpacing()[2]))
    print('{:16s} = {:16.10f} {:16.10f} {:16.10f}'.format('origin',im.GetOrigin()[0],im.GetOrigin()[1],im.GetOrigin()[2]))
    print('{:16s} = {:16.10f} {:16.10f} {:16.10f}'.format('direction',im.GetDirection()[0],im.GetDirection()[1],im.GetDirection()[2]))
    print('{:16s}   {:16.10f} {:16.10f} {:16.10f}'.format('',im.GetDirection()[3],im.GetDirection()[4],im.GetDirection()[5]))
    print('{:16s}   {:16.10f} {:16.10f} {:16.10f}'.format('',im.GetDirection()[6],im.GetDirection()[7],im.GetDirection()[8]))

    # Create output image
    nda = sitk.GetArrayFromImage(im)
    
    im_new = sitk.GetImageFromArray(nda)
    im_new.SetSpacing(im.GetSpacing())
    im_new.SetOrigin(im.GetOrigin())
    im_new.SetDirection(im.GetDirection())
    
    if (debug): print(im_new)
    
    print('-------------------------------------------------------------------------------')
    print('{:16s} = {:16s}'.format('output_filename',output_filename))
    print('{:16s} = {:16.10f} {:16.10f} {:16.10f}'.format('spacing',im_new.GetSpacing()[0],im_new.GetSpacing()[1],im_new.GetSpacing()[2]))
    print('{:16s} = {:16.10f} {:16.10f} {:16.10f}'.format('origin',im_new.GetOrigin()[0],im_new.GetOrigin()[1],im_new.GetOrigin()[2]))
    print('{:16s} = {:16.10f} {:16.10f} {:16.10f}'.format('direction',im_new.GetDirection()[0],im_new.GetDirection()[1],im_new.GetDirection()[2]))
    print('{:16s}   {:16.10f} {:16.10f} {:16.10f}'.format('',im_new.GetDirection()[3],im_new.GetDirection()[4],im_new.GetDirection()[5]))
    print('{:16s}   {:16.10f} {:16.10f} {:16.10f}'.format('',im_new.GetDirection()[6],im_new.GetDirection()[7],im_new.GetDirection()[8]))
      
    sitk.WriteImage(im_new, output_filename)

def main():
    # Setup description
    description='''
A utility to repair a known issue with some nifti headers that
make it incompatible with SimpleITK. The classic error message is:

  ITK ERROR: ITK only supports orthonormal direction cosines.  
  No orthonormal definition found!

This utility will correct the issue.
'''
    epilog = '''
Example calls: 
  ogoRepairNifti ./MODELS/RETRO_00037.nii --output_filename \
                 ./MODELS/RETRO_00037_FIXED.nii
  ogoRepairNifti ./MODELS/RETRO_00037.nii.gz
'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="RepairNifti",
        description=description,
        epilog=epilog
    )
    parser.add_argument('input_filename', help='Input NIFTI file')
    parser.add_argument('--output_filename',default=None,metavar='FILE', help='Output NIFTI file (default: %(default)s)')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite output without asking')

    # Parse and display
    args = parser.parse_args()
    print(echo_arguments('RepairNifti', vars(args)))
        
    # Run program
    repair_nifti(**vars(args))

if __name__ == '__main__':
    main()
