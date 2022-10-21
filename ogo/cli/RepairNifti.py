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

from bonelab.util.echo_arguments import echo_arguments
from bonelab.util.time_stamp import message
from vtk.util.numpy_support import vtk_to_numpy
from bonelab.io.vtk_helpers import get_vtk_reader, get_vtk_writer, handle_filetype_writing_special_cases

def CheckExt(choices):
    class Act(argparse.Action):
        def __call__(self,parser,namespace,fname,option_string=None):
            ext = os.path.splitext(fname)[1][1:]
            if ext not in choices:
                option_string = '({})'.format(option_string) if option_string else ''
                parser.error("File extension must be {}{}".format(choices,option_string))
            else:
                setattr(namespace,self.dest,fname)

    return Act

def repair_nifti(input_filename, output_filename, overwrite=False):

    debug = False
    
    # Read input
    if not os.path.isfile(input_filename):
        os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(input_filename))

    message('Inputs:',
            '{:16s} = {:16s}'.format('input_filename',input_filename),
            '{:16s} = {:16s}'.format('output_filename',output_filename),
            '{:16s} = {:16s}'.format('overwrite',('True' if overwrite else 'False'))
            )

    # Check if output exists
    if os.path.isfile(output_filename) and not overwrite:
      result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(output_filename))
      if result.lower() not in ['y', 'yes']:
        print('Not overwriting. Exiting...')
        os.sys.exit()

    im = sitk.ReadImage(input_filename)
    if (debug): print(im)
    print('-------------------------------------------------------------------------------')
    print('{:16s} = {:16s}'.format('input_filename',input_filename))
    print('{:16s} = {:16.10f} {:16.10f} {:16.10f}'.format('spacing',im.GetSpacing()[0],im.GetSpacing()[1],im.GetSpacing()[2]))
    print('{:16s} = {:16.10f} {:16.10f} {:16.10f}'.format('origin',im.GetOrigin()[0],im.GetOrigin()[1],im.GetOrigin()[2]))
    print('{:16s} = {:16.10f} {:16.10f} {:16.10f}'.format('direction',im.GetDirection()[0],im.GetDirection()[1],im.GetDirection()[2]))
    print('{:16s}   {:16.10f} {:16.10f} {:16.10f}'.format('',im.GetDirection()[3],im.GetDirection()[4],im.GetDirection()[5]))
    print('{:16s}   {:16.10f} {:16.10f} {:16.10f}'.format('',im.GetDirection()[6],im.GetDirection()[7],im.GetDirection()[8]))

    # If user wants to output an image
    if (output_filename != 'None'):
      # Create a new image
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
A utility to repair a known issues with some nifti headers.
The classic error seems to be:

  ITK ERROR: ITK only supports orthonormal direction cosines.  
  No orthonormal definition found!

This utility will correct the issue.

Example call:
  python ogoRepairNifty.py /Users/skboyd/Documents/projects/OppScreening/work/MODELS/RETRO_00037.nii
  python ogoRepairNifty.py ./MODELS/RETRO_00037.nii
'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="RepairNifty",
        description=description,
    )
    parser.add_argument('input_filename', help='Input NIFTI file')
    parser.add_argument('--output_filename',default="None",action=CheckExt({'nii','NII'}), metavar='FILE', help='Output NIFTI file (default: %(default)s)')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite output without asking')

    # Parse and display
    args = parser.parse_args()
    print(echo_arguments('DebugNifty', vars(args)))
        
    # Run program
    repair_nifti(**vars(args))

if __name__ == '__main__':
    main()
