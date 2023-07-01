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
import nibabel as nib

from ogo.util.echo_arguments import echo_arguments
import ogo.util.Helper as ogo
from vtk.util.numpy_support import vtk_to_numpy

def print_info(ct):
    ogo.message('  {:>12s}'.format('spacing:')+', '.join('{:12.6f}'.format(i) for i in ct.GetSpacing()))
    ogo.message('  {:>12s}'.format('origin:')+', '.join('{:12.6f}'.format(i) for i in ct.GetOrigin()))
    ogo.message('  {:>12s}'.format('direction:')+'{:12.6f}, {:12.6f}, {:12.6f},'.format(ct.GetDirection()[0],ct.GetDirection()[1],ct.GetDirection()[2]))
    ogo.message('  {:>12s}'.format(' ')         +'{:12.6f}, {:12.6f}, {:12.6f},'.format(ct.GetDirection()[3],ct.GetDirection()[4],ct.GetDirection()[5]))
    ogo.message('  {:>12s}'.format(' ')         +'{:12.6f}, {:12.6f}, {:12.6f} '.format(ct.GetDirection()[6],ct.GetDirection()[7],ct.GetDirection()[8]))

def print_info_nib(ct):
    hdr = ct.header
    ogo.message('  {:>10s}'.format('spacing:')+'{:12.6f}, {:12.6f}, {:12.6f}'.format(hdr['pixdim'][1],hdr['pixdim'][2],hdr['pixdim'][3]))
    ogo.message('  {:>10s}'.format('origin:')+'{:12.6f}, {:12.6f}, {:12.6f}'.format(hdr['qoffset_x'],hdr['qoffset_y'],hdr['qoffset_z']))
    ogo.message('  {:>10s}'.format('srow_x:')+'{:12.6f}, {:12.6f}, {:12.6f}'.format(hdr['srow_x'][0],hdr['srow_x'][1],hdr['srow_x'][2]))
    ogo.message('  {:>10s}'.format('srow_y:')+'{:12.6f}, {:12.6f}, {:12.6f}'.format(hdr['srow_y'][0],hdr['srow_y'][1],hdr['srow_y'][2]))
    ogo.message('  {:>10s}'.format('srow_z:')+'{:12.6f}, {:12.6f}, {:12.6f}'.format(hdr['srow_z'][0],hdr['srow_z'][1],hdr['srow_z'][2]))
    #print(hdr)
    
def repair_NIfTI(infile, match, outfile, mode, overwrite=False):

    debug = False
    
    # Set output filename to input filename if not already set
    if outfile is None:
        outfile = infile
        ogo.message('Output is set to input filename:\n         \"{}\"'.format(outfile))
        
    # Check if output exists and should overwrite
    if os.path.isfile(outfile) and not overwrite:
        result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(outfile))
        if result.lower() not in ['y', 'yes']:
            ogo.message('Not overwriting. Exiting...')
            os.sys.exit()
    
    if not (outfile.lower().endswith('.nii') or outfile.lower().endswith('.nii.gz')):
        os.sys.exit('[ERROR] Output must be type NIfTI file: \"{}\"'.format(infile))
    
    ogo.message('Mode is: {}'.format(mode))

    # Read input image
    if not os.path.isfile(infile):
        os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(infile))
    
    if not (infile.lower().endswith('.nii') or infile.lower().endswith('.nii.gz')):
        os.sys.exit('[ERROR] Input must be type NIfTI file: \"{}\"'.format(infile))
    
    
    if mode == 'orthonormal':
            
        ct = sitk.ReadImage(infile)
            
        if (debug): print(ct)
        ogo.message('-----------------------------------------------------------------------')
        ogo.message('{:16}  {}'.format('Input:',infile))        
        print_info(ct)
        ogo.message('-----------------------------------------------------------------------')
        
        # Create output image
        nda = sitk.GetArrayFromImage(ct)
        
        ct_new = sitk.GetImageFromArray(nda)
        ct_new.SetSpacing(ct.GetSpacing())
        ct_new.SetOrigin(ct.GetOrigin())
        ct_new.SetDirection(ct.GetDirection())
        
        if (debug): print(ct_new)
        ogo.message('-----------------------------------------------------------------------')
        ogo.message('{:16}  {}'.format('Output:',outfile))        
        print_info(ct_new)
        ogo.message('-----------------------------------------------------------------------')
          
        sitk.WriteImage(ct_new, outfile)
        
    elif mode == 'affine':
        
        #ct = sitk.ReadImage(infile)
        ct = nib.load(infile)
        
        if match is None:
            ogo.message('The matching NIfTI file must be set with --match. Exiting...')
            os.sys.exit()
        
        if not (match.lower().endswith('.nii') or match.lower().endswith('.nii.gz')):
            os.sys.exit('[ERROR] Match must be type NIfTI file: \"{}\"'.format(match))
        
        # Uses SITK
        #match_ct = sitk.ReadImage(match)
        match_ct = nib.load(match)
        
        if (debug): print(ct)
        ogo.message('-----------------------------------------------------------------------')
        ogo.message('{:16}  {}'.format('Input:',infile))        
        print_info_nib(ct)
        ogo.message('-----------------------------------------------------------------------')
        
        if (debug): print(match_ct)
        ogo.message('-----------------------------------------------------------------------')
        ogo.message('{:16}  {}'.format('Match:',match))        
        print_info_nib(match_ct)
        ogo.message('-----------------------------------------------------------------------')
        
        # Create output image
        data = ct.get_data()
                
        ct_new = nib.Nifti1Image(data, match_ct.affine, match_ct.header.copy())
        #ct_new = nib.NIfTI1Image(data, None, header=new_header)
        

        if (debug): print(ct_new)
        ogo.message('-----------------------------------------------------------------------')
        ogo.message('{:16}  {}'.format('Output:',outfile))        
        print_info_nib(ct_new)
        ogo.message('-----------------------------------------------------------------------')
          
        nib.save(ct_new, outfile)
        
    else:
        ogo.message('[ERROR] User has not defined a valid mode.')
        os.sys.exit()
        
def main():
    # Setup description
    description='''
A utility to repair NIfTI files. There are various modes of
operation:

orthonormal

A known issue with some NIfTI headers that make it incompatible 
with SimpleITK. The typical error message is:

  ITK ERROR: ITK only supports orthonormal direction cosines.  
  No orthonormal definition found!

affine

This will copy the affine transform from one image to another.
It is useful in case a raw CT and matching labels have different
affine transforms.

'''
    epilog = '''
Example calls: 
  ogoRepairNIfTI ./MODELS/RETRO_00037.nii --outfile \
                 ./MODELS/RETRO_00037_FIXED.nii
  ogoRepairNIfTI ./MODELS/RETRO_00037.nii.gz
'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="RepairNIfTI",
        description=description,
        epilog=epilog
    )
    parser.add_argument('infile', metavar='FILE', help='Input NIfTI file')
    parser.add_argument('--match', metavar='FILE', default=None, help='Use with affine mode (default: %(default)s)')
    parser.add_argument('--outfile', metavar='FILE', default=None, help='Output NIfTI file (default: %(default)s)')
    parser.add_argument('--mode', default='orthonormal', choices=['affine', 'orthonormal'],help='Correction mode (default: %(default)s)')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite output without asking')

    # Parse and display
    args = parser.parse_args()
    print(echo_arguments('RepairNIfTI', vars(args)))
        
    # Run program
    repair_NIfTI(**vars(args))

if __name__ == '__main__':
    main()
