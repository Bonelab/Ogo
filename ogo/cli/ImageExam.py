# /------------------------------------------------------------------------------+
# | 19-JUL-2022                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+

# Imports
import argparse
import os
import vtk
import numpy as np
import math
from vtk.util.numpy_support import vtk_to_numpy
from ogo.util.echo_arguments import echo_arguments
import ogo.util.Helper as ogo
import SimpleITK as sitk


def ImageExam(input_filename, header, histogram, bins):
    # Check valid number of histogram bins
    if bins < 1:
        os.sys.exit('[ERROR] Number of histogram bins must be greater than zero: bins = {}'.format(bins))

    # ogo.message('Start ogoImageExam!')

    # Read input
    if not os.path.isfile(input_filename):
        os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(input_filename))

    if input_filename.lower().endswith('.nii'):
        reader = vtk.vtkNIFTIImageReader()
    elif input_filename.lower().endswith('.nii.gz'):
        reader = vtk.vtkNIFTIImageReader()
    else:
        os.sys.exit('[ERROR] Cannot find reader for file \"{}\"'.format(input_filename))

    print()
    ogo.message('Reading input image ' + input_filename)
    reader.SetFileName(input_filename)
    reader.Update()

    # ogo.message('Image dimensions and type')
    ogo.aix(input_filename, reader.GetOutput())

    if header:
        ogo.message('Header from NIFTI file (VTK reader)')
        ogo.infoNIFTI(reader)
        
        ogo.message('Header from NIFTI file –-> Using SimpleITK reader')
        ct = sitk.ReadImage(input_filename)
        for k in ct.GetMetaDataKeys(): 
            v = ct.GetMetaData(k) 
            print('  {:20s} = {}'.format(k,v))
        
    if histogram:
        ogo.message('Histogram of NIFTI image data')
        ogo.histogram(reader.GetOutput(), bins)
    
    
    
    # ogo.message('Done ogoImageExam!')


def main():
    # Setup description
    description = '''
Utility to read a NIFTI file and report characteristics such as image type, 
dimensions, histograms, and header.
'''
    epilog = '''
Example calls: 
ogoImageExam input.nii.gz
ogoImageExam --histogram input.nii.gz
'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoImageExam",
        description=description,
        epilog=epilog
    )
    parser.add_argument('input_filename', help='Input image file (*.nii, *.nii.gz)')
    parser.add_argument('--header', action='store_true', help='Show NIFTI header (default: %(default)s)')
    parser.add_argument('--histogram', action='store_true', help='Show histogram (default: %(default)s)')
    parser.add_argument('--bins', type=int, default=128, metavar='#',
                        help='Number of bins in histogram (default: %(default)s)')

    # Parse and display
    args = parser.parse_args()
    # print(echo_arguments('ImageExam', vars(args)))

    # Run program
    ImageExam(**vars(args))


if __name__ == '__main__':
    main()
