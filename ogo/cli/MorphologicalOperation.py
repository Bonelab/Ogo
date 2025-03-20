# /------------------------------------------------------------------------------+
# | 09-MAR-2023                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+

script_version = 1.0

# Imports
import argparse
import os
import sys
import vtk
import math
import numpy as np
import SimpleITK as sitk
from scipy.spatial import procrustes
from datetime import date
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from ogo.util.echo_arguments import echo_arguments
import ogo.util.Helper as ogo
import ogo.dat.OgoMasterLabels as lb

def get_labels(ct):
    filt = sitk.LabelShapeStatisticsImageFilter()
    filt.Execute(ct)
    labels = filt.GetLabels()
    return labels

# +------------------------------------------------------------------------------+
def MorphologicalOperation(input_image, output_image, operation_labels, operation, kernel, overwrite, kernels=None):
    # Check if output exists and should overwrite
    if os.path.isfile(output_image) and not overwrite:
        result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(output_image))
        if result.lower() not in ['y', 'yes']:
            ogo.message('Not overwriting. Exiting...')
            os.sys.exit()

    # Read input
    if not os.path.isfile(input_image):
        os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(input_image))

    if not (input_image.lower().endswith('.nii') or input_image.lower().endswith('.nii.gz')):
        os.sys.exit('[ERROR] Input must be type NIFTI file: \"{}\"'.format(input_image))
    
    ogo.message('Reading image: ')
    ogo.message('\"{}\"'.format(input_image))
    ct = sitk.ReadImage(input_image, sitk.sitkUInt8)
    
    # Get all the available labels in the input image
    labels = get_labels(ct)
    n_labels = len(labels)
    ogo.message('Input image contains the following labels:')
    ogo.message('  [' + ', '.join('{:d}'.format(i) for i in labels) + ']')
    
    # Set the labels on which we will perform morphological operations
    if not operation_labels:
        operation_labels = labels
    else:
        for label in operation_labels:
            if label not in labels:
                os.sys.exit('[ERROR] Label {} is not in input file.'.format(label))
    
    # Validate kernels if mulitple provided
    if kernels:
        if len(kernels) != len(operation_labels):
            os.sys.exit('[ERROR] Number of kernels must match number of operation labels.')
    else:
        kernels = [kernel] * len(operation_labels)
    
    # Set the morphological operation
    ogo.message('Setting morphological operation to: [{}]'.format(operation))
    if operation == 'dilate':
        morphFilt = sitk.BinaryDilateImageFilter()
    elif operation == 'erode':
        morphFilt = sitk.BinaryErodeImageFilter()
    elif operation == 'opening':
        morphFilt = sitk.BinaryMorphologicalOpeningImageFilter()
    elif operation == 'closing':
        morphFilt = sitk.BinaryMorphologicalClosingImageFilter()
    else:
        os.sys.exit('[ERROR] Unknown morphological operation: {}'.format(operation))
    
    # Cycle through the labels
    seg = ct<0 # Zero all labels
        
    for idx, label in enumerate(labels):
        ct_thres = ct==label
        if label in operation_labels:
            kernel_idx = operation_labels.index(label)
            current_kernel = kernels[kernel_idx]
            
            if current_kernel <= 0:
                os.sys.exit('[ERROR] Kernel must be positive: {}'.format(current_kernel))
            
            morphFilt.SetKernelRadius(current_kernel)
            ogo.message('  label {:3d}: {:s}, kernel={:d}'.format(label, operation, current_kernel))
            ct_part = morphFilt.Execute(ct_thres)
        else:
            ogo.message('  label {:3d}: {:s}'.format(label, 'skip'))
            ct_part = ct_thres
            
        bin_seg = seg>0
        bin_part = ct_part>0
        overlap = (bin_seg + bin_part)==2
        mask = 1 - overlap
        this_label = sitk.Mask(ct_part, mask)
        seg = seg + label*(this_label>0)
        
    ct_final = seg
    
    ogo.message('Writing output file:')
    ogo.message('  {}'.format(output_image))
    sitk.WriteImage(ct_final, output_image)            
    
    ogo.message('Done.')
    

def main():
    description = '''
Function performs morphological operations. Although it is possible to 
apply this to images, it is meant for application with labels. The operations
include:

- dilate
- erode
- closing (erosion of the dilation)
- opening (dilation of the erosion)

A list of specific labels can be defined to restrict the operation to those 
labels only. The labels are processed in the order they are defined in the
list.
'''

    epilog = '''
Example call: 
     
ogoMorphologicalOperation image_in.nii.gz image_out.nii.gz --operation erode --operation_labels 1 2 3 --kernels 2 3 4
ogoMorphologicalOperation image_in.nii.gz image_out.nii.gz --operation erode --operation_labels 1 2 3 --kernel 3

'''
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoMorphologicalOperation",
        description=description,
        epilog=epilog
    )

    parser.add_argument('input_image', help='Input image file (*.nii, *.nii.gz)')
    parser.add_argument('output_image', help='Output image file (*.nii, *.nii.gz)')
    parser.add_argument('--operation_labels', type=int, nargs='*', default=[], metavar='LABEL', help='List of labels to apply operation (default: all)')
    parser.add_argument('--operation', default='dilate', choices=['dilate', 'erode', 'opening', 'closing'],
                                                           help='Select morphological operation (default: %(default)s)')
    parser.add_argument('--kernel', type=int, default=1, metavar='KERNEL', help='Default kernel for morphological operation (default: %(default)s)')
    parser.add_argument('--kernels', type=int, nargs='*', metavar='KERNELS', help='List of kernels for each label (overrides --kernel)')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite output file without asking')

    args = parser.parse_args()
    print(echo_arguments('MorphologicalOperation', vars(args)))

    # Run program
    MorphologicalOperation(**vars(args))


if __name__ == '__main__':
    main()