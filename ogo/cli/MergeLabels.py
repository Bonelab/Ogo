# /------------------------------------------------------------------------------+
# | 27-NOV-2022                                                                  |
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
import ogo.dat.OgoMasterLabels as lb

from ogo.util.echo_arguments import echo_arguments
import ogo.util.Helper as ogo
from vtk.util.numpy_support import vtk_to_numpy

def merge_labels(input1_filename, input2_filename, output_filename, merge_method, overwrite=False):
    
    # Check if output exists and should overwrite
    if os.path.isfile(output_filename) and not overwrite:
        result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(output_filename))
        if result.lower() not in ['y', 'yes']:
            ogo.message('Not overwriting. Exiting...')
            os.sys.exit()

    # Read inputs 1 and 2
    if not os.path.isfile(input1_filename):
        os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(input1_filename))

    if not (input1_filename.lower().endswith('.nii') or input1_filename.lower().endswith('.nii.gz')):
        os.sys.exit('[ERROR] Input must be type NIFTI file: \"{}\"'.format(input1_filename))
    
    if not os.path.isfile(input2_filename):
        os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(input2_filename))

    if not (input2_filename.lower().endswith('.nii') or input2_filename.lower().endswith('.nii.gz')):
        os.sys.exit('[ERROR] Input must be type NIFTI file: \"{}\"'.format(input2_filename))    
    
    ogo.message('Reading image 1: ')
    ogo.message('  {}'.format(input1_filename))
    ct1 = sitk.ReadImage(input1_filename)

    ogo.message('Reading image 2: ')
    ogo.message('  {}'.format(input2_filename))
    ct2 = sitk.ReadImage(input2_filename)

    merge = sitk.MergeLabelMapFilter()
    print(method)
    match method.split(): # 'Aggregate', 'Strict', 'Keep', 'Pack'
        case ['Aggregate']:
            sitk.MergeLabelMapFilter.Aggregate
        case ['Strict']:
            sitk.MergeLabelMapFilter.Strict
        case ['Keep']:
            sitk.MergeLabelMapFilter.Keep
        case ['Pack']:
            sitk.MergeLabelMapFilter.Pack
        case _:
            ogo.message('unknown')
    
    
    # >>> def file_handler_v1(command):
    # ...     match command.split():
    # ...         case ['show']:
    # ...             print('List all files and directories: ')
    # ...             # code to list files
    # ...         case ['remove', *files]:
    # ...             print('Removing files: {}'.format(files))
    # ...             # code to remove files
    # ...         case _:
    # ...             print('Command not recognized')
    
    #merge.SetMethod(sitk.MergeLabelMapFilter.Aggregate)
    #merge.SetMethod(1)
    #merge.Execute(ct1,ct2)
    
    #print(merge.GetMethod())
    #merger->SetInput(0, objectByObjectLabelMapFilter->GetOutput(0));
    #merger->SetInput(1, selector->GetOutput(1));
    #merger->SetMethod(itk::MergeLabelMapFilterEnums::ChoiceMethod::KEEP);
    
    
    #merger->SetInput(0, objectByObjectLabelMapFilter->GetOutput(0));
    #merger->SetInput(1, selector->GetOutput(1));
    #merger->SetMethod(itk::MergeLabelMapFilterEnums::ChoiceMethod::KEEP);
    
    # thres = sitk.BinaryThresholdImageFilter() # Used to isolate labels and find # of parts
    # thres.SetInsideValue(1)
    # thres.SetOutsideValue(0)
    # thres.SetUpperThreshold(label)
    # thres.SetLowerThreshold(label)
    # ct_thres = thres.Execute(ct)

    #if (debug): print(im)
    #print('-------------------------------------------------------------------------------')
    #print('{:16s} = {:16s}'.format('input_filename',input_filename))
    #print('{:16s} = {:16.10f} {:16.10f} {:16.10f}'.format('spacing',im.GetSpacing()[0],im.GetSpacing()[1],im.GetSpacing()[2]))
    #print('{:16s} = {:16.10f} {:16.10f} {:16.10f}'.format('origin',im.GetOrigin()[0],im.GetOrigin()[1],im.GetOrigin()[2]))
    #print('{:16s} = {:16.10f} {:16.10f} {:16.10f}'.format('direction',im.GetDirection()[0],im.GetDirection()[1],im.GetDirection()[2]))
    #print('{:16s}   {:16.10f} {:16.10f} {:16.10f}'.format('',im.GetDirection()[3],im.GetDirection()[4],im.GetDirection()[5]))
    #print('{:16s}   {:16.10f} {:16.10f} {:16.10f}'.format('',im.GetDirection()[6],im.GetDirection()[7],im.GetDirection()[8]))
    #
    ## Create output image
    #nda = sitk.GetArrayFromImage(im)
    #
    #im_new = sitk.GetImageFromArray(nda)
    #im_new.SetSpacing(im.GetSpacing())
    #im_new.SetOrigin(im.GetOrigin())
    #im_new.SetDirection(im.GetDirection())
    #
    #if (debug): print(im_new)
    #
    #print('-------------------------------------------------------------------------------')
    #print('{:16s} = {:16s}'.format('output_filename',output_filename))
    #print('{:16s} = {:16.10f} {:16.10f} {:16.10f}'.format('spacing',im_new.GetSpacing()[0],im_new.GetSpacing()[1],im_new.GetSpacing()[2]))
    #print('{:16s} = {:16.10f} {:16.10f} {:16.10f}'.format('origin',im_new.GetOrigin()[0],im_new.GetOrigin()[1],im_new.GetOrigin()[2]))
    #print('{:16s} = {:16.10f} {:16.10f} {:16.10f}'.format('direction',im_new.GetDirection()[0],im_new.GetDirection()[1],im_new.GetDirection()[2]))
    #print('{:16s}   {:16.10f} {:16.10f} {:16.10f}'.format('',im_new.GetDirection()[3],im_new.GetDirection()[4],im_new.GetDirection()[5]))
    #print('{:16s}   {:16.10f} {:16.10f} {:16.10f}'.format('',im_new.GetDirection()[6],im_new.GetDirection()[7],im_new.GetDirection()[8]))
    #  
    #sitk.WriteImage(im_new, output_filename)
    
    ogo.message('Done.')

def main():
    # Setup description
    description='''
A utility to merge labels and ensuring correct labels
for each component. 

This is a useful utility if using TotalSegmentator
https://arxiv.org/abs/2208.05868

Merge methods are:
    Aggregate: If the same label is found several times in the label 
               maps, the label objects with the same label are merged.
    Strict:    Keeps the labels unchanged and raises an exception if 
               the same label is found in several images.
    Keep:      Does its best to keep the label unchanged, but if a label
               is already used in a previous label map, a new label is 
               assigned.
    Pack:      Relabel all the label objects by order of processing.

'''
    epilog = '''
Example calls: 
  ogoMergeLabels femur_right.nii.gz femur_left.nii.gz \
                 femurs.nii.gz
                   
python /Users/skboyd/Desktop/ML/code/Ogo/ogo/cli/MergeLabels.py \
/Users/skboyd/Desktop/ML/test/CTDXAICI_0053_V00/femur_left.nii.gz \
/Users/skboyd/Desktop/ML/test/CTDXAICI_0053_V00/femur_right.nii.gz \
/Users/skboyd/Desktop/combined.nii.gz \
--overwrite

'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="MergeLabels",
        description=description,
        epilog=epilog
    )
    parser.add_argument('input1_filename', metavar='FILE IN', help='Input 1 NIFTI file')
    parser.add_argument('input2_filename', metavar='FILE IN', help='Input 2 NIFTI file')
    parser.add_argument('output_filename', metavar='FILE OUT', help='Output NIFTI file (default: %(default)s)')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite output without asking')
    parser.add_argument('--merge_method', default='Strict',choices=['Aggregate', 'Strict', 'Keep', 'Pack'],
                                                           help='Select merge type (default: %(default)s)')

    # Parse and display
    args = parser.parse_args()
    print(echo_arguments('MergeLabels', vars(args)))
        
    # Run program
    merge_labels(**vars(args))

if __name__ == '__main__':
    main()
