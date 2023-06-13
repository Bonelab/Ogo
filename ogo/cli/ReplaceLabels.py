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
import ogo.dat.OgoMasterLabels as lb

import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy
from ogo.util.echo_arguments import echo_arguments
import ogo.util.Helper as ogo

def get_labels(ct):
    filt = sitk.LabelShapeStatisticsImageFilter()
    filt.Execute(ct)
    labels = filt.GetLabels()
    return labels

def get_label_description(label):
    try:
        desc = lb.labels_dict[label]['LABEL']
    except KeyError:
        desc = 'unknown label'
    return desc
    
def ReplaceLabels(input_filename, output_filename, inputLabels, outputLabels, keepOnlyLabels, overwrite=False):

    #ogo.message('Start ogoReplaceLabels!')
    
    # Check if output exists and should overwrite
    if os.path.isfile(output_filename) and not overwrite:
        result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(output_filename))
        if result.lower() not in ['y', 'yes']:
            ogo.message('Not overwriting. Exiting...')
            os.sys.exit()

    # Check that length of input labels equals length of output labels
    if len(inputLabels) is not len(outputLabels):
        ogo.message('ERROR: Number of input labels defined must equal number of output labels.')
        ogo.message('       [input #{:d} != output #{:d}]'.format(len(inputLabels),len(outputLabels)))
        os.sys.exit()
    
    # Check that either inputLabels or keepOnlyLabels are defined (one or the other, or both is OK)
    if (not inputLabels) and (not keepOnlyLabels):
        ogo.message('ERROR: No labels defined.')
        ogo.message('       Define --keepOnlyLabels or --inputLabels & --outputLabels')
        os.sys.exit()

    # Feedback on labels to be replaced
    ogo.message('')
    ogo.message('Labels to be replaced:')
    if inputLabels:
        for idx,ids in enumerate(inputLabels):
            desc1 = get_label_description(inputLabels[idx])
            desc2 = get_label_description(outputLabels[idx])
            ogo.message('!> {:3d} ({:s})  --> {:3d} ({:s})'.format(inputLabels[idx],desc1,outputLabels[idx],desc2))
            if (inputLabels[idx]<0 or inputLabels[idx]>255 or outputLabels[idx]<0 or outputLabels[idx]>255):
                ogo.message('ERROR: Labels out of range. Must be 0-255.')
                ogo.message('       [Suspects are input label {:d} or output label {:d}]'.format(inputLabels[idx],outputLabels[idx]))
                os.sys.exit()
    else:
        ogo.message('!> none defined')
    map_for_labels = dict(zip(inputLabels, outputLabels))
    
    # Feedback on labels to keep
    ogo.message('')
    ogo.message('Labels to keep:')
    if keepOnlyLabels:
        for idx,ids in enumerate(keepOnlyLabels):
            ogo.message('!> {:3d} ({:s})'.format(keepOnlyLabels[idx],get_label_description(keepOnlyLabels[idx])))
            if (keepOnlyLabels[idx]<0 or keepOnlyLabels[idx]>255):
                ogo.message('ERROR: Labels out of range. Must be 0-255.')
                ogo.message('       [Culprit is {:d}]'.format(keepOnlyLabels[idx]))
                os.sys.exit()
    else:
        ogo.message('!> none defined')
    ogo.message('')
    
    # Read input
    if not os.path.isfile(input_filename):
        os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(input_filename))

    if not (input_filename.lower().endswith('.nii') or input_filename.lower().endswith('.nii.gz')):
        os.sys.exit('[ERROR] Input must be type NIFTI file: \"{}\"'.format(input_filename))
        
    ogo.message('Reading image: ')
    ogo.message('\"{}\"'.format(input_filename))
    ct = sitk.ReadImage(input_filename, sitk.sitkUInt8)
    ogo.message('')
    
    # Get all the available labels in the input image
    labels = get_labels(ct)
    n_labels = len(labels)

    ogo.message('Input image contains the following labels:')
    for label in labels:
        ogo.message('!> {:3d} ({:s})'.format(label,get_label_description(label)))
    
    # Create a copy
    #ct_base = sitk.Image(ct) # Do we need to make a deep copy
    ct_base = ct<0
    
    # Remove all labels except the ones defined to keep
    ogo.message('')
    ogo.message('!> Keeping labels...')
    if keepOnlyLabels:
        for label in labels:
            if label in keepOnlyLabels:
                
                bin_part = ct==label
                mask = 1 - bin_part
                ct_base = sitk.Mask(ct_base, mask)
                ct_base = ct_base + label*bin_part
                
                ogo.message('!>   Passing label {} ({}) to output.'.format(label,get_label_description(label)))
                
    else: # pass all labels through if keepOnlyLabels not set
        ct_base = ct
        ogo.message('!>   Passing all labels to output.')
        
    # Replace labels
    ogo.message('')
    ogo.message('!> Changing labels...')
    if inputLabels:
        for idx,label in enumerate(inputLabels):
            new_label = outputLabels[idx]
            #print(idx,label,new_label)
            if label in get_labels(ct_base):
                bin_part = ct==label
                mask = 1 - bin_part
                ct_base = sitk.Mask(ct_base, mask)
                ct_base = ct_base + new_label*bin_part
                ogo.message('!>   Changing label {} ({}) to {} ({}) in output.'.format(label,get_label_description(label),new_label,get_label_description(new_label)))
            else:
                ogo.message('!>   Label {} ({}) not available.'.format(label,get_label_description(label)))

    # Write output
    labels = get_labels(ct_base)
    ogo.message('')
    ogo.message('Output image contains the following labels:')
    for label in labels:
        ogo.message('!> {:3d} ({:s})'.format(label,get_label_description(label)))
    
    ogo.message('')
    ogo.message('Writing output image:')
    ogo.message('  {}'.format(output_filename))
    sitk.WriteImage(ct_base, output_filename)            
    
def main():
    # Setup description
    description='''
Utility to read a segmented image and replace existing labels with new labels.

It can be used to remove all but the labels you want to keep (keepOnlyLabels)
or to change a label (inputLabels/outputLabels). If changing labels then the
number of inputLabels must equal the number of outputLabels.

You can do both keep labels and change labels at the same time. If the labels
you want to change are NOT in the list of labels you keep (default is to keep
all) then it cannot be swapped. 

Input labels must be between 0 and 255.
'''
    epilog='''
Example call: 
ogoReplaceLabels input.nii.gz output.nii.gz --inputLabels 4 5 --outputLabels 0 7
ogoReplaceLabels input.nii.gz output.nii.gz --inputLabels 4 5 --outputLabels 0 7 --keepOnlyLabels 1 2 3 4 5
ogoReplaceLabels input.nii.gz output.nii.gz --keepOnlyLabels 1 2 3 4 5 6 7
'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoReplaceLabels",
        description=description,
        epilog=epilog
    )
    parser.add_argument('input_filename', help='Input image file (*.nii, *.nii.gz)')
    parser.add_argument('output_filename', help='Output image file (*.nii, *.nii.gz)')
    parser.add_argument('--inputLabels', type=int, nargs='*', default=[], metavar='ID', help='Target label input IDs; space separated (e.g. 1 2 3)')
    parser.add_argument('--outputLabels', type=int, nargs='*', default=[], metavar='ID', help='Target label output IDs; space separated (e.g. 4 5 6)')
    parser.add_argument('--keepOnlyLabels', type=int, nargs='*', default=[], metavar='ID', help='Labels to keep (others removed); (e.g. 7 8 9)')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite output without asking')

    print()

    # Parse and display
    args = parser.parse_args()
    print(echo_arguments('ReplaceLabels', vars(args)))

    # Run program
    ReplaceLabels(**vars(args))

if __name__ == '__main__':
    main()
