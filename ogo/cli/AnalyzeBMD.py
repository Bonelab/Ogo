# /------------------------------------------------------------------------------+
# | 19-JUL-2022                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+
# Import the required modules
import os
import sys
import vtk
import numpy as np
import argparse
import time
from datetime import date
from collections import OrderedDict
from ogo.util.echo_arguments import echo_arguments
import ogo.util.Helper as ogo
import ogo.dat.OgoMasterLabels as lb
from vtk.util.numpy_support import vtk_to_numpy

script_version = 2.0


# Start Script
def AnalyzeBMD(image_filename, mask_filename, labels, output_filename, custom_labels=None, noheader=False, overwrite=False):
    if output_filename:
        ogo.message("Start of Script...")

    if custom_labels:
      [labelsDict,labelsHdr] = ogo.custom_labelsDict(custom_labels)
    else:
      # Master list of labels
      labelsDict = lb.labels_dict

    # Check if output exists and should overwrite
    ogo.pass_check_if_output_exists(output_filename,overwrite)

    # Set up to read image and mask inputs
    ogo.pass_check_if_file_exists(image_filename)

    # Check endings
    ogo.pass_check_file_ending(image_filename)
    ogo.pass_check_file_ending(mask_filename)
    
    # Read inputs
    if output_filename:
        ogo.message("Reading input calibrated image...")
    reader_image = vtk.vtkNIFTIImageReader()
    reader_image.SetFileName(image_filename)
    reader_image.Update()
    if output_filename:
        ogo.message("Reading input mask image...")
    reader_mask = vtk.vtkNIFTIImageReader()
    reader_mask.SetFileName(mask_filename)
    reader_mask.Update()

    array = vtk_to_numpy(reader_mask.GetOutput().GetPointData().GetScalars())

    # Check that at least one label is defined
    if len(labels) < 1:
        print('ERROR: At least one label must be defined.')
        os.sys.exit()

    # If user does not define labels, then select all labels available
    if labels[0] == 0:
        if output_filename:
            ogo.message('No labels specified. Building list...')
        labels = np.unique(array).tolist()
        print(labels)
        #labels = np.unique(array)

    # Remove label 0 from list (background)
    if labels[0] == 0:
        labels = labels[1:]

    # Check that labels are in range
    if output_filename:
        ogo.message('Labels to be applied:')
    for idx, ids in enumerate(labels):
        if labels[idx] < 0 or labels[idx] > 255:
            print('ERROR: Labels out of range. Must be 0-255.')
            print('       Invalid label is {:d}.'.format(labels[idx]))
            os.sys.exit()
        #        if (str(labels[idx]) not in labelsDict):
        if labels[idx] not in labelsDict:
            print("ERROR: Label {:d} does not exist in list of labels.".format(labels[idx]))
            os.sys.exit()
        if output_filename:
            ogo.message('!> label {:3d} {:12s} --> {:d} voxels'.format(labels[idx], labelsDict[labels[idx]]['LABEL'],
                                                                       np.count_nonzero(array == labels[idx])))

    parameters_dict = OrderedDict()

    if output_filename:
        ogo.message('Writing output to {:s}'.format(output_filename))
        txt_file = open(output_filename, "w")

    # Loop through each of the valid labels and calculate BMD
    for idx, lab in enumerate(labels):

        bone_mask = ogo.maskThreshold(reader_mask.GetOutput(), lab)
        bone_VOI = ogo.applyMask(reader_image.GetOutput(), bone_mask)
        bmd_outcomes = ogo.bmd_metrics(bone_VOI)

        # parameters_dict['ID'] = os.path.basename(image_filename)
        parameters_dict['Script'] = os.path.basename(sys.argv[0])
        parameters_dict['Version'] = script_version
        parameters_dict['Created'] = str(date.today())
        parameters_dict['Image'] = os.path.basename(image_filename)
        # parameters_dict['ImageDir'] = os.path.dirname(image_filename)
        parameters_dict['Mask'] = os.path.basename(mask_filename)
        # parameters_dict['MaskDir'] = os.path.dirname(mask_filename)
        parameters_dict['Label'] = lab
        parameters_dict['LabelDesc'] = labelsDict[lab]['LABEL']

        parameters_dict['Integral BMD [mg/cc]'] = '{:.3f}'.format(bmd_outcomes['Integral BMD [mg/cc]'])
        parameters_dict['Integral BMC [mg]'] = '{:.3f}'.format(bmd_outcomes['Integral BMC [mg]'])
        parameters_dict['Bone Volume [mm^3]'] = '{:.3f}'.format(bmd_outcomes['Bone Volume [mm^3]'])
        parameters_dict['Bone Volume [cm^3]'] = '{:.3f}'.format(bmd_outcomes['Bone Volume [cm^3]'])

        if output_filename:
            if idx == 0 and not noheader:
                txt_file.write(",".join("{}".format(k) for k, v in parameters_dict.items()))
                txt_file.write("\n")
            txt_file.write(",".join("{}".format(v) for k, v in parameters_dict.items()))
            txt_file.write("\n")

            print("!> ---------------------------------------------------- RECORD")
            print("\n".join("{} = {}".format(k, v) for k, v in parameters_dict.items()))
        else:
            if idx == 0 and not noheader:
                print(",".join("{}".format(k) for k, v in parameters_dict.items()))
            print(",".join("{}".format(v) for k, v in parameters_dict.items()))

    if output_filename:
        txt_file.close()
        print()
        ogo.message("End of Script.")

    sys.exit()


def main():
    description = '''
This script computes BMD by applying a mask to the calibrated input image. 
Multiple labels representing different bones in the mask image can be assessed
at once. However, if only a subset of labels are needed to be analysed they
should be defined (see example calls).

If no output text file is defined then results are output to the screen.

WARNING: If your input image is not calibrated then the results here will be 
incorrect. There is no calibration done as part of this application.

'''
    epilog = '''
Example calls:

ogoAnalyzeBMD image_k2hpo4.nii mask.nii.gz
ogoAnalyzeBMD image_k2hpo4.nii mask.nii.gz >> results.txt
ogoAnalyzeBMD image_k2hpo4.nii mask.nii.gz --output_filename output.txt
ogoAnalyzeBMD image_k2hpo4.nii mask.nii.gz --labels 7 8 9 10
'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoAnalyzeBMD",
        description=description,
        epilog=epilog
    )

    parser.add_argument('image_filename', help='Input image file (*.nii, *.nii.gz)')
    parser.add_argument('mask_filename', help='Input image mask file (*.nii, *.nii.gz)')
    parser.add_argument('--labels', type=int, nargs='*', default=[0], metavar='ID',
                        help='Specify list of labels; space separated (e.g. 1 2 3) (Default: report all available)')
    parser.add_argument('--output_filename', default=None, metavar='TEXTFILE',
                        help='Output file name (*.txt) (Default: dump to screen)')
    parser.add_argument('--custom_labels', default=None, metavar='TEXTFILE',
                        help='Custom labels of type ITKSnap (*.txt) (Default: OgoMasterLabels)')
    parser.add_argument('--noheader', action='store_true', help='Suppress printing header (default: %(default)s)')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite output without asking')
    print()

    # Parse and display
    args = parser.parse_args()

    if args.output_filename:
        print(echo_arguments('AnalyzeBMD', vars(args)))

    # Run program
    AnalyzeBMD(**vars(args))


if __name__ == '__main__':
    main()
