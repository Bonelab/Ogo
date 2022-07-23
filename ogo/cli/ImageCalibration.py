# /------------------------------------------------------------------------------+
# | 19-JUL-2022                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+

#####
# ogo_phantom_calibration.py
#
# This script performs phantom-based calibration from an input image and mask of the
# calibration phantom. Density calibration phantom is required to be the Mindways Model 3
# CT phantom.
#
#####
#
# Andrew Michalski
# University of Calgary
# Biomedical Engineering Graduate Program
# April 23, 2019
# Modified to Py3: March 25, 2020
#####

script_version = 1.0

##
# Import the required modules
import os
import sys
import argparse
import vtk
import time
from datetime import date
from collections import OrderedDict

import ogo.cli.Helper as ogo
from ogo.util.echo_arguments import echo_arguments

##
# Start Script
def mindwaysModel3PhantomCalib(args):
    ogo.message("Start of Script...")

    ##
    # Collect the input arguments
    image = args.image_input
    mask = args.calibration_mask_input
    phantom_h2o_densities = args.phantom_h2o_densities
    phantom_k2hpo4_densities = args.phantom_k2hpo4_densities

    ##
    # Determine image locations and names of files
    image_pathname = os.path.dirname(image)
    image_basename = os.path.basename(image)
    mask_pathname = os.path.dirname(mask)
    mask_basename = os.path.basename(mask)

    ##
    # Read input image with correct reader
    ogo.message("Reading input image...")
    if (os.path.isdir(image)):
        ogo.message("Input image is DICOM")
        imageData = ogo.readDCM(image)
        imageType = 1

    else:
        ext = os.path.splitext(image)[1]
        if (ext == ".nii" or ext == ".nifti"):
            ogo.message("Input image is NIFTI")
            imageData = ogo.readNii(image)
            imageType = 2
        else:
            print(("ERROR: image format not recognized for " + image))
            sys.exit()

    ##
    # Read mask image with correct reader
    ogo.message("Reading input mask image...")
    if (os.path.isdir(mask)):
        ogo.message("Input mask is DICOM")
        maskData = ogo.readDCM(mask)

    else:
        ext = os.path.splitext(mask)[1]
        if (ext == ".nii" or ext == ".nifti"):
            ogo.message("Input mask is NIFTI")
            maskData = ogo.readNii(mask)
        else:
            print(("ERROR: image format not recognized for " + mask))
            sys.exit()

    ##
    # Extract calibration rods from calibration mask image
    ogo.message("Extracting calibration rod A ROI...")
    rod_A = ogo.maskThreshold(maskData, 1)
    rod_A_mask = ogo.applyMask(imageData, rod_A)
    rod_A_mean = ogo.imageHistogramMean(rod_A_mask)

    ogo.message("Extracting calibration rod B ROI...")
    rod_B = ogo.maskThreshold(maskData, 2)
    rod_B_mask = ogo.applyMask(imageData, rod_B)
    rod_B_mean = ogo.imageHistogramMean(rod_B_mask)

    ogo.message("Extracting calibration rod C ROI...")
    rod_C = ogo.maskThreshold(maskData, 3)
    rod_C_mask = ogo.applyMask(imageData, rod_C)
    rod_C_mean = ogo.imageHistogramMean(rod_C_mask)

    ogo.message("Extracting calibration rod D ROI...")
    rod_D = ogo.maskThreshold(maskData, 4)
    rod_D_mask = ogo.applyMask(imageData, rod_D)
    rod_D_mean = ogo.imageHistogramMean(rod_D_mask)

    ogo.message("Extracting calibration rod E ROI...")
    rod_E = ogo.maskThreshold(maskData, 5)
    rod_E_mask = ogo.applyMask(imageData, rod_E)
    rod_E_mean = ogo.imageHistogramMean(rod_E_mask)

    phantom_HU = [rod_A_mean[0], rod_B_mean[0], rod_C_mean[0], rod_D_mean[0], rod_E_mean[0]]

    ##
    # determine the phantom calibration parameters
    ogo.message("Determining the phantom calibration parameters...")
    cali_parameters = ogo.phantomParameters(phantom_h2o_densities, phantom_k2hpo4_densities, phantom_HU)

    ##
    # Apply calibration parameters the image
    ogo.message("Applying Image Calibration...")
    calibrated_image = ogo.applyPhantomParameters(imageData, cali_parameters)

    ##
    # write out calibrated image
    if imageType == 1:
        fileName = image_basename + "_K2HPO4.nii"
    elif imageType == 2:
        fileName = image_basename.replace(".nii", "_PC_K2HPO4.nii")
    else:
        print(("ERROR: image format not recognized for " + image))
        sys.exit()

    ogo.message("Writing out the calibrated image: %s" % fileName)
    ogo.writeNii(calibrated_image, fileName, image_pathname)

    ##
    # Write parameters to output txt file
    if imageType == 1:
        org_fileName = image_basename
    elif imageType == 2:
        org_fileName = image_basename.replace(".nii", "")

    # Define the dictionary for output parameters
    parameters_dict = OrderedDict()
    parameters_dict['ID'] = org_fileName
    parameters_dict['Output File'] = fileName
    parameters_dict['Python Script'] = sys.argv[0]
    parameters_dict['Version'] = script_version
    parameters_dict['Date Created'] = str(date.today())
    parameters_dict['Image Directory'] = image_pathname
    parameters_dict['Image'] = image_basename
    parameters_dict['Mask Directory'] = mask_pathname
    parameters_dict['Mask'] = mask_basename
    parameters_dict['Phantom H2O Densities'] = phantom_h2o_densities
    parameters_dict['Phantom K2HPO4 Densities'] = phantom_k2hpo4_densities
    parameters_dict['Phantom Density Rod HU'] = phantom_HU
    parameters_dict['Calibration Slope'] = cali_parameters['Calibration Slope']
    parameters_dict['Calibration Y-Intercept'] = cali_parameters['Calibration Y-Intercept']

    # Write the output text file
    txt_fileName = org_fileName + "_PhanCalibParameters.txt"
    ogo.message("Writing parameters to output text file: %s" % txt_fileName)
    ogo.writeTXTfile(parameters_dict, txt_fileName, image_pathname)

    ##
    # End of script
    ogo.message("End of Script.")
    sys.exit()

# PHANTOM CALIBARATION ------------------------------------------------------------------
def phantom(input_image,input_mask,output_image,async_image,overwrite,func):
    ogo.message('Starting phantom based calibration.')
    
    # Check if output exists and should overwrite
    if os.path.isfile(output_image) and not overwrite:
        result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(output_image))
        if result.lower() not in ['y', 'yes']:
            ogo.message('Not overwriting. Exiting...')
            os.sys.exit()

    # Read input image
    if not os.path.isfile(input_image):
        os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(input_image))

    if input_image.lower().endswith('.nii'):
        reader_image = vtk.vtkNIFTIImageReader()
    elif input_image.lower().endswith('.nii.gz'):
        reader_image = vtk.vtkNIFTIImageReader()
    else:
        os.sys.exit('[ERROR] Cannot find reader for file \"{}\"'.format(input_image))

    ogo.message('Reading input image to be calibrated...')
    ogo.message('      \"{}\"'.format(input_image))
    reader_image.SetFileName(input_image)
    reader_image.Update()

    # Read input mask
    if not os.path.isfile(input_mask):
        os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(input_mask))

    if input_mask.lower().endswith('.nii'):
        reader_mask = vtk.vtkNIFTIImageReader()
    elif input_mask.lower().endswith('.nii.gz'):
        reader_mask = vtk.vtkNIFTIImageReader()
    else:
        os.sys.exit('[ERROR] Cannot find reader for file \"{}\"'.format(input_mask))

    ogo.message('Reading input mask used for calibration...')
    ogo.message('      \"{}\"'.format(input_mask))
    reader_mask.SetFileName(input_mask)
    reader_mask.Update()
    
    # Read asynchronous input image, if available
    if (async_image):
        if not os.path.isfile(async_image):
            os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(async_image))
        if async_image.lower().endswith('.nii'):
            reader_async = vtk.vtkNIFTIImageReader()
        elif async_image.lower().endswith('.nii.gz'):
            reader_async = vtk.vtkNIFTIImageReader()
        else:
            os.sys.exit('[ERROR] Cannot find reader for file \"{}\"'.format(async_image))
        
        ogo.message('Reading asynchronous image of phantom...')
        ogo.message('      \"{}\"'.format(async_image))
        reader_async.SetFileName(async_image)
        reader_async.Update()

    if (async_image):
        ogo.message('NOTE: Asynchronous calibration selected.')
        ogo.message('      Using mask file against asynchronous file:')
        ogo.message('      \"{}\"'.format(input_mask))
        ogo.message('      \"{}\"'.format(async_image))
    
    
    ogo.message('Done phantom based calibration.')
    
# INTERNAL CALIBARATION ------------------------------------------------------------------
def internal(input_image,input_mask,output_image,overwrite,func):
    ogo.message('Starting internal calibration.')

    # Check if output exists and should overwrite
    if os.path.isfile(output_image) and not overwrite:
        result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(output_image))
        if result.lower() not in ['y', 'yes']:
            ogo.message('Not overwriting. Exiting...')
            os.sys.exit()

    # Read input image
    if not os.path.isfile(input_image):
        os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(input_image))

    if input_image.lower().endswith('.nii'):
        reader_image = vtk.vtkNIFTIImageReader()
    elif input_image.lower().endswith('.nii.gz'):
        reader_image = vtk.vtkNIFTIImageReader()
    else:
        os.sys.exit('[ERROR] Cannot find reader for file \"{}\"'.format(input_image))

    ogo.message('Reading input image to be calibrated...')
    ogo.message('      \"{}\"'.format(input_image))
    reader_image.SetFileName(input_image)
    reader_image.Update()

    # Read input mask
    if not os.path.isfile(input_mask):
        os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(input_mask))

    if input_mask.lower().endswith('.nii'):
        reader_mask = vtk.vtkNIFTIImageReader()
    elif input_mask.lower().endswith('.nii.gz'):
        reader_mask = vtk.vtkNIFTIImageReader()
    else:
        os.sys.exit('[ERROR] Cannot find reader for file \"{}\"'.format(input_mask))

    ogo.message('Reading input mask used for calibration...')
    ogo.message('      \"{}\"'.format(input_mask))
    reader_mask.SetFileName(input_mask)
    reader_mask.Update()

    ogo.message('Done internal calibration.')
    
def main():
    # Setup description
    description='''
Density calibration of clinical CT by either use of a phantom in the scan, 
asynchronous phantom or internal calibration.

Current phantoms are:
  Mindways Model 3 CT phantom
  B-MAS200
  
Valid file formats are NIFTII: .nii, .nii.gz

Outputs are the K2HPO4 density calibrated image and a text file of calibration
parameters and results.

'''

    epilog='''
USAGE: 
ogoImageCalibration internal input.nii.gz output.nii.gz  
ogoImageCalibration phantom  input.nii.gz output.nii.gz  
ogoImageCalibration phantom  input.nii.gz output.nii.gz --async_image asynch.nii.gz

python ImageCalibration.py phantom \
  /Users/skboyd/Desktop/ML/test/qct.nii \
  /Users/skboyd/Desktop/ML/test/qct_mask.nii.gz \
  /Users/skboyd/Desktop/ML/test/test.nii
  
python ImageCalibration.py phantom \
  /Users/skboyd/Desktop/ML/test/qct.nii \
  /Users/skboyd/Desktop/ML/test/qct_mask.nii.gz \
  /Users/skboyd/Desktop/ML/test/test.nii \
  --async_image /Users/skboyd/Desktop/ML/test/async.nii.gz

python ImageCalibration.py internal \
  /Users/skboyd/Desktop/ML/test/qct.nii \
  /Users/skboyd/Desktop/ML/test/qct_mask.nii.gz \
  /Users/skboyd/Desktop/ML/test/test.nii

'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoImageCalibration",
        description=description,
        epilog=epilog
    )
    subparsers = parser.add_subparsers()
    
    # Phantom
    parser_phantom = subparsers.add_parser('phantom')
    parser_phantom.add_argument('input_image', help='Input image file (*.nii, *.nii.gz)')
    parser_phantom.add_argument('input_mask', help='Input image mask for either input image or asynchronous image, if present (*.nii, *.nii.gz)')
    parser_phantom.add_argument('output_image', help='Output image file (*.nii, *.nii.gz)')
    parser_phantom.add_argument('--async_image', default='', metavar='IMAGE', help='Asynchronous image of phantom (*.nii, *.nii.gz)')
    parser_phantom.add_argument('--overwrite', action='store_true', help='Overwrite output without asking')
    parser_phantom.set_defaults(func=phantom)

    # Internal
    parser_internal = subparsers.add_parser('internal')
    parser_internal.add_argument('input_image', help='Input image file (*.nii, *.nii.gz)')
    parser_internal.add_argument('input_mask', help='Input image mask file (*.nii, *.nii.gz)')
    parser_internal.add_argument('output_image', help='Output image file (*.nii, *.nii.gz)')
    parser_internal.add_argument('--overwrite', action='store_true', help='Overwrite output without asking')
    parser_internal.set_defaults(func=internal)
    
    #parser.add_argument("image_input",
    #    help = "DICOM directory or *.nii image file")
    #parser.add_argument("phantom_input",
    #    help = "DICOM directory or *.nii phantom image file")
    #parser.add_argument("calibration_mask_input",
    #    help = "DICOM direcory or *.nii mask image of calibration phantom")
    #
    #parser.add_argument("--phantom_h2o_densities", type = float, nargs = 5,
    #    default = [1012.25, 1056.95, 1103.57, 1119.52, 923.20],
    #    help = "Set the H2O equivalent density values for the  Mindways Calibration Phantom. Set Rod A to Rod E in mg/cc. (Default: %(default)s mg/cc)")
    #parser.add_argument("--phantom_k2hpo4_densities", type = float, nargs = 5,
    #    default = [-51.83, -53.40, 58.88, 157.05, 375.83],
    #    help = "Set the K2HPO4 equivalent density values for the  Mindways Calibration Phantom. Set Rod A to Rod E in mg/cc. (Default: %(default)s mg/cc)")

    # Parse and display
    args = parser.parse_args()
    print(echo_arguments('ImageCalibration', vars(args)))
    
    # Run program
    args.func(**vars(args))

if __name__ == '__main__':
    main()