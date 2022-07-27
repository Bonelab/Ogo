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
import numpy as np
from scipy import stats
from datetime import date
from collections import OrderedDict
from vtk.util.numpy_support import vtk_to_numpy

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
def phantom(input_image,input_mask,output_image,async_image,phantom,overwrite,func):
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
    
    # Get the correct reference calibration phantom
    phantom_dict = ogo.get_phantom(phantom)
    #print("\n".join("{:33s} = {}".format(k,v) for k, v in phantom_dict.items()))
    ogo.message('Calibration phantom is \"{}\".'.format(phantom_dict['name']))
    if phantom_dict['type'] in 'CHA':
        cha_phantom = phantom_dict['densities']
        ogo.message("  CHA Phantom [mg/cc] --> ["+" ".join("{:8.3f}".format(i) for i in cha_phantom)+"]")
    else:
        cha_phantom = phantom_dict['densities']
        h2o_phantom = phantom_dict['h2o_densities']
        ogo.message("  K2HPO4 Phantom [mg/cc] --> ["+" ".join("{:8.3f}".format(i) for i in cha_phantom)+"]")
        ogo.message("     H2O Phantom [mg/cc] --> ["+" ".join("{:8.3f}".format(i) for i in h2o_phantom)+"]")
        
    # Gather the image and mask
    if (async_image):
        image_array = vtk_to_numpy(reader_async.GetOutput().GetPointData().GetScalars())
    else:
        image_array = vtk_to_numpy(reader_image.GetOutput().GetPointData().GetScalars())
    mask_array = vtk_to_numpy(reader_mask.GetOutput().GetPointData().GetScalars())
    
    # Calculate the mean and SD of each phantom rod in the image
    ogo.message("Extracting calibration data from image...")
    ogo.message("       {:>8s} {:>8s} {:>8s}".format('Mean','StdDev','#Voxels'))
    rod_labels = phantom_dict['rod_labels']
    rod_names = phantom_dict['rod_names']
    phantom_HU = []

    for rod in range(phantom_dict['number_rods']):
        rod_label = rod_labels[rod]
        if rod_label not in mask_array:
            os.sys.exit('[ERROR] Mask image does not contain label \"{}\". Are you using the wrong mask?'.format(rod_label))
        rod_name = rod_names[rod]
        rod_mean = np.mean(image_array[mask_array == rod_label])
        rod_std = np.std(image_array[mask_array == rod_label])
        rod_count = len(image_array[mask_array == rod_label])
        ogo.message('Rod {:s}: {:8.3f} {:8.3f} {:8d}'.format(rod_name,rod_mean,rod_std,rod_count)) # Report mean, SD, # voxels
        phantom_HU.append(rod_mean)
        
    ogo.message("  Image [HU]--> ["+" ".join("{:8.3f}".format(i) for i in phantom_HU)+"]")

    exit()
    # See Figure 2.1 on page 26 of Andy's thesis
    if (phantom in 'Mindways Model 3'):
        y_values = np.subtract(phantom_HU,h2o_phantom_densities)
        x_values = k2hpo4_phantom_densities
        regression_parameters = stats.linregress(x_values, y_values)
        print(regression_parameters)
        # convert slope and y-intercept to CT parameters
        sigma_ct = regression_parameters[0] - 0.2174
        beta_ct = regression_parameters[1] + 999.6
        calibration_slope = 1/sigma_ct
        calibration_yint = 1*beta_ct/sigma_ct
        print('calibration_slope = {}'.format(calibration_slope))
        print('calibration_yint = {}'.format(calibration_yint))
    elif (phantom in 'B-MAS 200'):
        y_values = 0
    else:
        os.sys.exit('[ERROR] Cannot find appropriate phantom density for \"{}\"'.format(phantom))
    
    exit()
#def phantomParameters(h2o_density, k2hpo4_density, phantom_HU):
#    """Determine the slope and y-intercept for the phantom calibration.
#    The first argument are the phantom specific H2O equivalent density values. The second
#    argument are the phantom specific K2HPO4 equivalent density values. The third
#    argument are the phantom rod mean HU values from the image and mask.
#    Returns the slope and y-intercept for the calibration as a float list.
#    """
#    y_values = np.subtract(phantom_HU, h2o_density)
#    x_values = k2hpo4_density
#    regression_parameters = stats.linregress(x_values, y_values)
#
#    # convert slope and y-intercept to CT parameters
#    sigma_ct = regression_parameters[0] - 0.2174
#    beta_ct = regression_parameters[1] + 999.6
#
#    # Determine calibration parameters
#    calibration_slope = 1/sigma_ct
#    calibration_yint = 1*beta_ct/sigma_ct
#    return {
#    'Calibration Slope':calibration_slope,
#    'Calibration Y-Intercept':calibration_yint
#    }



    ##
    # determine the phantom calibration parameters
    ogo.message("Determining the phantom calibration parameters...")
    cali_parameters = ogo.phantomParameters(phantom_h2o_densities, phantom_k2hpo4_densities, phantom_HU)

    ##
    # Apply calibration parameters the image
    ogo.message("Applying Image Calibration...")
    calibrated_image = ogo.applyPhantomParameters(imageData, cali_parameters)
    
#    if (async_image):
    
    
    
    #Mindways Model 3 CT phantom
    #B-MAS 200 CT phantom
    
    # Get the correct phantom parameters (Mindways, B200)
    # Calculate the phantom values from the CT image (direct or asynchronous)
    # Calculate the linear relation and apply it to the input image
    # Write out some report information
    
    #labels = np.unique(array)
    ## Remove label 0 from list (background)
    #if (labels[0] == 0):
    #    labels = labels[1:]
    
    ## Loop through each of the valid labels and calculate BMD
    #for idx,lab in enumerate(labels):
    #
    #    bone_mask = ogo.maskThreshold(reader_mask.GetOutput(), lab)
    #    bone_VOI = ogo.applyMask(reader_image.GetOutput(), bone_mask)
    #    bmd_outcomes = ogo.bmd_metrics(bone_VOI)
    #    
    #    #parameters_dict['ID'] = os.path.basename(image_filename)
    #    parameters_dict['Script'] = os.path.basename(sys.argv[0])
    #    parameters_dict['Version'] = script_version
    #    parameters_dict['Created'] = str(date.today())
    #    parameters_dict['Image'] = os.path.basename(image_filename)
    #    #parameters_dict['ImageDir'] = os.path.dirname(image_filename)
    #    parameters_dict['Mask'] = os.path.basename(mask_filename)
    #    #parameters_dict['MaskDir'] = os.path.dirname(mask_filename)
    #    parameters_dict['Label'] = lab
    #    parameters_dict['LabelDesc'] = labelsDict[str(labels[idx])]
    #    
    #    parameters_dict['Integral BMD [mg/cc]'] = bmd_outcomes['Integral BMD [mg/cc]']
    #    parameters_dict['Integral BMC [mg]'] = bmd_outcomes['Integral BMC [mg]']
    #    parameters_dict['Bone Volume [mm^3]'] = bmd_outcomes['Bone Volume [mm^3]']
    #    parameters_dict['Bone Volume [cm^3]'] = bmd_outcomes['Bone Volume [cm^3]']
    #
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
  B-MAS 200 CT phantom
  
Valid file formats are NIFTII: .nii, .nii.gz

Outputs are the K2HPO4 density calibrated image and a text file of calibration
parameters and results.

'''

    epilog='''
USAGE: 
ogoImageCalibration internal input_image.nii.gz input_mask.nii.gz output.nii.gz  
ogoImageCalibration phantom  input_image.nii.gz asynch_mask.nii.gz --async_image asynch.nii.gz output.nii.gz  

python ImageCalibration.py phantom \
  /Users/skboyd/Desktop/ML/test/kub.nii.gz \
  /Users/skboyd/Desktop/ML/test/kub_mask.nii.gz \
  /Users/skboyd/Desktop/ML/test/test.nii \
  --phantom 'Mindways Model 3 CT'
    
python ImageCalibration.py phantom \
  /Users/skboyd/Desktop/ML/test/kub.nii.gz \
  /Users/skboyd/Desktop/ML/test/async_mask_mindways.nii.gz \
  /Users/skboyd/Desktop/ML/test/test.nii \
  --async_image /Users/skboyd/Desktop/ML/test/asynch.nii.gz \
  --phantom 'Mindways Model 3 CT' 

python ImageCalibration.py internal \
  /Users/skboyd/Desktop/ML/test/qct.nii \
  /Users/skboyd/Desktop/ML/test/qct_mask.nii.gz \
  /Users/skboyd/Desktop/ML/test/test.nii

Citation:
Michalski AS, Besler BA, Michalak GJ, Boyd SK, 2020. CT-based internal density 
calibration for opportunistic skeletal assessment using abdominal CT scans. 
Med Eng Phys 78, 55-63.
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
    parser_phantom.add_argument('--phantom', default='Mindways Model 3 CT', choices=['Mindways Model 3 CT','Mindways Model 3 QA','QRM-BDC 3-rod','QRM-BDC 6-rod','Image Analysis QCT-3D Plus','B-MAS 200'], help='Specify phantom used (default: %(default)s)')
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