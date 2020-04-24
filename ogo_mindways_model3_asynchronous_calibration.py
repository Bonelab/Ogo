#####
# ogo_asynchronous_calibration.py
#
# This script performs asychroous calibration from an input image and input calibration
# image and mask of the calibration phantom. Density calibration phantom is required to
# be the Mindways Model 3 CT phantom.
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
import ogo_helper as ogo
import os
import sys
import argparse
import time
from datetime import date
from collections import OrderedDict

##
# Start Script
ogo.message("Start of Script...")

##
# Set up the argument parser for the Script
parser = argparse.ArgumentParser(
    description="""This script performs asynchronous density calibration for an input image and input phantom image using an input mask of the density calibration rods. This script is only compatible with the Mindways Model 3 CT phantom. INPUT: Image (as Dicom or NIFTI), Calibration Phantom Image (Dicom or NIFTI), and Calibration Mask (as Dicom or NIFTI) OUTPUT: K2HPO4 Density Calibrated Image (NIFTI), Text file of Calibration Parameters and associated information""")

parser.add_argument("image_input",
    help = "DICOM directory or *.nii image file")
parser.add_argument("phantom_input",
    help = "DICOM directory or *.nii phantom image file")
parser.add_argument("calibration_mask_input",
    help = "DICOM direcory or *.nii mask image of calibration phantom")

parser.add_argument("--phantom_h2o_densities", type = float, nargs = 5,
    default = [1012.25, 1056.95, 1103.57, 1119.52, 923.20],
    help = "Set the H2O equivalent density values for the  Mindways Calibration Phantom. Set Rod A to Rod E in mg/cc. (Default: %(default)s mg/cc)")
parser.add_argument("--phantom_k2hpo4_densities", type = float, nargs = 5,
    default = [-51.83, -53.40, 58.88, 157.05, 375.83],
    help = "Set the K2HPO4 equivalent density values for the  Mindways Calibration Phantom. Set Rod A to Rod E in mg/cc. (Default: %(default)s mg/cc)")

##
# Collect the input arguments
args = parser.parse_args()
image = args.image_input
phantom = args.phantom_input
mask = args.calibration_mask_input
phantom_h2o_densities = args.phantom_h2o_densities
phantom_k2hpo4_densities = args.phantom_k2hpo4_densities

##
# Determine image locations and names of files
image_pathname = os.path.dirname(image)
image_basename = os.path.basename(image)
phantom_pathname = os.path.dirname(phantom)
phantom_basename = os.path.basename(phantom)
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
# Read input image with correct reader
ogo.message("Reading input phantom image...")
if (os.path.isdir(phantom)):
    ogo.message("Input Phantom image is DICOM")
    phantomData = ogo.readDCM(phantom)

else:
    ext = os.path.splitext(phantom)[1]
    if (ext == ".nii" or ext == ".nifti"):
        ogo.message("Input Phantom image is NIFTI")
        phantomData = ogo.readNii(phantom)
    else:
        print(("ERROR: Phantom image format not recognized for " + phantom))
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
rod_A_mask = ogo.applyMask(phantomData, rod_A)
rod_A_mean = ogo.imageHistogramMean(rod_A_mask)

ogo.message("Extracting calibration rod B ROI...")
rod_B = ogo.maskThreshold(maskData, 2)
rod_B_mask = ogo.applyMask(phantomData, rod_B)
rod_B_mean = ogo.imageHistogramMean(rod_B_mask)

ogo.message("Extracting calibration rod C ROI...")
rod_C = ogo.maskThreshold(maskData, 3)
rod_C_mask = ogo.applyMask(phantomData, rod_C)
rod_C_mean = ogo.imageHistogramMean(rod_C_mask)

ogo.message("Extracting calibration rod D ROI...")
rod_D = ogo.maskThreshold(maskData, 4)
rod_D_mask = ogo.applyMask(phantomData, rod_D)
rod_D_mean = ogo.imageHistogramMean(rod_D_mask)

ogo.message("Extracting calibration rod E ROI...")
rod_E = ogo.maskThreshold(maskData, 5)
rod_E_mask = ogo.applyMask(phantomData, rod_E)
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
    fileName = image_basename + "_AC_K2HPO4.nii"
elif imageType == 2:
    fileName = image_basename.replace(".nii", "_AC_K2HPO4.nii")
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
parameters_dict['Phantom Directory'] = phantom_pathname
parameters_dict['Phantom'] = phantom_basename
parameters_dict['Mask Directory'] = mask_pathname
parameters_dict['Mask'] = mask_basename
parameters_dict['Phantom H2O Densities'] = phantom_h2o_densities
parameters_dict['Phantom K2HPO4 Densities'] = phantom_k2hpo4_densities
parameters_dict['Phantom Density Rod HU'] = phantom_HU
parameters_dict['Calibration Slope'] = cali_parameters['Calibration Slope']
parameters_dict['Calibration Y-Intercept'] = cali_parameters['Calibration Y-Intercept']

# Write the output text file
txt_fileName = org_fileName + "_AsyncCalibParameters.txt"
ogo.message("Writing parameters to output text file: %s" % txt_fileName)
ogo.writeTXTfile(parameters_dict, txt_fileName, image_pathname)

##
# End of script
ogo.message("End of Script.")
sys.exit()
