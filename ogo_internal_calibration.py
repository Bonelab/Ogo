#####
# ogo_internal_calibration.py
#
# This script performs internal calibration from an input image and tissue reference mask.
# This method has been published as Michalski et al. (2020) "CT-based internal density calibration for opportunistic skeletal assessment using abdominal CT scans" Med Eng Phys
# DOI: https://doi.org/10.1016/j.medengphy.2020.01.009
#####
#
# Andrew Michalski
# University of Calgary
# Biomedical Engineering Graduate Program
# April 25, 2019
# Modified to Py3: March 25, 2020
#####

script_version = 1.0

##
# Import the required modules
import ogo_helper as ogo
import MassAttenuationTables as mat
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
    description="""This script performs internal density calibration for an input image using an input calibration mask of the reference tissues. The reference tissues required in the mask with values are 1 = Adipose, 2 = Air, 3 = Blood, 4 = Cortical Bone, 5 = Skeletal Muscle. INPUT: Image (as Dicom or NIFTI), Calibration Mask (as Dicom or NIFTI) OUTPUT: K2HPO4 Density Calibrated Image (NIFTI), Text file of Calibration Parameters and associated information""")

parser.add_argument("image_input",
    help = "DICOM directory or *.nii image file")
parser.add_argument("calibration_mask_input",
    help = "DICOM direcory or *.nii mask image of calibration phantom")

##
# Collect the input arguments
args = parser.parse_args()
image = args.image_input
mask = args.calibration_mask_input

##
# Determine image locations and names of files
image_pathname = os.path.dirname(image)
image_basename = os.path.basename(image)
mask_pathname = os.path.dirname(mask)
mask_basename = os.path.basename(mask)
org_fileName = image_basename.replace(".nii","")
fileName = image_basename.replace(".nii","_IC_K2HPO4.nii")

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
# Extract reference tissues from the mask
ogo.message("Extracting reference tissue: Adipose...")
adipose_roi = ogo.maskThreshold(maskData, 1)
adipose_mask = ogo.applyMask(imageData, adipose_roi)
adipose_HU = ogo.imageHistogramMean(adipose_mask)
# ogo.message("Adipose ROI Mean HU: %8.4f " % adipose_HU[0])

ogo.message("Extracting reference tissue: Air...")
air_roi = ogo.maskThreshold(maskData, 2)
air_mask = ogo.applyMask(imageData, air_roi)
air_HU = ogo.imageHistogramMean(air_mask)
# ogo.message("Air ROI Mean HU: %8.4f " % air_HU[0])

ogo.message("Extracting reference tissue: Blood...")
blood_roi = ogo.maskThreshold(maskData, 3)
blood_mask = ogo.applyMask(imageData, blood_roi)
blood_HU = ogo.imageHistogramMean(blood_mask)
# ogo.message("Blood ROI Mean HU: %8.4f " % blood_HU[0])

ogo.message("Extracting reference tissue: Cortical Bone...")
bone_roi = ogo.maskThreshold(maskData, 4)
bone_mask = ogo.applyMask(imageData, bone_roi)
bone_HU = ogo.imageHistogramMean(bone_mask)
# ogo.message("Cortical Bone ROI Mean HU: %8.4f " % bone_HU[0])

ogo.message("Extracting reference tissue: Skeletal Muscle...")
muscle_roi = ogo.maskThreshold(maskData, 5)
muscle_mask = ogo.applyMask(imageData, muscle_roi)
muscle_HU = ogo.imageHistogramMean(muscle_mask)
# ogo.message("Skeletal Muscle ROI Mean HU: %8.4f " % muscle_HU[0])

mean_hu = [adipose_HU[0], air_HU[0], blood_HU[0], bone_HU[0], muscle_HU[0]]

##
# Prep Reference Material Tables with interpolation over energy levels 1-200 keV
ogo.message("Deriving material tables...")
adipose_interp = ogo.icInterpolation(mat.adipose_table)
air_interp = ogo.icInterpolation(mat.air_table)
blood_interp = ogo.icInterpolation(mat.blood_table)
bone_interp = ogo.icInterpolation(mat.bone_table)
muscle_interp = ogo.icInterpolation(mat.muscle_table)
k2hpo4_interp = ogo.icInterpolation(mat.k2hpo4_table)
cha_interp = ogo.icInterpolation(mat.cha_table)
triglyceride_interp = ogo.icInterpolation(mat.triglyceride_table)
water_interp = ogo.icInterpolation(mat.water_table)

##
# Determine scan effective energy
ogo.message("Determining the scan effective energy...")
ic_parameters = ogo.icEffectiveEnergy(mean_hu, adipose_interp, air_interp, blood_interp, bone_interp, muscle_interp, k2hpo4_interp, cha_interp, triglyceride_interp, water_interp)
attenuation_values = [ic_parameters['Adipose u/p'], ic_parameters['Air u/p'], ic_parameters['Blood u/p'], ic_parameters['Cortical Bone u/p'], ic_parameters['Skeletal Muscle u/p']]

##
# Determine the HU-Mass Attenuation Relationship
ogo.message("Determining the HU-Mass Attenuation Relationship...")
hu_MUrho = ogo.icLinearRegression(mean_hu, attenuation_values, 'HU-u/p Slope', 'HU-u/p Y-Intercept')

##
# Determine the material densities
ogo.message("Determining the Material Densities...")
adipose_den = ogo.icMaterialDensity(adipose_HU[0], ic_parameters['Adipose u/p'], ic_parameters['Water u/p'], 1.0)
air_den = ogo.icMaterialDensity(air_HU[0], ic_parameters['Air u/p'], ic_parameters['Water u/p'], 1.0)
blood_den = ogo.icMaterialDensity(blood_HU[0], ic_parameters['Blood u/p'], ic_parameters['Water u/p'], 1.0)
bone_den = ogo.icMaterialDensity(bone_HU[0], ic_parameters['Cortical Bone u/p'], ic_parameters['Water u/p'], 1.0)
muscle_den = ogo.icMaterialDensity(muscle_HU[0], ic_parameters['Skeletal Muscle u/p'], ic_parameters['Water u/p'], 1.0)

material_densities = [adipose_den, air_den, blood_den, bone_den, muscle_den]

##
# Determine the HU-density relationship
ogo.message("Determining the HU-Material Density Relationship...")
hu_rho = ogo.icLinearRegression(mean_hu, material_densities, 'HU-Material Density Slope', 'HU-Material Density Y-Intercept')

##
# Compile the calibration parameters
ogo.message("Compiling the internal calibration parameters...")
cali_parameters = OrderedDict()
cali_parameters['ID'] = org_fileName
cali_parameters['Output File'] = fileName
cali_parameters['Python Script'] = sys.argv[0]
cali_parameters['Version'] = script_version
cali_parameters['Date Created'] = str(date.today())
cali_parameters['Image Directory'] = image_pathname
cali_parameters['Image'] = image_basename
cali_parameters['Mask Directory'] = mask_pathname
cali_parameters['Mask'] = mask_basename
cali_parameters['+++++'] = '+++++'
cali_parameters['Effective Energy [keV]'] = ic_parameters['Effective Energy [keV]']
cali_parameters['Max R^2'] = ic_parameters['Max R^2']
cali_parameters['HU-u/p Slope'] = hu_MUrho['HU-u/p Slope']
cali_parameters['HU-u/p Y-Intercept'] = hu_MUrho['HU-u/p Y-Intercept']
cali_parameters['HU-Material Density Slope'] = hu_rho['HU-Material Density Slope']
cali_parameters['HU-Material Density Y-Intercept'] = hu_rho['HU-Material Density Y-Intercept']
cali_parameters['Adipose u/p'] = ic_parameters['Adipose u/p']
cali_parameters['Air u/p'] = ic_parameters['Air u/p']
cali_parameters['Blood u/p'] = ic_parameters['Blood u/p']
cali_parameters['Cortical Bone u/p'] = ic_parameters['Cortical Bone u/p']
cali_parameters['Skeletal Muscle u/p'] = ic_parameters['Skeletal Muscle u/p']
cali_parameters['K2HPO4 u/p'] = ic_parameters['K2HPO4 u/p']
cali_parameters['CHA u/p'] = ic_parameters['CHA u/p']
cali_parameters['Triglyceride u/p'] = ic_parameters['Triglyceride u/p']
cali_parameters['Water u/p'] = ic_parameters['Water u/p']

# Write the output text file
txt_fileName = org_fileName + "_IntCalibParameters.txt"
ogo.message("Writing parameters to output text file: %s" % txt_fileName)
ogo.writeTXTfile(cali_parameters, txt_fileName, image_pathname)

##
# Apply the internal density calibration to the image
ogo.message("Applying the calibration to the image...")
calibrated_image = ogo.applyInternalCalibration(imageData, cali_parameters)

##
# Write out calibrated image
ogo.message("Writing out the calibrated image: %s" % fileName)
ogo.writeNii(calibrated_image, fileName, image_pathname)



##
# End of script
ogo.message("End of Script.")
ogo.message("Please cite 'Michalski et al. 2020 Med Eng Phys' when using this analysis.")
ogo.message("https://doi.org/10.1016/j.medengphy.2020.01.009")
sys.exit()
