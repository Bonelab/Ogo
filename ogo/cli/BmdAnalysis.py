#####
# ogo_bmd_analysis.py
#
# This script performs BMD analysis from an input calibrated image and a bone mask. This
# script outputs a text file containing the BMD analysis results.
#
#####
#
# Andrew Michalski
# University of Calgary
# Biomedical Engineering Graduate Program
# April 24, 2019
# Modified Aug 2, 2019
# Modified to Py3: March 25, 2020
#####

script_version = 1.0

##
# Import the required modules
import ogo.cli.Helper as ogo
import os
import sys
import argparse
import time
from datetime import date
from collections import OrderedDict


from ogo.util.echo_arguments import echo_arguments

##
# Start Script
def bmdAnalysis(args):
    ogo.message("Start of Script...")

    ##
    # Collect the input arguments
    image = args.calibrated_image
    mask = args.bone_mask
    mask_threshold = args.mask_threshold

    ##
    # Determine image locations and names of files
    image_pathname = os.path.dirname(image)
    image_basename = os.path.basename(image)
    mask_pathname = os.path.dirname(mask)
    mask_basename = os.path.basename(mask)

    ##
    # Read input image
    ogo.message("Reading calibrated image...")
    imageData = ogo.readNii(image)

    ##
    # Read bone mask
    ogo.message("Reading bone mask...")
    maskData = ogo.readNii(mask)

    ##
    # Apply the mask to the data
    ogo.message("Extracting bone VOI from the image...")
    bone_thres = ogo.maskThreshold(maskData, mask_threshold)

    if mask_threshold >= 6 and mask_threshold <= 10:
        image_resample = ogo.imageResample(imageData, 1.0)
        mask_resample = ogo.imageResample(bone_thres, 1.0)
        ogo.message("Extracting vertebral body from full vertebra...")
        vertebral_body_image, vertebral_body_mask = ogo.vertebralBodyExtract(image_resample, mask_resample)
        imageData = vertebral_body_image
        # bone_thres = vertebral_body_mask

    bone_VOI = ogo.applyMask(imageData, bone_thres)

    ##
    # Compute the BMD metrics
    ogo.message("Computing the BMD metrics for the bone...")
    bmd_outcomes = ogo.bmd_metrics(bone_VOI)

    ##
    # Write output Text file
    fileName = mask_basename.replace("_PERI_CORR.nii", "")
    if mask_threshold == 1:
        fileName = fileName + "_RT_FEMUR"
    if mask_threshold == 2:
        fileName = fileName + "_LT_FEMUR"
    if mask_threshold == 3:
        fileName = fileName + "_RT_PELVIS"
    if mask_threshold == 4:
        fileName = fileName + "_LT_PELVIS"
    if mask_threshold == 5:
        fileName = fileName + "_SACRUM"
    if mask_threshold == 6:
        fileName = fileName + "_L5"
    if mask_threshold == 7:
        fileName = fileName + "_L4"
    if mask_threshold == 8:
        fileName = fileName + "_L3"
    if mask_threshold == 9:
        fileName = fileName + "_L2"
    if mask_threshold == 10:
        fileName = fileName + "_L1"

    txt_fileName = fileName + "_BMD_Results.txt"
    parameters_dict = OrderedDict()
    parameters_dict['ID'] = fileName
    parameters_dict['Python Script'] = sys.argv[0]
    parameters_dict['Version'] = script_version
    parameters_dict['Date Created'] = str(date.today())
    parameters_dict['Calibrated Image Directory'] = image_pathname
    parameters_dict['Calibrated Image'] = image_basename
    parameters_dict['Bone Mask Directory'] = mask_pathname
    parameters_dict['Bone Mask'] = mask_basename
    parameters_dict['Bone Mask Threshold Value'] = mask_threshold
    parameters_dict['+++++++++++'] = '++++++++++++'
    parameters_dict['Integral BMD [mg/cc]'] = bmd_outcomes['Integral BMD [mg/cc]']
    parameters_dict['Integral BMC [mg]'] = bmd_outcomes['Integral BMC [mg]']
    parameters_dict['Bone Volume [mm^3]'] = bmd_outcomes['Bone Volume [mm^3]']
    parameters_dict['Bone Volume [cm^3]'] = bmd_outcomes['Bone Volume [cm^3]']

    ogo.message("Writing parameters to output text file: %s" % txt_fileName)
    ogo.writeTXTfile(parameters_dict, txt_fileName, image_pathname)

    ##
    # End of script
    ogo.message("End of Script.")
    sys.exit()

def main():
    description = '''
    This script performs BMD analysis from an input calibrated image and a bone mask. 
    
    This script outputs a text file containing the BMD analysis results. 
    
    INPUT: Calibrated Image (as NIFTI), Bone Mask (as NIFTI) 
    OUTPUT: Text file of BMD analysis and associated information.
    
    '''


    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoBMDAnalysis",
        description=description
    )
    

    parser.add_argument("calibrated_image",
        help = "*_K2HPO4.nii image file")
    parser.add_argument("bone_mask",
        help = "*_MASK.nii mask image of bone")

    parser.add_argument("--mask_threshold", type = int,
        default = 1,
        help = "Set the threshold value to extract the bone of interest from the mask. (Default: %(default)s)")

    # Parse and display
    args = parser.parse_args()
    
    print(echo_arguments('InternalCalibration', vars(args)))

    # Run program
    bmdAnalysis(args)

if __name__ == '__main__':
    main()