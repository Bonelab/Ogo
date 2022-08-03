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
import ogo.dat.MassAttenuationTables as mat
import ogo.dat.OgoMasterLabels as lb
from scipy import stats
from datetime import date
from collections import OrderedDict
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

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
    ogo.message('Calibration phantom selected is:')
    ogo.message('  {:>14s} = {:s}'.format('name',phantom_dict['name']))
    ogo.message('  {:>14s} = {:s}'.format('type',phantom_dict['type']))
    ogo.message('  {:>14s} = {:s}'.format('serial',phantom_dict['serial']))
    ogo.message('  {:>14s} = {:d}'.format('number_rods',phantom_dict['number_rods']))
    ogo.message("  {:>14s} = ".format('rod_labels')+'['+', '.join("{:d}".format(i) for i in phantom_dict['rod_labels'])+"]")
    ogo.message("  {:>14s} = ".format('rod_names')+'['+', '.join("{:s}".format(i) for i in phantom_dict['rod_names'])+"]")
    ogo.message("  {:>14s} = ".format('densities')+'['+', '.join("{:8.3f}".format(i) for i in phantom_dict['densities'])+"]")
    if phantom_dict['h2o_densities'][0] != None:
        ogo.message("  {:>14s} = ".format('h2o_densities')+'['+', '.join("{:8.3f}".format(i) for i in phantom_dict['h2o_densities'])+"]")

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
        
    ogo.message("  {:>14s} = ".format('Image [HU]')+'['+', '.join("{:8.3f}".format(i) for i in phantom_HU)+"]")

    # Calculate slope and intercept
    ogo.message('Calculating calibration slope and intercept.')
    ogo.message('Calculating type is: \'{}\''.format(phantom_dict['type']))
    if phantom_dict['type'] is 'k2hpo4':
        x_values = phantom_dict['densities']
        y_values = np.subtract(phantom_HU,phantom_dict['h2o_densities'])

        regression_parameters = stats.linregress(x_values, y_values)
        sigma_ct = regression_parameters[0] - 0.2174
        beta_ct = regression_parameters[1] + 999.6
        calibration_slope = 1/sigma_ct
        calibration_yint = 1*beta_ct/sigma_ct
        calibration_rvalue = regression_parameters[2]
        calibration_stderr = regression_parameters[3]
        ogo.message('  {:>14s} = {:8.6f}'.format('slope',calibration_slope))
        ogo.message('  {:>14s} = {:8.6f}'.format('y-intercept',calibration_yint))
        ogo.message('  {:>14s} = {:8.6f}'.format('r-value',calibration_rvalue))
        ogo.message('  {:>14s} = {:8.6f}'.format('stderr',calibration_stderr))
        
    elif phantom_dict['type'] is 'CHA':
        x_values = phantom_HU
        y_values = phantom_dict['densities']

        regression_parameters = stats.linregress(x_values, y_values)
        calibration_slope = regression_parameters[0]
        calibration_yint = regression_parameters[1]
        calibration_rvalue = regression_parameters[2]
        calibration_stderr = regression_parameters[3]
        ogo.message('  {:>14s} = {:8.6f}'.format('slope',calibration_slope))
        ogo.message('  {:>14s} = {:8.6f}'.format('y-intercept',calibration_yint))
        ogo.message('  {:>14s} = {:8.6f}'.format('r-value',calibration_rvalue))
        ogo.message('  {:>14s} = {:8.6f}'.format('stderr',calibration_stderr))

    else:
        os.sys.exit('[ERROR] Unknown calibration phantom type: \"{}\"'.format(phantom_dict['type']))
    
    # Check that the calibration results are within a reasonable range
    threshold_warning = 2 # percent
    threshold_error =  10 # percent
    if (abs(1.0-calibration_rvalue)*100.0>threshold_error):
        ogo.message('[ERROR] R-value is too low ({:8.6f}).'.format(calibration_rvalue))
        ogo.message('          * did you select the correct phantom?')
        ogo.message('          * do all the rods have a reasonable value?')
        ogo.message('          * is the StdDev of the rods high (i.e. missed the target)?')
        os.sys.exit()
    if (abs(1.0-calibration_rvalue)*100.0>threshold_warning):
        ogo.message('[WARNING] R-value is low ({:8.6f}). Checking following:'.format(calibration_rvalue))
        ogo.message('          * did you select the correct phantom?')
        ogo.message('          * do all the rods have a reasonable value?')
        ogo.message('          * is the StdDev of the rods high (i.e. missed the target)?')
        
    # Apply calibration parameters to image
    ogo.message("Applying image calibration.")
    
    np_image = vtk_to_numpy(reader_image.GetOutput().GetPointData().GetScalars())
    np_image = np_image * calibration_slope + calibration_yint
    vtk_image_scalars = numpy_to_vtk(num_array=np_image, deep=True, array_type=vtk.VTK_SHORT) # Cast to short

    calibrated_image = vtk.vtkImageData()
    calibrated_image.DeepCopy(reader_image.GetOutput())
    calibrated_image.GetPointData().SetScalars(vtk_image_scalars)
    
    # Write image
    if output_image.lower().endswith('.nii'):
        writer = vtk.vtkNIFTIImageWriter()
    elif output_image.lower().endswith('.nii.gz'):
        writer = vtk.vtkNIFTIImageWriter()
    else:
        os.sys.exit('[ERROR] Cannot find writer for file \"{}\"'.format(output_image))
          
    ogo.message('Saving output image ' + output_image)

    writer.SetInputData(calibrated_image)
    writer.SetFileName(output_image)
    writer.SetTimeDimension(reader_image.GetTimeDimension())
    writer.SetTimeSpacing(reader_image.GetTimeSpacing())
    writer.SetRescaleSlope(reader_image.GetRescaleSlope())
    writer.SetRescaleIntercept(reader_image.GetRescaleIntercept())
    writer.SetQFac(reader_image.GetQFac())
    writer.SetQFormMatrix(reader_image.GetQFormMatrix())
    writer.SetNIFTIHeader(reader_image.GetNIFTIHeader())
    writer.Update()
    
    # Write text file
    output_text = os.path.splitext(output_image)[0] + '.txt'
    ogo.message('Writing parameters to output text file...')
    ogo.message('      \"{}\"'.format(output_text))
    
    txt_file = open(output_text, "w")
    
    txt_file.write('Phantom-based calibration:\n')
    txt_file.write('  {:>14s} = {:s}\n'.format('ID',os.path.basename(output_image)))
    txt_file.write('  {:>14s} = {:s}\n'.format('python script',os.path.splitext(os.path.basename(sys.argv[0]))[0]))
    txt_file.write('  {:>14s} = {:.2f}\n'.format('version',script_version))
    txt_file.write('  {:>14s} = {:s}\n'.format('creation date',str(date.today())))
    txt_file.write('\n')
    txt_file.write('Files:\n')
    txt_file.write('  {:>14s} = {:s}\n'.format('input image',input_image))
    txt_file.write('  {:>14s} = {:s}\n'.format('input mask',input_mask))
    txt_file.write('  {:>14s} = {:s}\n'.format('output image',output_image))
    if async_image:
        txt_file.write('  {:>14s} = {:s}\n'.format('async image',async_image))
    else:
        txt_file.write('  {:>14s} = {:s}\n'.format('async image','n/a'))
    txt_file.write('\n')
    txt_file.write('Calibration phantom:\n')
    txt_file.write('  {:>14s} = {:s}\n'.format('name',phantom_dict['name']))
    txt_file.write('  {:>14s} = {:s}\n'.format('type',phantom_dict['type']))
    txt_file.write('  {:>14s} = {:s}\n'.format('serial',phantom_dict['serial']))
    txt_file.write('  {:>14s} = {:d}\n'.format('number_rods',phantom_dict['number_rods']))
    txt_file.write('  {:>14s} = '.format('rod_labels')+'['+', '.join("{:d}".format(i) for i in phantom_dict['rod_labels'])+"]\n")
    txt_file.write('  {:>14s} = '.format('rod_names')+'['+', '.join("{:s}".format(i) for i in phantom_dict['rod_names'])+"]\n")
    txt_file.write('  {:>14s} = '.format('densities')+'['+', '.join("{:8.3f}".format(i) for i in phantom_dict['densities'])+"]\n")
    if phantom_dict['h2o_densities'][0] != None:
        txt_file.write("  {:>14s} = ".format('h2o_densities')+'['+', '.join("{:8.3f}".format(i) for i in phantom_dict['h2o_densities'])+"]\n")
    txt_file.write('\n')
    txt_file.write('Calibration results:\n')
    txt_file.write('  {:>14s} = {:8.6f}\n'.format('slope',calibration_slope))
    txt_file.write('  {:>14s} = {:8.6f}\n'.format('y-intercept',calibration_yint))
    txt_file.write('  {:>14s} = {:8.6f}\n'.format('r-value',calibration_rvalue))
    txt_file.write('  {:>14s} = {:8.6f}\n'.format('stderr',calibration_stderr))
    
    txt_file.close()
    
    ogo.message('Done ImageCalibration.')
    
    
# INTERNAL CALIBARATION ------------------------------------------------------------------
def internal(input_image,input_mask,output_image,excludeLabels,overwrite,func):
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

    # Create list of valid internal calibration labels
    # Define the dictionary of possible valid labels. If in the future new labels
    # are defined, add them to this dictionary.
    valid_labels_dict = OrderedDict()
    valid_labels_dict[91] = 'Adipose'
    valid_labels_dict[92] = 'Air'
    valid_labels_dict[93] = 'Blood'
    valid_labels_dict[94] = 'Cortical Bone'
    valid_labels_dict[95] = 'Skeletal Muscle'
    all_valid_labels = valid_labels_dict.keys() # We need a 'list' of valid labels

    array = vtk_to_numpy(reader_mask.GetOutput().GetPointData().GetScalars()).ravel()
    labels_found = np.unique(array)
    valid_labels = []
    for idx,label in enumerate(all_valid_labels): # Final list only includes labels present in mask image
        if label in labels_found:
            valid_labels.append(label)
    
    ogo.message('All sample labels available:')
    for idx,label in enumerate(valid_labels):
        ogo.message('  {:>16s} = {:d}'.format(valid_labels_dict.get(label),label))
        
    if excludeLabels:
        ogo.message('Excluding:')
        for idx,label in enumerate(excludeLabels):
            if label in valid_labels:
                ogo.message('  {:>16s} = {:d}'.format(valid_labels_dict.get(label),label))
                valid_labels.pop(valid_labels.index(label))
    else:
        ogo.message('No labels to be excluded.')
        
    ogo.message('Final labels:')
    for idx,label in enumerate(valid_labels):
        ogo.message('  {:>16s} = {:d}'.format(valid_labels_dict.get(label),label))

    if len(valid_labels)<3:
        os.sys.exit('[ERROR] A minimum of three samples are needed for internal calibration.')
    
    
    # Gather the image and mask
    image_array = vtk_to_numpy(reader_image.GetOutput().GetPointData().GetScalars())
    mask_array = vtk_to_numpy(reader_mask.GetOutput().GetPointData().GetScalars())
    
    # Calculate the mean and SD of each phantom rod in the image
    ogo.message('Extracting sample calibration data from image.')
    ogo.message('  {:>22s} {:>8s} {:>8s} {:>8s}'.format(' ','Mean','StdDev','#Voxels'))
    sample_HU = []
    
    voxel_warning = False
    voxel_warning_thres = 100 # Minimum number of voxels
    std_warning = False
    std_warning_thres = 100 # Maximum stdev
    
    for idx,label in enumerate(valid_labels):
        #rod_name = rod_names[rod]
        sample_mean = np.mean(image_array[mask_array == label])
        sample_std = np.std(image_array[mask_array == label])
        sample_count = len(image_array[mask_array == label])
        if (sample_count < voxel_warning_thres):
            voxel_warning = True
        if (sample_std > std_warning_thres):
            std_warning = True
        ogo.message('  {:>22s} {:8.3f} {:8.3f} {:8d}'.format(valid_labels_dict.get(label)+' ('+str(label)+'):',sample_mean,sample_std,sample_count)) # Report mean, SD, # voxels
        sample_HU.append(sample_mean)

    ogo.message("  {:>16s} = ".format('Image [HU]')+'['+', '.join("{:8.3f}".format(i) for i in sample_HU)+"]")    
    if (voxel_warning):
        ogo.message('[WARNING] At least one sample has less than {} voxels. Caution!'.format(voxel_warning_thres))
    if (std_warning):
        ogo.message('[WARNING] At least one sample SD greater than {:.2f}. Caution!'.format(std_warning_thres))
            
    #print(mat.adipose_table)
    print(lb.master_labels_dict)
    
    ogo.message('Done internal calibration.')
    
def main():
    # Setup description
    description='''
Performs quantitative density calibration of a clinical CT image either
by phantom calibration or internal calibration. 

Phantom calibration can be performed with the phantom in the image or
in an asynchronous scan of the phantom. The image mask is used to define
the rods in the image. If the rods are in an asynchronous scan, then set
the argument --async_image and the mask will be applied to that image.

Current phantoms are:
  Mindways Model 3 CT
  Mindways Model 3 QA
  QRM-BDC 3-rod
  QRM-BDC 6-rod
  Image Analysis QCT-3D Plus
  B-MAS 200

Internal calibration is typically based on five samples in the image, 
which are adipose tissue, skeletal muscle, cortical bone, blood, and air.
However, the user is free to define any number of samples, and samples of
any type. The only requirement is to use the correct image labels. At least
three labels must be defined, and current valid labels are:

   91 Adipose tissue
   92 Air
   93 Blood (artery)
   94 Cortical Bone
   95 Skeletal Muscle
 
Valid file formats are NIFTII: .nii, .nii.gz

Outputs include the calibrated density image and a text file that includes
the calibration parameters and results.

Citation:
Michalski AS, Besler BA, Michalak GJ, Boyd SK, 2020. CT-based internal density 
calibration for opportunistic skeletal assessment using abdominal CT scans. 
Med Eng Phys 78, 55-63.

'''

    epilog='''
USAGE: 
ogoImageCalibration phantom image.nii.gz rod_mask.nii.gz \\
                            image_qct.nii.gz  
ogoImageCalibration phantom image.nii.gz rod_mask.nii.gz \\
                            image_qct.nii.gz --phantom 'QRM-BDC 3-rod' 
ogoImageCalibration phantom image.nii.gz rod_mask.nii.gz \\
                            image_qct.nii.gz --async_image asynch_image.nii.gz 
ogoImageCalibration internal image.nii.gz samples_mask.nii.gz \\
                            image_qct.nii.gz  

ogoImageCalibration phantom \
  /Users/skboyd/Desktop/ML/test/kub.nii.gz \
  /Users/skboyd/Desktop/ML/test/kub_mask.nii.gz \
  /Users/skboyd/Desktop/ML/test/test.nii \
  --phantom 'Mindways Model 3 CT' \
  --overwrite
    
ogoImageCalibration phantom \
  /Users/skboyd/Desktop/ML/test/kub.nii.gz \
  /Users/skboyd/Desktop/ML/test/async_mask_mindways.nii.gz \
  /Users/skboyd/Desktop/ML/test/test.nii \
  --async_image /Users/skboyd/Desktop/ML/test/async.nii.gz \
  --phantom 'Mindways Model 3 CT' 

ogoImageCalibration phantom \
  /Users/skboyd/Desktop/ML/test/kub.nii.gz \
  /Users/skboyd/Desktop/ML/test/async_mask_bmas200.nii.gz \
  /Users/skboyd/Desktop/ML/test/test.nii \
  --async_image /Users/skboyd/Desktop/ML/test/async.nii.gz \
  --phantom 'B-MAS 200' 

ogoImageCalibration internal \
  /Users/skboyd/Desktop/ML/test/retro.nii \
  /Users/skboyd/Desktop/ML/test/retro_mask.nii.gz \
  /Users/skboyd/Desktop/ML/test/test.nii \
  --overwrite

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
    parser_phantom.add_argument('input_mask', help='Image mask of rods for input image or asynchronous image (*.nii, *.nii.gz)')
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
    parser_internal.add_argument('--excludeLabels', type=int, nargs='*', default=[], metavar='ID', help='Labels to be excluded from internal calibration; space separated (e.g. 93 94)')
    parser_internal.add_argument('--overwrite', action='store_true', help='Overwrite output without asking')
    parser_internal.set_defaults(func=internal)
    
    # Parse and display
    args = parser.parse_args()
    print(echo_arguments('ImageCalibration', vars(args)))
    
    # Run program
    args.func(**vars(args))

if __name__ == '__main__':
    main()