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
import SimpleITK as sitk
from scipy import stats
from datetime import date
from collections import OrderedDict
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

import ogo.dat.MassAttenuationTables as mat
import ogo.dat.OgoMasterLabels as lb
from ogo.calib.internal_calibration import InternalCalibration
import ogo.util.Helper as ogo
from ogo.util.echo_arguments import echo_arguments
from ogo.util.write_txt import write_txt

# PHANTOM CALIBARATION ------------------------------------------------------------------
def phantom(input_image,input_mask,output_image,calib_file_name,async_image,phantom,overwrite,func):
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
        # QCT PRO User Guide, Mindways Software, Inc. v5.0, rev 20110801
        
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
    if calib_file_name:
        ogo.message('Saving calibration parameters to file:')
        ogo.message('      \"{}\"'.format(calib_file_name))
        
        txt_file = open(calib_file_name, "w")
        
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
        txt_file.write('  {:>14s} = {:s}\n'.format('calibration file name',calib_file_name))
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
def internal(input_image,input_mask,output_image,calib_file_name,useLabels,useL4,overwrite,func):
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

    if not (input_image.lower().endswith('.nii') or input_image.lower().endswith('.nii.gz')):
        os.sys.exit('[ERROR] Input must be type NIFTI file: \"{}\"'.format(input_image))

    ogo.message('Reading input CT image to be calibrated:')
    ogo.message('      \"{}\"'.format(input_image))
    ct = sitk.ReadImage(input_image)

    # Read input mask
    if not os.path.isfile(input_mask):
        os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(input_mask))

    if not (input_image.lower().endswith('.nii') or input_image.lower().endswith('.nii.gz')):
        os.sys.exit('[ERROR] Input must be type NIFTI file: \"{}\"'.format(input_mask))

    ogo.message('Reading input mask used for calibration:')
    ogo.message('      \"{}\"'.format(input_mask))
    mask = sitk.ReadImage(input_mask)

    # Dictionary of valid labels for internal calibration.
    labels = OrderedDict()
    for k in (91,92,93,94,95): # if adding a new label, append the ID to the list here
        labels[lb.labels_dict[k].get('LABEL')] = k

    # Search for labels in mask image
    ogo.message('Computing calibration data from image.')
    filt = sitk.LabelStatisticsImageFilter()
    filt.Execute(ct, mask)

    label_list = filt.GetLabels()
    n_labels = len(label_list)
    ogo.message('Found {} labels'.format(n_labels))
    #print('\n'.join('{:8d} ({})'.format(k,lb.labels_dict[k]['LABEL']) for k in np.sort(label_list)))
    for k in np.sort(label_list):
        ogo.message(' {:>22s}'.format(lb.labels_dict[k]['LABEL']+' ('+str(k)+')'))

    # Calculate the values for each of the five possible valid labels (returns 0 if label unavailable)
    labels_data = OrderedDict()
    for label, value in labels.items():

        if not filt.HasLabel(value):
            ogo.message('[WARNING] Could not find values for label \"{}\" ({})'.format(label, value))

        labels_data[label] = {'ID': value, 'mean': filt.GetMean(value), 'stdev': filt.GetVariance(value), 'count': filt.GetCount(value), 'marker': ''}

        if (useL4 and label in 'Cortical Bone'):
            ogo.message('Calculating \"{}\" ({}) from L4'.format(label, value))
            L4_label = 7
            if (not filt.HasLabel(L4_label)):
                os.sys.exit('[ERROR] No L4 found in image.')
            else:
                bone = sitk.MaskImageFilter()
                bone.SetMaskingValue(L4_label)
                array = bone.Execute(ct, mask)
                
                [bone_mean, bone_std, bone_count] = ogo.get_cortical_bone((sitk.GetArrayFromImage(array).ravel()))
                labels_data[label] = {'ID': value, 'mean': bone_mean, 'stdev': bone_std, 'count': bone_count, 'marker': '(from L4)'}
    
    # Print the values for each label
    ogo.message('  {:>22s} {:>8s} {:>8s} {:>8s}'.format(' ','Mean','StdDev','#Voxels'))
    for label, value in labels.items():
        ogo.message('  {:>22s} {:8.3f} {:8.3f} {:8d} {:s}'.format(label+' ('+str(value)+'):',\
            labels_data[label]['mean'],labels_data[label]['stdev'],labels_data[label]['count'],labels_data[label]['marker'])) # Report mean, SD, # voxels
    
    # Finalize the labels to be used (user may select subset)
    labelList = []
    if (useLabels): # User explicitly defines which labels to use
        for labelID in useLabels:
            if (labelID not in labels.values()):
                os.sys.exit('[ERROR] Invalid label selected: ({})'.format(labelID))
            labelList.append(labelID)
    else:
        labelList = [91, 92, 93, 94, 95]
    
    if (len(labelList)<3):
        os.sys.exit('[ERROR] A minimum of three sample tissues needed.')
        
    ogo.message('')
    ogo.message('Labels used for internal calibration:')
    for label,value in labels.items():
        if (value in labelList):
            ogo.message(' {:>22s}'.format(label+' ('+str(value)+')'))
            if (labels_data[label]['mean'] == 0.0):
                os.sys.exit('[ERROR] Invalid HU for {} ({}). Explicitly define labels to use \nor define --useL4.'.format(label,value))
    ogo.message('')

    # Perform the internal calibration fit
    ogo.message('Computing calibration parameters.')
    calib = InternalCalibration(
        adipose_hu=labels_data['Adipose']['mean'],
        air_hu=labels_data['Air']['mean'],
        blood_hu=labels_data['Blood']['mean'],
        bone_hu=labels_data['Cortical Bone']['mean'],
        muscle_hu=labels_data['Skeletal Muscle']['mean'],
        label_list=labelList
    )
    calib.fit()

    ogo.message('  {:>27s} {:8s}'.format('---------------------------','--------'))
    ogo.message('  {:>27s} {:8.3f}'.format('Energy [keV]:',calib.effective_energy))
    ogo.message('  {:>27s} {:8.3f}'.format('Max R^2:',calib.max_r2))
    ogo.message('  {:>27s} {:8s}'.format('---------------------------','--------'))
    ogo.message('  {:>27s}'.format('Mass Attenuation [cm2/g]:'))
    ogo.message('  {:>27s} {:8.3f}'.format('Adipose ',calib.adipose_mass_attenuation))
    ogo.message('  {:>27s} {:8.3f}'.format('Air ',calib.air_mass_attenuation))
    ogo.message('  {:>27s} {:8.3f}'.format('Blood ',calib.blood_mass_attenuation))
    ogo.message('  {:>27s} {:8.3f}'.format('Cortical Bone ',calib.bone_mass_attenuation))
    ogo.message('  {:>27s} {:8.3f}'.format('Skeletal Muscle ',calib.muscle_mass_attenuation))
    ogo.message('  {:>27s} {:8s}'.format('---------------------------','--------'))

    ogo.message('Calibrating input file.')
    voxel_volume = np.prod(ct.GetSpacing())
#    ogo.message('Voxel_volume = {:.3f} mm^3'.format(voxel_volume))
    den = calib.predict(sitk.Cast(ct, sitk.sitkFloat64), voxel_volume)
    den = sitk.Cast(den, ct.GetPixelID())

    ogo.message('  {:>27s} {:8s}'.format('---------------------------','--------'))
    imfilt = sitk.StatisticsImageFilter()
    imfilt.Execute(ct)
    ogo.message('  {:>27s} {:s}'.format('INPUT Image Information:',''))
    ogo.message('  {:>27s} {:8.1f}'.format('Minimum ',imfilt.GetMinimum()))
    ogo.message('  {:>27s} {:8.1f}'.format('Maximum ',imfilt.GetMaximum()))
    ogo.message('  {:>27s} {:8.1f}'.format('Mean ',imfilt.GetMean()))
    ogo.message('  {:>27s} {:8.1f}'.format('Variance ',imfilt.GetVariance()))
    imfilt.Execute(den)
    ogo.message('  {:>27s} {:s}'.format('OUTPUT Image Information:',''))
    ogo.message('  {:>27s} {:8.1f}'.format('Minimum ',imfilt.GetMinimum()))
    ogo.message('  {:>27s} {:8.1f}'.format('Maximum ',imfilt.GetMaximum()))
    ogo.message('  {:>27s} {:8.1f}'.format('Mean ',imfilt.GetMean()))
    ogo.message('  {:>27s} {:8.1f}'.format('Variance ',imfilt.GetVariance()))
    ogo.message('  {:>27s} {:8s}'.format('---------------------------','--------'))
    
    ogo.message('Writing result to ' + output_image)
    sitk.WriteImage(den, output_image)

    #print('!> {:30s} = {}'.format('TimeDimension',reader.GetTimeDimension()))
    #print('!> {:30s} = {}'.format('TimeSpacing',reader.GetTimeSpacing()))
    #print('!> {:30s} = {}'.format('RescaleSlope',reader.GetRescaleSlope()))
    #print('!> {:30s} = {}'.format('RescaleIntercept',reader.GetRescaleIntercept()))
    #print('!> {:30s} = {}'.format('QFac',reader.GetQFac()))
    #print('!> {:30s} = {}'.format('QFormMatrix',reader.GetQFormMatrix()))
    #print('!> {:30s} = {}'.format('NIFTIHeader',reader.GetNIFTIHeader()))
    
    if calib_file_name:
        ogo.message('Saving calibration parameters to file:')
        ogo.message('      \"{}\"'.format(calib_file_name))
        
        txt_file = open(calib_file_name, "w")
        
        txt_file.write('Internal calibration:\n')
        txt_file.write('  {:>27s} {:8s}\n'.format('---------------------------','--------'))
        txt_file.write('  {:>27s} {:s}\n'.format('ID:',os.path.basename(output_image)))
        txt_file.write('  {:>27s} {:s}\n'.format('python script:',os.path.splitext(os.path.basename(sys.argv[0]))[0]))
        txt_file.write('  {:>27s} {:.2f}\n'.format('version:',script_version))
        txt_file.write('  {:>27s} {:s}\n'.format('creation date:',str(date.today())))
        txt_file.write('\n')
        txt_file.write('Files:\n')
        txt_file.write('  {:>27s} {:8s}\n'.format('---------------------------','--------'))
        txt_file.write('  {:>27s} {:s}\n'.format('input image:',input_image))
        txt_file.write('  {:>27s} {:s}\n'.format('input mask:',input_mask))
        txt_file.write('  {:>27s} {:s}\n'.format('output image:',output_image))
        txt_file.write('  {:>27s} {:s}\n'.format('calibration file name:',calib_file_name))
        txt_file.write('\n')
        txt_file.write('Calibration parameters:\n')
        txt_file.write('  {:>27s} {:8s}\n'.format('---------------------------','--------'))
        txt_file.write('  {:>27s} {}\n'.format('Fit:',calib._is_fit))
        txt_file.write('  {:>27s} {:8.3f}\n'.format('Energy [keV]:',calib.effective_energy))
        txt_file.write('  {:>27s} {:8.3f}\n'.format('Max R^2:',calib.max_r2))
        txt_file.write('  {:>27s} {:8s}\n'.format('---------------------------','--------'))
        txt_file.write('  {:>27s}\n'.format('Density [HU]:'))
        txt_file.write('  {:>27s} {:8.3f}\n'.format('Adipose ',calib.adipose_hu))
        txt_file.write('  {:>27s} {:8.3f}\n'.format('Air ',calib.air_hu))
        txt_file.write('  {:>27s} {:8.3f}\n'.format('Blood ',calib.blood_hu))
        txt_file.write('  {:>27s} {:8.3f}\n'.format('Cortical Bone ',calib.bone_hu))
        txt_file.write('  {:>27s} {:8.3f}\n'.format('Skeletal Muscle ',calib.muscle_hu))
        txt_file.write('  {:>27s} {:8s}\n'.format('---------------------------','--------'))
        txt_file.write('  {:>27s}\n'.format('Mass Attenuation [cm2/g]:'))
        txt_file.write('  {:>27s} {:8.3f}\n'.format('Adipose ',calib.adipose_mass_attenuation))
        txt_file.write('  {:>27s} {:8.3f}\n'.format('Air ',calib.air_mass_attenuation))
        txt_file.write('  {:>27s} {:8.3f}\n'.format('Blood ',calib.blood_mass_attenuation))
        txt_file.write('  {:>27s} {:8.3f}\n'.format('Cortical Bone ',calib.bone_mass_attenuation))
        txt_file.write('  {:>27s} {:8.3f}\n'.format('Skeletal Muscle ',calib.muscle_mass_attenuation))
        txt_file.write('\n')
        txt_file.write('  {:>27s} {:8.3f}\n'.format('K2HPO4 ',calib.K2HPO4_mass_attenuation))
        txt_file.write('  {:>27s} {:8.3f}\n'.format('CHA ',calib.CHA_mass_attenuation))
        txt_file.write('  {:>27s} {:8.3f}\n'.format('Triglyceride ',calib.triglyceride_mass_attenuation))
        txt_file.write('  {:>27s} {:8.3f}\n'.format('Water ',calib.water_mass_attenuation))
        txt_file.write('\n')
        txt_file.write('  {:>27s} {:8.3f}\n'.format('Water ',calib.water_mass_attenuation))
        txt_file.write('  {:>27s} {:8.3f}\n'.format('HU-u/p Slope',calib.hu_to_mass_attenuation_slope))
        txt_file.write('  {:>27s} {:8.3f}\n'.format('HU-u/p Y-Intercept',calib.hu_to_mass_attenuation_intercept))
        txt_file.write('  {:>27s} {:8.3f}\n'.format('HU-Material Density Slope',calib.hu_to_density_slope))
        txt_file.write('  {:>27s} {:8.3f}\n'.format('HU-Material Density Y-Intercept',calib.hu_to_density_intercept))
        
        txt_file.write('  {:>27s} {:8s}\n'.format('---------------------------','--------'))
        txt_file.write('\n')
        txt_file.write('Unformatted calibration parameters:\n')
        txt_file.write('  {:>27s} {:8s}\n'.format('---------------------------','--------'))
        for label,value in calib.get_dict().items():
            txt_file.write('  {:>27s} {}\n'.format(label,value))
        txt_file.write('  {:>27s} {:8s}\n'.format('---------------------------','--------'))

        txt_file.close()
        
    ogo.message('Done internal calibration.')
    
def main():
    # Setup description
    description='''
Performs quantitative density calibration of a clinical CT image either
by phantom calibration or internal calibration. 

Phantom calibration can be performed with the phantom in the image or
by an asynchronous scan of the phantom. The image mask is used to define
the rods in the image. If the rods are in an asynchronous scan, then set
the appropriate argument.

Currently implemented phantoms are:
  Mindways Model 3 CT
  Mindways Model 3 QA
  QRM-BDC 3-rod
  QRM-BDC 6-rod
  Image Analysis QCT-3D Plus
  B-MAS 200

Internal calibration is typically based on five samples in the image, 
which are adipose tissue, skeletal muscle, cortical bone, blood, and air.
However, the user is free to define any number of samples, and samples of
any type. The requirement is to use the correct image labels. At least
two labels must be defined, and current valid labels are:

   91 Adipose tissue
   92 Air
   93 Blood (artery)
   94 Cortical Bone
   95 Skeletal Muscle
 
Outputs include the calibrated density image and a text file that includes
the calibration parameters and results.

Citation:
Michalski AS, Besler BA, Michalak GJ, Boyd SK, 2020. CT-based internal density 
calibration for opportunistic skeletal assessment using abdominal CT scans. 
Med Eng Phys 78, 55-63.

'''

    epilog='''
Example calls: 
ogoImageCalibration phantom image.nii.gz rod_mask.nii.gz \\
                            image_qct.nii.gz  
ogoImageCalibration phantom image.nii.gz rod_mask.nii.gz \\
                            image_qct.nii.gz --phantom 'QRM-BDC 3-rod' 
ogoImageCalibration phantom image.nii.gz rod_mask.nii.gz \\
                            image_qct.nii.gz --async_image asynch_image.nii.gz 
ogoImageCalibration internal image.nii.gz samples_mask.nii.gz \\
                            image_qct.nii.gz --useL4
ogoImageCalibration internal image.nii.gz samples_mask.nii.gz \\
                            image_qct.nii.gz --useLabels 91 92 93 95

ogoImageCalibration phantom \
  /Users/skboyd/Desktop/ML/test/kub.nii.gz \
  /Users/skboyd/Desktop/ML/test/kub_mask.nii.gz \
  /Users/skboyd/Desktop/ML/test/test.nii \
  --phantom 'Mindways Model 3 CT' \
  --calib_file_name /Users/skboyd/Desktop/ML/test/test.txt \
  --overwrite
    
ogoImageCalibration phantom \
  /Users/skboyd/Desktop/ML/test/kub.nii.gz \
  /Users/skboyd/Desktop/ML/test/async_mask_mindways.nii.gz \
  /Users/skboyd/Desktop/ML/test/test.nii \
  --async_image /Users/skboyd/Desktop/ML/test/async.nii.gz \
  --calib_file_name /Users/skboyd/Desktop/ML/test/test.txt \
  --phantom 'Mindways Model 3 CT' 

ogoImageCalibration phantom \
  /Users/skboyd/Desktop/ML/test/kub.nii.gz \
  /Users/skboyd/Desktop/ML/test/async_mask_bmas200.nii.gz \
  /Users/skboyd/Desktop/ML/test/test.nii \
  --async_image /Users/skboyd/Desktop/ML/test/async.nii.gz \
  --calib_file_name /Users/skboyd/Desktop/ML/test/test.txt \
  --phantom 'B-MAS 200' 

ogoImageCalibration internal \
  /Users/skboyd/Desktop/ML/test/retro.nii \
  /Users/skboyd/Desktop/ML/test/retro_mask.nii.gz \
  /Users/skboyd/Desktop/ML/test/test.nii \
  --calib_file_name /Users/skboyd/Desktop/ML/test/test.txt \
  --overwrite --useL4

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
    parser_phantom.add_argument('--calib_file_name', help='Calibration results file (*.txt)')
    parser_phantom.add_argument('--async_image', default='', metavar='IMAGE', help='Asynchronous image of phantom (*.nii, *.nii.gz)')
    parser_phantom.add_argument('--phantom', default='Mindways Model 3 CT', choices=['Mindways Model 3 CT','Mindways Model 3 QA','QRM-BDC 3-rod','QRM-BDC 6-rod','Image Analysis QCT-3D Plus','B-MAS 200'], help='Specify phantom used (default: %(default)s)')
    parser_phantom.add_argument('--overwrite', action='store_true', help='Overwrite output without asking')
    parser_phantom.set_defaults(func=phantom)

    # Internal
    parser_internal = subparsers.add_parser('internal')
    parser_internal.add_argument('input_image', help='Input image file (*.nii, *.nii.gz)')
    parser_internal.add_argument('input_mask', help='Input image mask file (*.nii, *.nii.gz)')
    parser_internal.add_argument('output_image', help='Output image file (*.nii, *.nii.gz)')
    parser_internal.add_argument('--calib_file_name', help='Calibration results file (*.txt)')
    parser_internal.add_argument('--useLabels', type=int, nargs='*', default=[], metavar='ID', help='Explicitly define labels for internal calibration; space separated (e.g. 91 92 93 94 95) (default: all)')
    parser_internal.add_argument('--useL4', action='store_true', help='Use when label 94 (cortical bone) is not available (it samples L4).')
    parser_internal.add_argument('--overwrite', action='store_true', help='Overwrite output without asking')
    parser_internal.set_defaults(func=internal)
    
    # Parse and display
    args = parser.parse_args()
    print(echo_arguments('ImageCalibration', vars(args)))
    
    # Run program
    args.func(**vars(args))

if __name__ == '__main__':
    main()