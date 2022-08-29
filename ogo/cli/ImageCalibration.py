# /------------------------------------------------------------------------------+
# | 19-JUL-2022                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+
#
# This script is based on work by Andrew Michalski as part of his PhD. In 2022 it was
# updated to incorporate a templated version of both internal and phantom calibration.
#
# Steve Boyd, 2022
#
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
from ogo.calib.mindways_calibration import MindwaysCalibration
from ogo.calib.standard_calibration import StandardCalibration
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

    if not (input_image.lower().endswith('.nii') or input_image.lower().endswith('.nii.gz')):
        os.sys.exit('[ERROR] Input must be type NIFTI file: \"{}\"'.format(input_image))

    ogo.message('Reading input CT image to be calibrated:')
    ogo.message('      \"{}\"'.format(input_image))
    ct = sitk.ReadImage(input_image)

    # Read input mask
    if not os.path.isfile(input_mask):
        os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(input_mask))

    if not (input_mask.lower().endswith('.nii') or input_mask.lower().endswith('.nii.gz')):
        os.sys.exit('[ERROR] Input must be type NIFTI file: \"{}\"'.format(input_mask))

    ogo.message('Reading input mask used for calibration:')
    ogo.message('      \"{}\"'.format(input_mask))
    rods = sitk.ReadImage(input_mask)

    # Read asynchronous input image, if available
    if (async_image):
        if not os.path.isfile(async_image):
            os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(async_image))
    
        if not (async_image.lower().endswith('.nii') or async_image.lower().endswith('.nii.gz')):
            os.sys.exit('[ERROR] Input must be type NIFTI file: \"{}\"'.format(async_image))
    
        ogo.message('Reading asynchronous image of phantom...')
        ogo.message('      \"{}\"'.format(async_image))
        async_ct = sitk.ReadImage(async_image)
    
        ogo.message('NOTE: Asynchronous calibration selected.')
        ogo.message('      Using mask file against asynchronous file:')
        ogo.message('      \"{}\"'.format(input_mask))
        ogo.message('      \"{}\"'.format(async_image))
    
    # Get the correct reference calibration phantom
    phantom_dict = ogo.get_phantom(phantom)
    ogo.message('Calibration phantom:')
    ogo.message('  {:>14s} = {:s}'.format('name',phantom_dict['name']))
    ogo.message('  {:>14s} = {:s}'.format('type',phantom_dict['type']))
    ogo.message('  {:>14s} = {:s}'.format('serial',phantom_dict['serial']))
    ogo.message('  {:>14s} = {:d}'.format('number_rods',phantom_dict['number_rods']))
    ogo.message("  {:>14s} = ".format('rod_labels')+'['+', '.join("{:d}".format(i) for i in phantom_dict['rod_labels'])+"]")
    ogo.message("  {:>14s} = ".format('rod_names')+'['+', '.join("{:s}".format(i) for i in phantom_dict['rod_names'])+"]")
    ogo.message("  {:>14s} = ".format('densities')+'['+', '.join("{:8.3f}".format(i) for i in phantom_dict['densities'])+"]")
    if phantom_dict['h2o_densities'][0] != None:
        ogo.message("  {:>14s} = ".format('h2o_densities')+'['+', '.join("{:8.3f}".format(i) for i in phantom_dict['h2o_densities'])+"]")

    filt = sitk.LabelStatisticsImageFilter()
    if (async_image):
        filt.Execute(async_ct, rods)
    else:
        filt.Execute(ct, rods)
    
    all_labels = np.sort(filt.GetLabels())
    #all_labels = all_labels[all_labels!=0] # remove zero label from list
    n_labels = len(all_labels)
    label_list = []
    
    ogo.message('Found {} labels in\n'
                '               {}'.format(n_labels,input_mask))
                
    for label in all_labels:
        ogo.message(' {:>14s} {}'.format('('+str(label)+')',lb.labels_dict[label]['LABEL']))
        if label in phantom_dict['rod_labels']: # take only the labels for rods and ignore others
            label_list.append(label)

    if len(label_list)<2:
        os.sys.exit('[ERROR] Need at least 2 rods to fit a function.\n'
                    '        Found {} rods.'.format(len(label_list)))
    
    ogo.message('Determining mean intensities.')
    ogo.message('Image HU for rods:')
    ogo.message("       {:>8s} {:>8s} {:>8s}".format('Mean','StdDev','#Voxels'))

    HU = []
    for idx,label in enumerate(label_list):
        rod_name = phantom_dict['rod_names'][idx]
        rod_mean = filt.GetMean(int(label))
        rod_std = filt.GetVariance(int(label))
        rod_count = filt.GetCount(int(label))
        ogo.message('Rod {:s}: {:8.3f} {:8.3f} {:8d}'.format(rod_name,rod_mean,rod_std,rod_count)) # Report mean, SD, # voxels
        HU.append(rod_mean)
        
    # Perform calibration
    ogo.message('Fitting parameters:')
    
    if (phantom_dict['type'] == 'k2hpo4'):     # Mindways calibration (k2hpo4)
        densities = phantom_dict['densities']
        water = phantom_dict['h2o_densities']
        
        if len(water) != len(densities) or len(water) != len(label_list):
            os.sys.exit('[ERROR] Number of water, K2HPO4, and segmented labels are\n'
                        '        not the same: {}, {}, and {}'.format(len(water), len(densities), len(label_list)))
                        
        ogo.message('  {:>14s} = '.format('Water')+'['+', '.join("{:8.3f}".format(i) for i in water)+"]")
        ogo.message('  {:>14s} = '.format('K2HPO4')+'['+', '.join("{:8.3f}".format(i) for i in densities)+"]")
        ogo.message('  {:>14s} = '.format('HU')+'['+', '.join("{:8.3f}".format(i) for i in HU)+"]")
    
        calibrator = MindwaysCalibration()
        calibrator.fit(HU, densities, water)
    
    else:                                      # Calcium hydroxyapetite calibration (CHA)
        densities = phantom_dict['densities']
        
        if len(densities) != len(label_list):
            os.sys.exit('[ERROR] Number of CHA and segmented labels are\n'
                        '        not the same: {}, {}, and {}'.format(len(densities), len(label_list)))
                        
        ogo.message('  {:>14s} = '.format('K2HPO4')+'['+', '.join("{:8.3f}".format(i) for i in densities)+"]")
        ogo.message('  {:>14s} = '.format('HU')+'['+', '.join("{:8.3f}".format(i) for i in HU)+"]")
    
        calibrator = StandardCalibration()
        calibrator.fit(HU, densities)

    ogo.message('Found fit:')
    ogo.message('  Slope:     {:8.6f}'.format(calibrator.slope))
    ogo.message('  Intercept: {:8.6f}'.format(calibrator.intercept))
    ogo.message('  R^2:       {:8.6f}'.format(calibrator.r_value**2))
    
    ogo.message('Calibrating CT file.')
    density = calibrator.predict(sitk.Cast(ct, sitk.sitkFloat64))
    density = sitk.Cast(density, ct.GetPixelID())
    
    # Check that the calibration results are within a reasonable range
    threshold_warning = 2 # percent
    threshold_error =  10 # percent
    if (abs(1.0-(calibrator.r_value**2))*100.0>threshold_error):
        ogo.message('[ERROR] R-value is too low ({:8.6f}).'.format(calibrator.r_value**2))
        ogo.message('          * did you select the correct phantom?')
        ogo.message('          * do all the rods have a reasonable value?')
        ogo.message('          * is the StdDev of the rods high (i.e. missed the target)?')
        os.sys.exit()
    if (abs(1.0-(calibrator.r_value**2))*100.0>threshold_warning):
        ogo.message('[WARNING] R-value is below {:.2f} ({:8.6f}). Checking following:'.format((100-threshold_warning)/100.0,calibrator.r_value**2))
        ogo.message('          * did you select the correct phantom?')
        ogo.message('          * do all the rods have a reasonable value?')
        ogo.message('          * is the StdDev of the rods high (i.e. missed the target)?')
    
    ogo.message('Writing calibrated output image to file:')
    ogo.message('      \"{}\"'.format(output_image))
    sitk.WriteImage(density, output_image)

    # Write text file
    if calib_file_name:
        ogo.message('Saving calibration parameters to file:')
        ogo.message('      \"{}\"'.format(calib_file_name))

        txt_file = open(calib_file_name, "w")
        
        txt_file.write('Phantom calibration:\n')
        txt_file.write('  {:>27s} {:8s}\n'.format('---------------------------','--------'))
        txt_file.write('  {:>27s} {:s}\n'.format('ID:',os.path.basename(output_image)))
        txt_file.write('  {:>27s} {:s}\n'.format('python script:',os.path.splitext(os.path.basename(sys.argv[0]))[0]))
        txt_file.write('  {:>27s} {:.2f}\n'.format('version:',script_version))
        txt_file.write('  {:>27s} {:s}\n'.format('creation date:',str(date.today())))
        txt_file.write('\n')
        txt_file.write('Files:\n')
        txt_file.write('  {:>27s} {:8s}\n'.format('---------------------------','--------'))
        txt_file.write('  {:>14s} = {:s}\n'.format('input image',input_image))
        txt_file.write('  {:>14s} = {:s}\n'.format('input mask',input_mask))
        txt_file.write('  {:>14s} = {:s}\n'.format('output image',output_image))
        txt_file.write('  {:>14s} = {:s}\n'.format('async image',async_image))
        txt_file.write('  {:>14s} = {:s}\n'.format('calibration file name',calib_file_name))
        txt_file.write('\n')
        txt_file.write('Calibration parameters:\n')
        txt_file.write('  {:>27s} {:8s}\n'.format('---------------------------','--------'))
        txt_file.write('  {:>27s} {}\n'.format('Fit:',calibrator._is_fit))
        txt_file.write('  {:>27s} {:12.6f}\n'.format('Slope:',calibrator.slope))
        txt_file.write('  {:>27s} {:12.6f}\n'.format('Intercept:',calibrator.intercept))
        txt_file.write('  {:>27s} {:12.6f}\n'.format('R^2:',calibrator.r_value**2))
        txt_file.write('  {:>27s} {:8s}\n'.format('---------------------------','--------'))
        txt_file.write('  {:>27s}\n'.format('Density [HU]:'))
        #txt_file.write('  {:>27s} {:8.3f}'.format('Adipose ',calib.adipose_hu))
        for idx,label in enumerate(label_list):
            txt_file.write('  {:45s} = {:8.3f}\n'.format(lb.labels_dict[label]['LABEL'],filt.GetMean(int(label))))
        txt_file.write('  {:>27s} {:8s}\n'.format('---------------------------','--------'))
        txt_file.write('Phantom:\n')
        txt_file.write('  {:>14s} = {:s}\n'.format('name',phantom_dict['name']))
        txt_file.write('  {:>14s} = {:s}\n'.format('type',phantom_dict['type']))
        txt_file.write('  {:>14s} = {:s}\n'.format('serial',phantom_dict['serial']))
        txt_file.write('  {:>14s} = {:d}\n'.format('number_rods',phantom_dict['number_rods']))
        for idx,j in enumerate(phantom_dict['densities']):
            txt_file.write('  Rod {} (label {}): {:8.3f} k2hpo4 mg/cc\n'.format(phantom_dict['rod_names'][idx],phantom_dict['rod_labels'][idx],phantom_dict['densities'][idx]))
        if (phantom_dict['type'] == 'k2hpo4'):     # Mindways calibration (k2hpo4)
            txt_file.write('\n')
            for idx,j in enumerate(phantom_dict['h2o_densities']):
                txt_file.write('  Rod {} (label {}): {:8.3f} water mg/cc\n'.format(phantom_dict['rod_names'][idx],phantom_dict['rod_labels'][idx],phantom_dict['h2o_densities'][idx]))
        
        txt_file.write('  {:>27s} {:8s}\n'.format('---------------------------','--------'))
        txt_file.write('\n')
        txt_file.write('Unformatted calibration parameters:\n')
        txt_file.write('  {:>27s} {:8s}\n'.format('---------------------------','--------'))
        for label,value in calibrator.get_dict().items():
            txt_file.write('  {:>27s} {}\n'.format(label,value))
        txt_file.write('  {:>27s} {:8s}\n'.format('---------------------------','--------'))
        
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

    if not (input_mask.lower().endswith('.nii') or input_mask.lower().endswith('.nii.gz')):
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
        txt_file.write('  {:>27s} {:8.6f}\n'.format('Energy [keV]:',calib.effective_energy))
        txt_file.write('  {:>27s} {:8.6f}\n'.format('Max R^2:',calib.max_r2))
        txt_file.write('  {:>27s} {:8s}\n'.format('---------------------------','--------'))
        txt_file.write('  {:>27s}\n'.format('Density [HU]:'))
        txt_file.write('  {:>27s} {:8.3f}'.format('Adipose ',calib.adipose_hu))
        if 91 in labelList:
            txt_file.write(' [used]\n')
        else:
            txt_file.write(' [not used]\n')
        txt_file.write('  {:>27s} {:8.3f}'.format('Air ',calib.air_hu))
        if 92 in labelList:
            txt_file.write(' [used]\n')
        else:
            txt_file.write(' [not used]\n')
        txt_file.write('  {:>27s} {:8.3f}'.format('Blood ',calib.blood_hu))
        if 93 in labelList:
            txt_file.write(' [used]\n')
        else:
            txt_file.write(' [not used]\n')
        txt_file.write('  {:>27s} {:8.3f}'.format('Cortical Bone ',calib.bone_hu))
        if 94 in labelList:
            txt_file.write(' [used]\n')
        else:
            txt_file.write(' [not used]\n')
        txt_file.write('  {:>27s} {:8.3f}'.format('Skeletal Muscle ',calib.muscle_hu))
        if 95 in labelList:
            txt_file.write(' [used]\n')
        else:
            txt_file.write(' [not used]\n')
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

Please cite:
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