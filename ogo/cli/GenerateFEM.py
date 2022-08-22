# /------------------------------------------------------------------------------+
# | 22-AUG-2022                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+
#
script_version = 1.0

##
# Import the required modules
import os
import sys
import argparse
import vtk
import vtkbone
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

# L4 VERTEBRA ---------------------------------------------------------------------------
def vertebra(input_image,input_mask,output_model,vertebra,iso_resolution,overwrite,func):
    ogo.message('Generate an FE model of the '+vertebra+' vertebra.')

    # Check if output exists and should overwrite
    if os.path.isfile(output_model) and not overwrite:
        result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(output_model))
        if result.lower() not in ['y', 'yes']:
            ogo.message('Not overwriting. Exiting...')
            os.sys.exit()
    
    if not output_model.lower().endswith('.n88model'):
        os.sys.exit('[ERROR] Input must be .n88model file: \"{}\"'.format(output_model))
    
    # Read input image
    if not os.path.isfile(input_image):
        os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(input_image))

    if not (input_image.lower().endswith('.nii') or input_image.lower().endswith('.nii.gz')):
        os.sys.exit('[ERROR] Input must be type NIFTI file: \"{}\"'.format(input_image))

    ogo.message('Reading calibrated input CT image:')
    ogo.message('      \"{}\"'.format(input_image))
    ct = sitk.ReadImage(input_image)
    
    # Read input mask
    if not os.path.isfile(input_mask):
        os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(input_mask))

    if not (input_mask.lower().endswith('.nii') or input_mask.lower().endswith('.nii.gz')):
        os.sys.exit('[ERROR] Input must be type NIFTI file: \"{}\"'.format(input_mask))

    ogo.message('Reading input mask:')
    ogo.message('      \"{}\"'.format(input_mask))
    mask = sitk.ReadImage(input_mask)

    # Check for the appropriate label in mask
    filt = sitk.LabelStatisticsImageFilter()
    filt.Execute(ct, mask)

    vertebra_label = ogo.get_label(vertebra)
    if vertebra_label not in filt.GetLabels():
        os.sys.exit('[ERROR] Mask does not contain label {} for \"{}\" vertebra.'.format(vertebra_label,vertebra))
        
    # Resample to isotropic voxel size
    ogo.message('Resample mask.')
    mask_iso = ogo.isotropicResampling(mask,iso_resolution,'mask')
    ogo.message('Resample ct.')
    ct_iso = ogo.isotropicResampling(ct,iso_resolution,'ct')
    
    # Write output n88model file
    ogo.message('Writing n88model file:')
    ogo.message('      \"{}\"'.format(output_model))
    writer = vtkbone.vtkboneN88ModelWriter()
    #writer.SetInputData(model)
    writer.SetFileName(output_model)
    #writer.Update()
    
    ogo.message('Done generating finite element file.')
    
# FEMUR ---------------------------------------------------------------------------------
def femur(input_image,input_mask,output_model,side,iso_resolution,overwrite,func):
    ogo.message('Generate an FE model of the '+side+' femur.')

    # Check if output exists and should overwrite
    if os.path.isfile(output_model) and not overwrite:
        result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(output_model))
        if result.lower() not in ['y', 'yes']:
            ogo.message('Not overwriting. Exiting...')
            os.sys.exit()

    if not output_model.lower().endswith('.n88model'):
        os.sys.exit('[ERROR] Input must be .n88model file: \"{}\"'.format(output_model))
    
    # Read input image
    if not os.path.isfile(input_image):
        os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(input_image))

    if not (input_image.lower().endswith('.nii') or input_image.lower().endswith('.nii.gz')):
        os.sys.exit('[ERROR] Input must be type NIFTI file: \"{}\"'.format(input_image))

    ogo.message('Reading calibrated input CT image:')
    ogo.message('      \"{}\"'.format(input_image))
    ct = sitk.ReadImage(input_image)
    
    # Read input mask
    if not os.path.isfile(input_mask):
        os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(input_mask))

    if not (input_mask.lower().endswith('.nii') or input_mask.lower().endswith('.nii.gz')):
        os.sys.exit('[ERROR] Input must be type NIFTI file: \"{}\"'.format(input_mask))

    ogo.message('Reading input mask:')
    ogo.message('      \"{}\"'.format(input_mask))
    mask = sitk.ReadImage(input_mask)

    # Check for the appropriate label in mask
    filt = sitk.LabelStatisticsImageFilter()
    filt.Execute(ct, mask)

    if side in 'left':
        femur = 'Femur Left'
    else:
        femur = 'Femur Right'
    femur_label = ogo.get_label(femur)
    if femur_label not in filt.GetLabels():
        os.sys.exit('[ERROR] Mask does not contain label {} for \"{}\".'.format(femur_label,femur))
        
    # Resample to isotropic voxel size
    ogo.message('Resample mask.')
    mask_iso = ogo.isotropicResampling(mask,iso_resolution,'mask')
    ogo.message('Resample ct.')
    ct_iso = ogo.isotropicResampling(ct,iso_resolution,'ct')
    
    #sitk.WriteImage(mask_iso, '/Users/skboyd/Desktop/mask.nii.gz')
    #sitk.WriteImage(ct_iso, '/Users/skboyd/Desktop/ct.nii')
    
    # Write output n88model file
    ogo.message('Writing n88model file:')
    ogo.message('      \"{}\"'.format(output_model))
    writer = vtkbone.vtkboneN88ModelWriter()
    #writer.SetInputData(model)
    writer.SetFileName(output_model)
    #writer.Update()
    
    ogo.message('Done generating finite element file.')

    #ogo.message('  {:>27s} {:8s}'.format('---------------------------','--------'))
    #ogo.message('  {:>27s} {:8.3f}'.format('Energy [keV]:',calib.effective_energy))
    #ogo.message('  {:>27s} {:8.3f}'.format('Max R^2:',calib.max_r2))
    #ogo.message('  {:>27s} {:8s}'.format('---------------------------','--------'))
    #ogo.message('  {:>27s}'.format('Mass Attenuation [cm2/g]:'))
    #ogo.message('  {:>27s} {:8.3f}'.format('Adipose ',calib.adipose_mass_attenuation))
    #ogo.message('  {:>27s} {:8.3f}'.format('Air ',calib.air_mass_attenuation))
    #ogo.message('  {:>27s} {:8.3f}'.format('Blood ',calib.blood_mass_attenuation))
    #ogo.message('  {:>27s} {:8.3f}'.format('Cortical Bone ',calib.bone_mass_attenuation))
    #ogo.message('  {:>27s} {:8.3f}'.format('Skeletal Muscle ',calib.muscle_mass_attenuation))
    #ogo.message('  {:>27s} {:8s}'.format('---------------------------','--------'))

def main():
    # Setup description
    description='''
Uses a calibrated CT scan and mask to generate an input file for finite 
element modelling. There are two main types of FE model supported: an L4
vertebral body or a femur (left or right).

Output is an .n88model file suitable for solving using FAIM software. 
For a free copy of FAIM software follow the installation instructions at:

https://bonelab.github.io/n88
 
Please cite:
Michalski AS, Besler BA, Burt LA, Boyd SK, 2021. Opportunistic CT screening 
predicts individuals at risk of major osteoporotic fracture. Osteoporos Int 
32, 1639-1649.

Michalski AS, Besler BA, Michalak GJ, Boyd SK, 2020. CT-based internal density 
calibration for opportunistic skeletal assessment using abdominal CT scans. 
Med Eng Phys 78, 55-63.

'''

    epilog='''
Example calls: 
ogoGenerateFEM vertebra image.nii.gz mask.nii.gz vertL4.n88model
ogoGenerateFEM femur image.nii.gz mask.nii.gz --side left femurL.n88model

'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoGenerateFEM",
        description=description,
        epilog=epilog
    )
    subparsers = parser.add_subparsers()
    
    # Vertebra
    parser_vertebra = subparsers.add_parser('vertebra')
    parser_vertebra.add_argument('input_image', help='Input image file (*.nii, *.nii.gz)')
    parser_vertebra.add_argument('input_mask', help='Image mask of rods for input image or asynchronous image (*.nii, *.nii.gz)')
    parser_vertebra.add_argument('output_model', help='Output N88 model file (*.n88model)')
    parser_vertebra.add_argument('--vertebra', default='L4', choices=['L1','L2','L3','L4','L5'], help='Specify vertebra for analysis (default: %(default)s)')
    parser_vertebra.add_argument("--iso_resolution", type = float, default = 1.0, help = "Set the isotropic voxel size [mm]. (default: %(default)s [mm])")
    parser_vertebra.add_argument('--overwrite', action='store_true', help='Overwrite output without asking')
    parser_vertebra.set_defaults(func=vertebra)

    # Femur
    parser_femur = subparsers.add_parser('femur')
    parser_femur.add_argument('input_image', help='Input image file (*.nii, *.nii.gz)')
    parser_femur.add_argument('input_mask', help='Input image mask file (*.nii, *.nii.gz)')
    parser_femur.add_argument('output_model', help='Output N88 model file (*.n88model)')
    parser_femur.add_argument('--side', default='left', choices=['left','right'], help='Specify left or right femur (default: %(default)s)')
    parser_femur.add_argument("--iso_resolution", type = float, default = 1.0, help = "Set the isotropic voxel size [mm]. (default: %(default)s [mm])")
    parser_femur.add_argument('--overwrite', action='store_true', help='Overwrite output without asking')
    parser_femur.set_defaults(func=femur)
    
    # Parse and display
    args = parser.parse_args()
    print(echo_arguments('GenerateFEM', vars(args)))
    
    # Run program
    args.func(**vars(args))

if __name__ == '__main__':
    main()