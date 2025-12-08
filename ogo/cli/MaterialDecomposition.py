# /------------------------------------------------------------------------------+
# | 5-DEC-2025                                                                   |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+

# Imports
import sys
import os
import argparse
import math
import copy
import SimpleITK as sitk
import numpy as np
from ogo.dat.MassAttenuationTables import mass_attenuation_tables
from ogo.util.echo_arguments import echo_arguments
import ogo.util.Helper as ogo
import ogo.util.spectral_util as md

# Just used to do some debugging
def print_pattern(pattern):
  basename,dirname,name,ext = md.parse_filename(pattern)
  print('{:>20s}: {:s}'.format('INPUT', pattern))
  print('{:>20s}: {:s}'.format('basename', basename))
  print('{:>20s}: {:s}'.format('dirname', dirname))
  print('{:>20s}: {:s}'.format('name', name))
  print('{:>20s}: {:s}'.format('ext', ext))
  
def MaterialDecomposition(images, energies, materials, pattern, suppress, quiet, overwrite):
  
  # Initialize variables
  choices=['adipose','air','blood','bone','calcium','cha','k2hpo4','muscle','water','softtissue','redmarrow','yellowmarrow','spongiosa','iodine']
  
  no_print_str = '--suppressed--'
  
  # Check that input images exist and are type NIfTI
  for image in images:
    if not os.path.isfile(image):
      os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(image))
    if not (image.lower().endswith('.nii') or image.lower().endswith('.nii.gz')):
      os.sys.exit('[ERROR] Input must be type NIFTI file: \"{}\"'.format(image))
  
  # Check that the energies are in a valid range and not equal
  for energy in energies:
    if energy<40 or energy>190:
      os.sys.exit('[ERROR] Energies must be between 40 and 190 keV.')
  if energies[0] == energies[1]:
    os.sys.exit('[ERROR] Energies cannot be equal.')

  # Check that the materials are valid
  if len(materials)<2 or len(materials)>3:
    os.sys.exit('[ERROR] Only 2- or 3-material decomposition is valid.')
  for material in materials:
    if material not in choices:
      os.sys.exit('[ERROR] Material is not valid: {}'.format(material))

  # Look up mass attenuations of materials
  mass_attenuations=[]
  for material in materials:
    for energy in energies:
      mass_attenuations.append(md.interpolate_mass_attenuation(material,energy))

  # Look up ICRU mass densitiesof materials
  mass_densities=[]
  for material in materials:
    mass_densities.append(md.mass_density()[material])

  if not quiet:
    ogo.message('Low energy image at {:.1f} keV:'.format(energies[0]))
    basename,dirname,name,ext = md.parse_filename(images[0])
    ogo.message('  {:s}'.format(name+ext))
    ogo.message('')
    ogo.message('High energy image at {:.1f} keV:'.format(energies[1]))
    basename,dirname,name,ext = md.parse_filename(images[1])
    ogo.message('  {:s}'.format(name+ext))
    ogo.message('')
    ogo.message('Materials:         [' + ', '.join('{:s}'.format(i) for i in materials) + ']')
    ogo.message('')
    ogo.message('Mass attenuations: [' + ', '.join('{:.4f}'.format(i) for i in mass_attenuations) + ']')
    ogo.message('')
    ogo.message('Mass densities:    [' + ', '.join('{:.4f}'.format(i) for i in mass_densities) + ']')
    ogo.message('')
    
  # Define output filenames and check for overwrite
  ofiles = []
  if pattern is not None:
    basename,dirname,name,ext = md.parse_filename(pattern)
    if not ext:
      ext='.nii.gz'
    if dirname and not os.path.isdir(os.path.dirname(pattern)):
      os.sys.exit('[ERROR] Output directory does not exist: {}'.format(os.path.dirname(pattern)))
  else:
    basename,dirname,name,ext = md.parse_filename(images[0])
  
  for material in materials:
    if material not in suppress:
      ofiles.append(name + '_' + '{}'.format(len(materials)) + material + ext)
    else:
      ofiles.append(no_print_str)

  for ofile in ofiles:
    if os.path.isfile(ofile) and not overwrite:
        result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(ofile))
        if result.lower() not in ['y', 'yes']:
            ogo.message('Not overwriting. Exiting...')
            os.sys.exit()
    
  # Read in CT data ---------------------------------------------------------------------------------
  if not quiet:
    ogo.message('Input files')
  ct1 = sitk.ReadImage(images[0])
  if not quiet:
    ogo.message('      low image: {:s}'.format(os.path.basename(images[0])))
    ogo.message('         energy: {:.1f}'.format(energies[0]))
    ogo.message('      directory: {:s}'.format(os.path.dirname(images[0])))
    ogo.message('      dimension: {:d}'.format(ct1.GetDimension()))
    ogo.message('           type: {:s}'.format(ct1.GetPixelIDTypeAsString()))
    ogo.message('           size: [' + ', '.join('{:d}'.format(i) for i in ct1.GetSize()) + ']')
    ogo.message('         origin: [' + ', '.join('{:.2f}'.format(i) for i in ct1.GetOrigin()) + ']')
    ogo.message('        spacing: [' + ', '.join('{:.2f}'.format(i) for i in ct1.GetSpacing()) + ']')
    ogo.message('      direction: [' + ', '.join('{:.1f}'.format(i) for i in ct1.GetDirection()) + ']')
    ogo.message('')
    
  ct2 = sitk.ReadImage(images[1])
  if not quiet:
    ogo.message('     high image: {:s}'.format(os.path.basename(images[1])))
    ogo.message('      directory: {:s}'.format(os.path.dirname(images[1])))
    ogo.message('      dimension: {:d}'.format(ct2.GetDimension()))
    ogo.message('           type: {:s}'.format(ct2.GetPixelIDTypeAsString()))
    ogo.message('           size: [' + ', '.join('{:d}'.format(i) for i in ct2.GetSize()) + ']')
    ogo.message('         origin: [' + ', '.join('{:.2f}'.format(i) for i in ct2.GetOrigin()) + ']')
    ogo.message('        spacing: [' + ', '.join('{:.2f}'.format(i) for i in ct2.GetSpacing()) + ']')
    ogo.message('      direction: [' + ', '.join('{:.1f}'.format(i) for i in ct2.GetDirection()) + ']')
    ogo.message('')
    
  # Check dimensions of images are the same and pixel type
  if ct1.GetDimension() != ct2.GetDimension():
    os.sys.exit('[ERROR] CT input files are not the same dimensions.')
  if ct1.GetSize() != ct2.GetSize():
    os.sys.exit('[ERROR] CT input files are not the same sizes.')
  if ct1.GetPixelID() != sitk.sitkInt16:
    ogo.message('[WARNING] Unexpected CT data type: {}'.format(ct1.GetPixelIDTypeAsString()))
  if ct2.GetPixelID() != sitk.sitkInt16:
    ogo.message('[WARNING] Unexpected CT data type: {}'.format(ct2.GetPixelIDTypeAsString()))
  if (ct1.GetOrigin() != ct2.GetOrigin() or 
      ct1.GetSpacing() != ct2.GetSpacing() or 
      ct1.GetDirection() != ct2.GetDirection()):
    ogo.message('[WARNING] Origin, Spacing or Direction are different for input images.')
    ogo.message('          Output settings are taken from image 1 input.')

  origin = ct1.GetOrigin()
  spacing = ct1.GetSpacing()
  direction = ct1.GetDirection()
  
  if not quiet:
    odir = os.path.dirname(ofiles[0])
    if not odir: odir='.'
    ogo.message('Output files')
    for idx, ofile in enumerate(ofiles):
      ogo.message('{:>15s}: {:s}'.format(materials[idx],os.path.basename(ofile)))
    ogo.message('{:>15s}: {:s}'.format('DIR',odir))
    ogo.message('')
    
  # Convert SimpleITK images to numpy arrays and flatten them to 1D
  np_low = sitk.GetArrayFromImage(ct1)
  np_low_flattened=np_low.reshape(-1)

  np_high = sitk.GetArrayFromImage(ct2)
  np_high_flattened=np_high.reshape(-1)
  
  original_shape = np_low.shape
  
  # Convert images from CT number to linear attenuation
  mu_water_low = md.interpolate_mass_attenuation('water',energies[0])
  mu_water_high = md.interpolate_mass_attenuation('water',energies[1])
  
  np_low_flattened_mu = md.ct_number_to_mu(np_low_flattened,mu_water_low)
  np_high_flattened_mu = md.ct_number_to_mu(np_high_flattened,mu_water_high)
  
  # 2-material decomposition ------------------------------------------------------------------------
  if len(materials)==2:
    ogo.message('Performing 2-material decomposition without calibration phantom.')
    ogo.message('  Materials are: [' + ', '.join('{}'.format(i) for i in materials) + ']')
    ogo.message('  No calibration phantom used.')
    ogo.message('')
    
    # AX = B, where
    # 
    # ⎡ a b ⎤⎡ rho_mat1 ⎤ = ⎡ e ⎤
    # ⎣ c d ⎦⎣ rho_mat2 ⎦   ⎣ f ⎦
    
    A_matrix = np.array(mass_attenuations).reshape(2, 2)
    B_matrix = np.vstack((np_low_flattened_mu,np_high_flattened_mu))
    
    solution = md.solve_system_equations(A_matrix,B_matrix,True)
    
    # Reshape
    np_ct_material1 = solution[0,:].reshape(original_shape)
    np_ct_material2 = solution[1,:].reshape(original_shape)
    
    # Multiply by mass densities
    np_ct_material1 = np_ct_material1 * mass_densities[0]
    np_ct_material2 = np_ct_material2 * mass_densities[1]
    
    # Create SimpleITK images
    ct_material1 = sitk.GetImageFromArray(np_ct_material1)
    ct_material1.SetOrigin(origin)
    ct_material1.SetSpacing(spacing)
    ct_material1.SetDirection(direction)
    
    ct_material2 = sitk.GetImageFromArray(np_ct_material2)
    ct_material2.SetOrigin(origin)
    ct_material2.SetSpacing(spacing)
    ct_material2.SetDirection(direction)
    
    # Write output
    ogo.message('Writing output')
    if ofiles[0] is not no_print_str:
      ogo.message('  {}'.format(ofiles[0]))
      sitk.WriteImage(ct_material1, ofiles[0])
    if ofiles[1] is not no_print_str:
      ogo.message('  {}'.format(ofiles[1]))
      sitk.WriteImage(ct_material2, ofiles[1])

  # 3-material decomposition ------------------------------------------------------------------------
  if len(materials)==3:
    ogo.message('Performing 3-material decomposition without calibration phantom.')
    ogo.message('  Materials are: [' + ', '.join('{}'.format(i) for i in materials) + ']')
    ogo.message('  No calibration phantom used.')
    ogo.message('')
    
    # AX = B, where
    # 
    # ⎡ a b c ⎤⎡ rho_mat1 ⎤   ⎡ j ⎤
    # ⎜ d e f ⎟⎜ rho_mat2 ⎟ = ⎜ k ⎟
    # ⎣ g h i ⎦⎣ rho_mat3 ⎦   ⎣ 1 ⎦
    
    A_matrix = np.array(mass_attenuations).reshape(2, 3) # Mass attenuations
    eq3 = np.array([[1.0/mass_densities[0], 1.0/mass_densities[1], 1.0/mass_densities[2]]]) # Inverse mass density
    A_matrix = np.append(A_matrix, eq3, axis=0)

    B_matrix = np.vstack((np_low_flattened_mu,np_high_flattened_mu,np.ones(np_low_flattened_mu.size))) # Add 1's for third equation
    
    solution = md.solve_system_equations(A_matrix,B_matrix,True)
    
    # Reshape
    np_ct_material1 = solution[0,:].reshape(original_shape)
    np_ct_material2 = solution[1,:].reshape(original_shape)
    np_ct_material3 = solution[2,:].reshape(original_shape)
    
    # Multiply by mass densities
    np_ct_material1 = np_ct_material1 * mass_densities[0]
    np_ct_material2 = np_ct_material2 * mass_densities[1]
    np_ct_material3 = np_ct_material3 * mass_densities[2]
    
    # Create SimpleITK images
    ct_material1 = sitk.GetImageFromArray(np_ct_material1)
    ct_material1.SetOrigin(origin)
    ct_material1.SetSpacing(spacing)
    ct_material1.SetDirection(direction)
    
    ct_material2 = sitk.GetImageFromArray(np_ct_material2)
    ct_material2.SetOrigin(origin)
    ct_material2.SetSpacing(spacing)
    ct_material2.SetDirection(direction)
    
    ct_material3 = sitk.GetImageFromArray(np_ct_material3)
    ct_material3.SetOrigin(origin)
    ct_material3.SetSpacing(spacing)
    ct_material3.SetDirection(direction)
    
    # Write output
    ogo.message('Writing output')
    if ofiles[0] is not no_print_str:
      ogo.message('  {}'.format(ofiles[0]))
      sitk.WriteImage(ct_material1, ofiles[0])
    if ofiles[1] is not no_print_str:
      ogo.message('  {}'.format(ofiles[1]))
      sitk.WriteImage(ct_material2, ofiles[1])
    if ofiles[2] is not no_print_str:
      ogo.message('  {}'.format(ofiles[2]))
      sitk.WriteImage(ct_material2, ofiles[2])
  
  ogo.message('')
  ogo.message('Done.')
  
  
def main():
  description = '''
  Performs image-space material decomposition based on a pair of 
  images representing LOW and HIGH acquisitions. These can be 
  from VMIs or any other LOW/HIGH pair.
  
  Typical use is 2-material decomposition into water and HA or a
  3-material decomposition into HA, water and adipose. However
  any combinations can be employed, although do this at your own
  risk.
  
  Input images are expected to be CT numbers in Hounsfield units.
  
  Valid ICRU materials are:
    adipose
    air
    blood
    bone
    calcium
    cha
    k2hpo4
    muscle
    water
    softtissue
    redmarrow
    yellowmarrow
    spongiosa
    iodine
    
  An output image is generated for two or three materials 
  depending on type type of composition. The default output
  filename has the suffix decomposition type (2 or 3) and 
  material added to it. For example, _2cha, _2water or possibly
  _3cha, _3water, _3adipose.

  Use the PATTERN option if a new directory and base filename 
  is preferred, otherwise output is written to same directory
  as input images and same base filename with suffix. 
  
  Suppress writing an output image by specifying which 
  material(s).
  
'''
  epilog = '''
Example calls: 
  ogoMaterialDecomposition im1.nii im2.nii 40 80 cha water
  
  ogoMaterialDecomposition \\
    PCD_HIPQC_120kV_060KEV_LOW_Qr40_04Th_10Slices.nii.gz \\
    PCD_HIPQC_120kV_084KEV_HGH_Qr40_04Th_10Slices.nii.gz \\
    60 84 cha water adipose --pattern PCD_HIPQC_120kV.nii.gz
  
  ogoMaterialDecomposition \\
    PCD_HIPQC_120kV_060KEV_LOW_Qr40_04Th_10Slices.nii.gz \\
    PCD_HIPQC_120kV_084KEV_HGH_Qr40_04Th_10Slices.nii.gz \\
    60 84 cha water --pattern PCD_HIPQC_120kV.nii.gz
  
  ogoMaterialDecomposition \\
    PCD_ESP_140kV_040keV_Br40_04Th_10Slices.nii.gz \\
    PCD_ESP_140kV_090keV_Br40_04Th_10Slices.nii.gz \\
    40 90 cha water --pattern ./PCD_ESP_140kV.nii.gz
  
  ogoMaterialDecomposition \\
    PCD_ESP_140kV_040keV_Br40_04Th_10Slices.nii.gz \\
    PCD_ESP_140kV_090keV_Br40_04Th_10Slices.nii.gz \\
    40 90 cha water adipose --pattern ./PCD_ESP_140kV.nii.gz
  
'''

  # Setup argument parsing
  parser = argparse.ArgumentParser(
      formatter_class=argparse.RawTextHelpFormatter,
      prog="ogoMaterialDecomposition",
      description=description,
      epilog=epilog
  )
  parser.add_argument('images',
                      metavar='Filename Filename', 
                      nargs=2, 
                      help='CT low and high images (*.nii.gz)')
  parser.add_argument('energies', 
                      type=float, 
                      default=[], 
                      nargs=2, 
                      metavar='keV keV', 
                      help='Energies of low and high images (40 to 190 keV)')
  parser.add_argument('materials', 
                      nargs='*', 
                      default=[], 
                      metavar='MAT1 MAT2 (MAT3)', 
                      help='Specify 2 or 3 ICRU materials (eg. cha, water, adipose)')
  parser.add_argument("--pattern", 
                      default=None, 
                      metavar='FILE', 
                      help='Directory and filename pattern.')
  parser.add_argument('--suppress', 
                      nargs='*', 
                      default=[], 
                      metavar='MAT1 MAT2 (MAT3)', 
                      help='File outputs to suppress (e.g., HA, water, adipose)')
  parser.add_argument("--quiet", 
                      action='store_true', 
                      help='Overwrite output without asking')
  parser.add_argument("--overwrite", 
                      action='store_true',
                      default=False,
                      help='Overwrite output files without asking (default: %(default)s)')
  # Parse and display
  args = parser.parse_args()
  if not args.quiet:
    print(echo_arguments('MaterialDecomposition', vars(args)))

  # Run program
  MaterialDecomposition(**vars(args))

if __name__ == '__main__':
    main()
