# /------------------------------------------------------------------------------+
# | 5-DEC-2025                                                                  |
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

def MaterialDecomposition(images, energies, materials, quiet, overwrite):
  
  # Initialize variables
  choices=['adipose','air','blood','bone','calcium','cha','k2hpo4','muscle','water','softtissue','redmarrow','yellowmarrow','spongiosa','iodine']
  
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

  if not quiet:
    ogo.message('Energies:          [' + ', '.join('{:.4f}'.format(i) for i in energies) + ']')
    ogo.message('Materials:         [' + ', '.join('{:s}'.format(i) for i in materials) + ']')
    ogo.message('Mass_attenuations: [' + ', '.join('{:.4f}'.format(i) for i in mass_attenuations) + ']')
    ogo.message('')
    
  # Read in CT data
  ct1 = sitk.ReadImage(images[0])
  if not quiet:
    ogo.message('Reading image 1: {:s}'.format(images[0]))
    ogo.message('      dimension: {:d}'.format(ct1.GetDimension()))
    ogo.message('           type: {:s}'.format(ct1.GetPixelIDTypeAsString()))
    ogo.message('           size: [' + ', '.join('{:d}'.format(i) for i in ct1.GetSize()) + ']')
    ogo.message('         origin: [' + ', '.join('{:.2f}'.format(i) for i in ct1.GetOrigin()) + ']')
    ogo.message('        spacing: [' + ', '.join('{:.2f}'.format(i) for i in ct1.GetSpacing()) + ']')
    ogo.message('      direction: [' + ', '.join('{:.1f}'.format(i) for i in ct1.GetDirection()) + ']')
    ogo.message('')
    
  ct2 = sitk.ReadImage(images[1])
  if not quiet:
    ogo.message('Reading image 2: {:s}'.format(images[1]))
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
  if ct1.GetOrigin() != ct2.GetOrigin():
    ogo.message('[WARNING] Origin of CT images differ.')
  if ct1.GetSpacing() != ct2.GetSpacing():
    ogo.message('[WARNING] Spacing of CT images differ.')
  if ct1.GetDirection() != ct2.GetDirection():
    ogo.message('[WARNING] Direction of CT images differ.')
    
  # Simple ITK image to numpy array
  np_ct1 = sitk.GetArrayFromImage(ct1)
  np_ct2 = sitk.GetArrayFromImage(ct2)
  
  dims = [ct1.GetSize()[2],ct1.GetSize()[1],ct1.GetSize()[0]]
  
  # Flatten
  np_ct1_flattened=np_ct1.reshape(-1)
  
  # Do material decomposition here...
  
  # Reshape
  reshaped_np_ct1 = np_ct1_flattened.reshape(dims)
  
  # Numpy array to simple ITK image
  ct_material1 = sitk.GetImageFromArray(reshaped_np_ct1)
  ct_material1.SetOrigin(ct1.GetOrigin())
  ct_material1.SetSpacing(ct1.GetSpacing())
  ct_material1.SetDirection(ct1.GetDirection())
  if not quiet:
    ogo.message('[IMPORTANT] Origin, spacing and direction are taken from image 1.')
  
  output_image = 'test.nii.gz'
  ogo.message('Writing output to file {}'.format(output_image))
  sitk.WriteImage(ct_material1, output_image)
  
  
  exit()
    
  # Print results to screen
  if not quiet:
    print('--------------------------------------------------------------------------------')
    print('Materials:')
    print('--------------------------------------------------------------------------------')
    print_dict(linear_attenuation_dict)
    print('--------------------------------------------------------------------------------')
    print('Summed material results:')
    print('--------------------------------------------------------------------------------')
    print('     {:>30s}: {:12.4f} [/cm]'.format('mu',mu_final))
    print('     {:>30s}: {:12.4f} [HU]'.format('ctn',ctn_final))
    print('     {:>30s}: {:s}'.format('description',description))
    if calculate_error:
      print('--------------------------------------------------------------------------------')
      print('Error calculations:')
      print('--------------------------------------------------------------------------------')
      print('     {:>30s}: {:12.1f} [HU]'.format('CT number',ctn_measured))
      print('     {:>30s}: {:12.4f} [HU]'.format('delta',ctn_delta))
      print('     {:>30s}: {:12.4f} [%]'.format('error',ctn_percent_error))
      print('')
      print('     {:>30s}: {:12.4f} [/cm]'.format('mu from CT number',mu_measured))
      print('     {:>30s}: {:12.4f} [/cm]'.format('delta',mu_delta))
      print('     {:>30s}: {:12.4f} [%]'.format('error',mu_percent_error))
      print('--------------------------------------------------------------------------------')
  
  # Primary outcomes only
  entry.append( ['description', '{:s}'.format(description),'TXT'] )
  entry.append( ['energy', '{:.4f}'.format(energy), '[keV]'] )
  for idx,mat in enumerate(material):
    entry.append( ['material'+'_{}'.format(idx), '{:s}'.format(mat), '[]'] )
    entry.append( ['mass_density'+'_{}'.format(idx), '{:.4f}'.format(linear_attenuation_dict[mat]['mass_density']), '[g/cm3]'] )
    entry.append( ['mass_concentration'+'_{}'.format(idx), '{:.4f}'.format(linear_attenuation_dict[mat]['mass_concentration']), '[g/cm3]'] )
    entry.append( ['mass_attenuation'+'_{}'.format(idx), '{:.4f}'.format(linear_attenuation_dict[mat]['mass_attenuation']), '[cm2/g]'] )
    entry.append( ['icru_mass_density'+'_{}'.format(idx), '{:.4f}'.format(linear_attenuation_dict[mat]['icru_mass_density']), '[g/cm3]'] )
    entry.append( ['mu'+'_{}'.format(idx), '{:.4f}'.format(linear_attenuation_dict[mat]['mu']), '[/cm]'] )
    entry.append( ['ctn'+'_{}'.format(idx), '{:.4f}'.format(linear_attenuation_dict[mat]['ctn']), '[HU]'] )
  entry.append( ['mu_final', '{:.4f}'.format(mu_final), '[/cm]'] )
  entry.append( ['ctn_final', '{:d}'.format(ctn_final), '[/cm]'] )
  if calculate_error:
    entry.append( ['ctn_measured', '{:.1f}'.format(ctn_measured), '[HU]'] )
    entry.append( ['ctn_delta', '{:.4f}'.format(ctn_delta), '[H]'] )
    if ctn_measured != 0:
      entry.append( ['ctn_percent_error', '{:.4f}'.format(ctn_percent_error), '[%]'] )  
    else:
      entry.append( ['ctn_percent_error', 'n/a', '[%]'] )  
    entry.append( ['mu_measured', '{:.4f}'.format(mu_measured), '[/cm]'] )
    entry.append( ['mu_delta', '{:.4f}'.format(mu_delta), '[/cm]'] )
    entry.append( ['mu_percent_error', '{:.4f}'.format(mu_percent_error), '[%]'] )  
  
  # Output results for stdout
  if quiet:
    # Print the output
    entry = list(zip(*entry))
    out = sys.stdout
    if header:
        out.write (delimiter.join(entry[0]) + "\n")
        out.write (delimiter.join(entry[2]) + "\n")
    out.write (delimiter.join(entry[1]) + "\n")
    
  
def main():
  description = '''
  Performs material decomposition based on a pair of CT virtual 
  monoenergetic images (VMIs). 
  
  Typical use is 2-material decomposition into water and HA or a
  3-material decomposition into HA, water and adipose.
  
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
  
'''
  epilog = '''
Example calls: 
  ogoMaterialDecomposition image1.nii.gz image2.nii.gz \\ 
                           50 80 \\
                           cha water
  ogoMaterialDecomposition im1.nii.gz im2.nii.gz \\ 
                           50 80 \\
                           cha water adipose
                
  cd /Users/skboyd/Desktop/erlangen/pcct/models/img/                         
  ogoMaterialDecomposition PCD_ESP_140kV_040keV_Br40_04Th_10Slices.nii.gz \\
                           PCD_ESP_140kV_090keV_Br40_04Th_10Slices.nii.gz \\
                           40 90 cha water
  
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
                      help='CT VMI images (*.nii.gz)')
  parser.add_argument('energies', 
                      type=float, 
                      default=[], 
                      nargs=2, 
                      metavar='keV keV', 
                      help='Energies of images (40 to 190 keV)')
  parser.add_argument('materials', 
                      nargs='*', 
                      default=[], 
                      metavar='MAT1 MAT2 (MAT3)', 
                      help='Specify 2 or 3 ICRU materials (cha, water, adipose)')
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
