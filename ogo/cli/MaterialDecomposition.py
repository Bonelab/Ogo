# /------------------------------------------------------------------------------+
# | 7-NOV-2025                                                                  |
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

def print_dict(d):
  keys = d.keys()
  for k in keys:
    print('  {:10s}: {:30s}'.format('material',k))
    print('     {:>30s}: {:12.4f} {}'.format('energy',d[k]['energy'],'keV'))
    print('     {:>30s}: {:12.4f} {}'.format('mass_density',d[k]['mass_density'],'g/cm3'))
    print('     {:>30s}: {:12.4f} {}'.format('mass_concentration',d[k]['mass_concentration'],'g/cm3'))
    print('     {:>30s}: {:12.4f} {}'.format('mass_attenuation',d[k]['mass_attenuation'],'cm2/g'))
    print('     {:>30s}: {:12.4f} {}'.format('icru_mass_density',d[k]['icru_mass_density'],'g/cm3'))
    print('     {:>30s}: {:12.4f} {}'.format('mu',d[k]['mu'],'/cm'))
    print('     {:>30s}: {:12.4f} {}'.format('ctn',d[k]['ctn'],'HU'))
  
def CalculateLinearAttenuation(energy, material, mass_density, ctn_measured, header=False, delimiter=',', description='', quiet=False):
  
  # Initialize variables
  entry = []
  choices=['adipose','air','blood','bone','calcium','cha','k2hpo4','muscle','water','softtissue','redmarrow','yellowmarrow','spongiosa','iodine']
  linear_attenuation_dict = {}
  material_dict = {
    'energy':0.0,
    'mass_density':0.0,
    'mass_concentration':0.0,
    'mass_attenuation':0.0,
    'icru_mass_density':0.0,
    'mu':0.0,
    'ctn':0.0
  }
  if ctn_measured is not None:
    calculate_error = True
  else:
    calculate_error = False
  
  # Check inputs
  if energy<30 or energy>200:
    os.sys.exit('[ERROR] Energy is expected to be between 30 and 200 keV.')
  if not material:
    os.sys.exit('[ERROR] At least one ICRU material must be defined.')
  if not mass_density:
    os.sys.exit('[ERROR] At least one mass concentration must be defined.')
  if len(material) != len(mass_density):
    os.sys.exit('[ERROR] N_material = {} and N_mass_density = {}.\n'.format(len(material),len(mass_density))+
                '        Every ICRU material needs a mass concentration defined.')
  for mat in material:
    if mat not in choices:
      os.sys.exit('[ERROR] Material is not valid: {}'.format(mat))
    linear_attenuation_dict[mat] = copy.deepcopy(material_dict)                   # initialize dictionary
    linear_attenuation_dict[mat]['energy'] = energy
   
  # Determine mass_concentration
  for idx,input_mass_density in enumerate(mass_density):
    mat = material[idx]
    if input_mass_density<0:
      os.sys.exit('[ERROR] Invalid mass density. Cannot be negative: {}'.format(input_mass_density))
    icru_density = md.mass_density()[material[idx]]
    mass_concentration = (input_mass_density / icru_density)
    linear_attenuation_dict[mat]['mass_density'] = input_mass_density
    linear_attenuation_dict[mat]['icru_mass_density'] = icru_density
    linear_attenuation_dict[mat]['mass_concentration'] = mass_concentration
    if mass_concentration > 1.0:
      os.sys.exit('[ERROR] Material = {}, ICRU density = {:.3f}, mass_density = {:.3f}\n'.format(mat,icru_density,input_mass_density) +
                  '        Mass concentration {} cannot be greater than its ICRU density.'.format(mass_concentration))

  # Look up mass attenuations from NIST
  for mat in material:
    mass_attenuation = md.interpolate_mass_attenuation(mat,energy)
    linear_attenuation_dict[mat]['mass_attenuation'] = mass_attenuation
  
  # Calculate linear attenuation
  mu_final = 0.0
  for mat in material:
    #mu = linear_attenuation_dict[mat]['mass_attenuation'] * linear_attenuation_dict[mat]['mass_concentration']
    mu = linear_attenuation_dict[mat]['mass_attenuation'] * linear_attenuation_dict[mat]['mass_density']
    linear_attenuation_dict[mat]['mu'] = mu
    mu_final += mu
  
  mu_water = md.interpolate_mass_attenuation('water',energy) # we need this to convert to CT numbers
  
  # Calculate the CT number (HU)
  for mat in material:
    mu = linear_attenuation_dict[mat]['mu']
    ctn = round(md.mu_to_ct_number(mu,mu_water))
    linear_attenuation_dict[mat]['ctn'] = ctn
  
  # CT number based on total linear attenuation
  ctn_final = round(md.mu_to_ct_number(mu_final,mu_water))

  # Calculate error if a CT value is provided?
  if calculate_error:
    mu_measured = md.ct_number_to_mu(ctn_measured,mu_water)
    mu_delta =  mu_final - mu_measured
    mu_percent_error = (mu_delta)/mu_measured * 100.0
    ctn_delta =  ctn_final - ctn_measured
    if ctn_measured != 0:
      ctn_percent_error = (ctn_delta)/ctn_measured * 100.0
    
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
  Given the mass density of an ICRU-defined material and the X-ray 
  energy of the virtual monoenergetic image (VMI) calculate the expected 
  linear attenuation.
  
  If multiple materials are defined (e.g., HA and water) then the 
  calculated linear attenuation will be the summed result. 
  
  If the CT number measured in Hounsfield units is provided from the
  scan then the error relative to the theoretical value will be output.
  
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
  
  Use quiet mode for output to STDOUT.

'''
  epilog = '''
Example calls: 
  ogoCalculateLinearAttenuation 80 --material cha 0.202 --quiet
  ogoCalculateLinearAttenuation 80 --material cha water \\
                                   --mass_density 0.202 0.972
  ogoCalculateLinearAttenuation 80 --material cha water \\
                                   --mass_density 0.202 0.972 \\
                                   --ctn_measured 236 \\
                                   -d "PCD_ESP_120kV_080keV_Br40_04Th"
  
'''

  # Setup argument parsing
  parser = argparse.ArgumentParser(
      formatter_class=argparse.RawTextHelpFormatter,
      prog="ogoCalculateLinearAttenuation",
      description=description,
      epilog=epilog
  )
  parser.add_argument('energy',
                      type=float, 
                      default=50.0, 
                      metavar='Energy', 
                      help='Energy of the scan (keV)')
  parser.add_argument('--material', 
                      nargs='*', 
                      default=[], 
                      metavar='MAT1 MAT2', 
                      help='Specify ICRU materials (bone, muscle, cha, etc)')
  parser.add_argument('--mass_density', 
                      type=float, nargs='*', 
                      default=[], 
                      metavar='MASS1 MASS2', 
                      help='Specify mass density for each material (g/cm3)')
  parser.add_argument('--ctn_measured', 
                      type=float, 
                      default=None, 
                      metavar='CTN', 
                      help='Measured CT number in HU')
  parser.add_argument('--header', '-H',
                      action="store_true",
                      help="""Print a header line first (default: %(default)s)""")
  parser.add_argument ("--delimiter",
                      metavar='VAL', 
                      default = "\t",
                      help="""Delimiter character (default: tab, '\\t')""")
  parser.add_argument ("--description", "-d",
                      metavar='TXT', 
                      default = '',
                      help="""Description to pass to output (e.g., filename) (default: %(default)s)""")
  parser.add_argument("--quiet", 
                      action='store_true', 
                      help='Overwrite output without asking')
  # Parse and display
  args = parser.parse_args()
  if not args.quiet:
    print(echo_arguments('CalculateLinearAttenuation', vars(args)))

  # Run program
  CalculateLinearAttenuation(**vars(args))

if __name__ == '__main__':
    main()
