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
import copy
from ogo.util.echo_arguments import echo_arguments
import ogo.util.spectral_util as md

def print_dict(d):
  keys = d.keys()
  for k in keys:
    print('  {:10s}: {:30s}'.format('material',k))
    print('     {:>30s}: {:12.4f} {}'.format('energy',d[k]['energy'],'keV'))
    print('     {:>30s}: {:12.4f} {}'.format('volume_fraction',d[k]['volume_fraction'],''))
    print('     {:>30s}: {:12.4f} {}'.format('mass_concentration',d[k]['mass_concentration'],'g/cm3'))
    print('     {:>30s}: {:12.4f} {}'.format('mass_attenuation',d[k]['mass_attenuation'],'cm2/g'))
    print('     {:>30s}: {:12.4f} {}'.format('icru_mass_density',d[k]['icru_mass_density'],'g/cm3'))
    print('     {:>30s}: {:12.4f} {}'.format('mu',d[k]['mu'],'/cm'))
    print('     {:>30s}: {:12.4f} {}'.format('ctn',d[k]['ctn'],'HU'))
  
def CalculateLinearAttenuation(energy, material, volume_fraction, ctn_measured, header=False, delimiter=',', description='', quiet=False):
  
  # Initialize variables
  entry = []
  inferred_volume_fraction = False
  choices=['adipose','air','blood','bone','calcium','cha','k2hpo4','muscle','water','softtissue','redmarrow','yellowmarrow','spongiosa','iodine']
  linear_attenuation_dict = {}
  material_dict = {
    'energy':0.0,
    'volume_fraction':0.0,
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
  if volume_fraction:
    if sum(volume_fraction)>1.0:
      os.sys.exit('[ERROR] Total volume fraction defined cannot be greater than 1.0.')
  if len(material) != len(volume_fraction):
    if (len(material) == 1) and (not volume_fraction):
      volume_fraction = [1.0] # only one material defined, so we know volume_fraction is 1.0
    elif len(material) > 1 and (len(material)-len(volume_fraction))==1:
      volume_fraction.append(1.0 - sum(volume_fraction)) # one less volume_fraction defined, so we assume its value
      inferred_volume_fraction = True
    else:
      os.sys.exit('[ERROR] N_material = {} and N_volume_fraction = {}.\n'.format(len(material),len(volume_fraction))+
                '        Every ICRU material needs a volume fraction defined.')
  for mat in material:
    if mat not in choices:
      os.sys.exit('[ERROR] Material is not valid: {}'.format(mat))
    linear_attenuation_dict[mat] = copy.deepcopy(material_dict)                   # initialize dictionary
    linear_attenuation_dict[mat]['energy'] = energy
   
  # Determine mass_concentration from known volume fractions of ICRU materials
  for idx,input_volume_fraction in enumerate(volume_fraction):
    mat = material[idx]
    if input_volume_fraction<0:
      os.sys.exit('[ERROR] Invalid volume fraction. Cannot be negative: {}'.format(input_volume_fraction))
    icru_density = md.mass_density()[material[idx]]
    mass_concentration = (input_volume_fraction * icru_density)
    linear_attenuation_dict[mat]['volume_fraction'] = input_volume_fraction
    linear_attenuation_dict[mat]['icru_mass_density'] = icru_density
    linear_attenuation_dict[mat]['mass_concentration'] = mass_concentration
    if mass_concentration > icru_density:
      os.sys.exit('[ERROR] Material = {}, ICRU density = {:.3f}, volume_fraction = {:.3f}\n'.format(mat,icru_density,input_volume_fraction) +
                  '        Mass concentration {} cannot be greater than its ICRU density.'.format(mass_concentration))

  # Look up mass attenuations from NIST
  for mat in material:
    mass_attenuation = md.interpolate_mass_attenuation(mat,energy)
    linear_attenuation_dict[mat]['mass_attenuation'] = mass_attenuation
  
  # Calculate linear attenuation
  mu_final = 0.0
  for mat in material:
    mu = linear_attenuation_dict[mat]['mass_attenuation'] * linear_attenuation_dict[mat]['icru_mass_density'] * linear_attenuation_dict[mat]['volume_fraction']
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
    if inferred_volume_fraction:
      print('     {:>44s}'.format('Note that volume fraction of this materials is inferred.'))
    print('--------------------------------------------------------------------------------')
    print('Summed material results:')
    print('--------------------------------------------------------------------------------')
    print('     {:>30s}: {:12.4f} [1]'.format('total_volume_fraction',sum(volume_fraction)))
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
    entry.append( ['volume_fraction'+'_{}'.format(idx), '{:.4f}'.format(100.0 * linear_attenuation_dict[mat]['volume_fraction']), '[%]'] )
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
  Given the volume fractions of ICRU-defined materials and the X-ray 
  energy of the virtual monoenergetic image (VMI) calculate the expected 
  linear attenuation.
  
  If multiple materials are defined (e.g., HA and water) then the 
  calculated linear attenuation sum of the parts.
  
  The total volume fractions must sum to 1. For an HA insert, the volume 
  fraction is based on its partial density of pure HA and the ICRU total
  density. Here are some examples:
  
    101.7 mg HA/cm3 insert: HA volume fraction is 101.7 / 3160.0 = 0.0322
    202.0 mg HA/cm3 insert: HA volume fraction is 202.0 / 3160.0 = 0.0639
  
  Knowing the volume fraction of the inserts from the certifications, 
  we can calculate the volume fraction of the water materials as:
  
    101.7 mg HA/cm3 insert: water volume fraction is 1 - 0.0322 = 0.9678
    202.0 mg HA/cm3 insert: water volume fraction is 1 - 0.0639 = 0.9361
    
  If the number of volume fractions is 1 less than the number of materials
  then we assume that the remaining volume fraction is for the material that
  was undefined. Usually this would be the water component of an HA phantom.
  
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
  ogoCalculateLinearAttenuation 140 --material cha 0.202 --quiet
  ogoCalculateLinearAttenuation 140 --material cha water \\
                                    --volume_fraction 0.202 0.972
  ogoCalculateLinearAttenuation 140 --material cha water \\
                                    --volume_fraction 0.202 \\
                                    --ctn_measured 135 \\
                                    -d "PCD_ESP_120kV_140keV_Br40_04Th"
  
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
  parser.add_argument('--volume_fraction', 
                      type=float, nargs='*', 
                      default=[], 
                      metavar='f1 f2', 
                      help='Specify volume fraction for each material (sum to 1)')
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
