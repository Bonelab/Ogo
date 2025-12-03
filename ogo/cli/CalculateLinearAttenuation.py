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
import SimpleITK as sitk
import numpy as np
from ogo.dat.MassAttenuationTables import mass_attenuation_tables
from ogo.util.echo_arguments import echo_arguments
import ogo.util.Helper as ogo
import ogo.util.spectral_util as md

def CalculateLinearAttenuation(energy, mass_concentration, material, measured_ctn, header, delimiter, description, quiet):
  
  entry = []
  
  # Check that input values are valid
  if energy<30 or energy>200:
    os.sys.exit('[ERROR] Energy must be between 30 and 200 keV.')
  if mass_concentration<=0:
    os.sys.exit('[ERROR] Mass concentration must be positive.')
  density = md.mass_density()[material]
  if mass_concentration > density and (not quiet):
    ogo.message('[WARNING] Mass concentration greater than density not physically possible.')
  
  # Calculate error if a CT value is provided?
  if measured_ctn is not None:
    calculate_error = True
  else:
    calculate_error = False
    
  # Look up mass attenuations from NIST of our material and reference material water
  mass_attenuation_of_material = md.interpolate_mass_attenuation(material,energy)
  mass_attenuation_of_water = md.interpolate_mass_attenuation('water',energy)
  
  mass_density_water = md.mass_density()['water'] # from NIST; rho_water = 1.0 g/cm3
  mass_concentration_water = md.mass_density()['water'] # from NIST; rho_water = 1.0 g/cm3
  
  # Calculate the linear attenuations knowing mass concentration
  mu_material = mass_attenuation_of_material * mass_concentration
  mu_water = mass_attenuation_of_water * mass_concentration_water
  
  # Calculate the CT number (HU)
  ctn_material = round(md.mu_to_ct_number(mu_material,mu_water))
  ctn_water = round(md.mu_to_ct_number(mu_water,mu_water))
  
  # Calculate total linear attenuation and CT number
  mu_final = mu_material + mu_water
  ctn_final = round(md.mu_to_ct_number(mu_final,mu_water))
  
  # Get the calculated mu (and CT number) 
  #mu_final = mass_attenuation_of_material * mass_concentration + mu_water * md.mass_density()['water'] * 0.97
  #mu_final = mass_attenuation_of_material * mass_concentration
  #ctn_final = round(md.mu_to_ct_number(mu_final,mu_water))

  # Calculate error
  if measured_ctn is not None:
    measured_mu = md.ct_number_to_mu(measured_ctn,mu_water)
    mu_delta =  mu_final - measured_mu
    mu_error = (mu_delta)/measured_mu * 100.0
    ct_delta =  ctn_final - measured_ctn
    if measured_ctn != 0:
      ct_error = (ct_delta)/measured_ctn * 100.0
    
  # Gather results
  if not quiet:
    ogo.message('Input:')
    ogo.message('  {:30s} = {:12.4f} [keV]'.format("energy",energy))
    ogo.message('  {:30s} = {:12.4f} [g/cm3]'.format("mass_concentration",mass_concentration))
    ogo.message('  {:30s} = {:>12s}'.format("material",material))
    ogo.message('Lookup and calculations:')
    ogo.message('  {:30s} = {:12.4f} [/cm]'.format('mass attenuation of {}'.format(material),mass_attenuation_of_material))
    ogo.message('  {:30s} = {:12.4f} [/cm]'.format('mass attenuation of {}'.format('water'),mass_attenuation_of_water))
    ogo.message('  {:30s} = {:12.4f} [g/cm3]'.format('mass density of {}'.format(material),density))
    ogo.message('  {:30s} = {:12.4f} [g/cm3]'.format('mass density of {}'.format('water'),mass_density_water))
    ogo.message('  {:30s} = {:12.4f} [g/cm3]'.format('mass concentration of {}'.format(material),mass_concentration))
    ogo.message('  {:30s} = {:12.4f} [g/cm3]'.format('mass concentration of {}'.format('water'),mass_concentration_water))
    ogo.message('  {:30s} = {:12.4f} [/cm]'.format('mu of {:s} at {} keV'.format(material,energy),mu_material))
    ogo.message('  {:30s} = {:12.4f} [/cm]'.format('mu of {:s} at {} keV'.format('water',energy),mu_water))
    ogo.message('  {:30s} = {:12.4f} [/cm]'.format('CTn of {:s} at {} keV'.format(material,energy),ctn_material))
    ogo.message('  {:30s} = {:12.4f} [/cm]'.format('CTn of {:s} at {} keV'.format('water',energy),ctn_water))
    ogo.message('Results:')
    ogo.message('  {:30s} = {:12.4f} [/cm]'.format('mu_final ({:s})'.format(material),mu_final))
    ogo.message('  {:30s} = {:12d} [HU]'.format('ctn_final ({:s})'.format(material),ctn_final))
    if calculate_error:
      ogo.message('Comparison to observed:')
      ogo.message('  {:30s} = {:12.1f} [HU]'.format('measured_ctn',measured_ctn))
      ogo.message('  {:30s} = {:12.4f} [/cm]'.format('converted to mu',measured_mu))
      ogo.message('  {:30s} = {:12.4f} [/cm]'.format('mu_delta',mu_delta))
      ogo.message('  {:30s} = {:12.4f} [HU]'.format('ct_delta',ct_delta))
      ogo.message('  {:30s} = {:12.4f} %'.format('mu_percent_error',mu_error))
      if measured_ctn != 0:
        ogo.message('  {:30s} = {:12.4f} %'.format('ct_percent_error',ct_error))
      else:
        ogo.message('  {:30s} = {:>12s}'.format('ct_percent_error','n/a'))
    if description is not None:
      ogo.message('  {:30s} = {:>12s}'.format('description',description))
  
  # Primary outcomes only
  if description is not None:
    entry.append( ['description', '{:s}'.format(description),'TXT'] )
  entry.append( ['energy', '{:.4f}'.format(energy), '[keV]'] )
  entry.append( ['mass_concentration', '{:.4f}'.format(mass_concentration), '[g/cm3]'] )
  entry.append( ['material', '{:s}'.format(material), '[name]'] )
  entry.append( ['mu_final', '{:.4f}'.format(mu_final), '[/cm]'] )
  entry.append( ['ctn_final', '{:d}'.format(ctn_final), '[/cm]'] )
  if calculate_error:
    entry.append( ['measured_ctn', '{:.1f}'.format(measured_ctn), '[HU]'] )
    entry.append( ['measured_mu', '{:.4f}'.format(measured_mu), '[/cm]'] )
    entry.append( ['mu_delta', '{:.4f}'.format(mu_delta), '[/cm]'] )
    entry.append( ['ct_delta', '{:.4f}'.format(ct_delta), '[H]'] )
    entry.append( ['mu_percent_error', '{:.4f}'.format(mu_error), '[%]'] )  
    if measured_ctn != 0:
      entry.append( ['ct_percent_error', '{:.4f}'.format(ct_error), '[%]'] )  
    else:
      entry.append( ['ct_percent_error', 'n/a', '[%]'] )  
  
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
  Given the mass concentration of an ICRU-defined material and the X-ray
  energy of the virtual monoenergetic image (VMI) calculate the expected
  linear attenuation.
  
  If the linear attenuation measured from a scan is provided then the 
  error relative to the calculated value will be provided.
  
  Use quiet mode for output to STDOUT.

'''
  epilog = '''
Example calls: 
  ogoCalculateLinearAttenuation cha 50.0 0.200
  ogoCalculateLinearAttenuation cha 100 0.202 --measured_ctn 186
  ogoCalculateLinearAttenuation water 90 1.0 --measured_ctn 0
  ogoCalculateLinearAttenuation cha 100 0.202 --measured_ctn 542 -d "PCD_ESP_120kV_080keV_Br40_04Th"
  
'''

  # Setup argument parsing
  parser = argparse.ArgumentParser(
      formatter_class=argparse.RawTextHelpFormatter,
      prog="ogoCalculateLinearAttenuation",
      description=description,
      epilog=epilog
  )
  parser.add_argument('material',
                      default='bone',
                      metavar='Material', 
                      choices=[
                        'adipose',
                        'air',
                        'blood',
                        'bone',
                        'calcium',
                        'cha',
                        'k2hpo4',
                        'muscle',
                        'water',
                        'softtissue',
                        'redmarrow',
                        'yellowmarrow',
                        'spongiosa',
                        'iodine'
                      ],
                      help='Specify ICRU material (bone, muscle, cha, etc)')
  parser.add_argument('energy',
                      type=float, 
                      default=50.0, 
                      metavar='Energy', 
                      help='Energy of the scan (keV)')
  parser.add_argument('mass_concentration',
                      type=float, 
                      default=1.0, 
                      metavar='Concentration', 
                      help='Mass concentration (g/cm3)')
  parser.add_argument('--measured_ctn', 
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
                      default = None,
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
