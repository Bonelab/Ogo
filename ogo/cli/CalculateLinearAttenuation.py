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
import dect_util as md

def CalculateLinearAttenuation(energy, mass_concentration, material, ct_number, header, delimiter, quiet):
  
  entry = []
  
  # Check that input values are valid
  if energy<30 or energy>200:
    os.sys.exit('[ERROR] Energy must be between 30 and 200 keV.')
  if mass_concentration<=0:
    os.sys.exit('[ERROR] Mass concentration must be positive.')
  mass_dens = md.mass_density()[material]
  if mass_concentration > mass_dens:
    os.sys.exit('[ERROR] Mass concentration cannot be greater than mass density.')
  
  # Calculate error if a CT value is provided?
  if ct_number:
    calculate_error = True
  else:
    calculate_error = False
    
  # Get the calculated mu (and CT number) 
  ma_material = md.interpolate_mass_attenuation(material,energy) # from NIST
  mu_water = md.interpolate_mass_attenuation('water',energy) * md.mass_density()['water']  # from NIST; rho_water = 1.0 g/cm3, so mu = ma * rho
  mu_calc = ma_material * mass_concentration + mu_water * md.mass_density()['water'] * 0.97
  #mu_calc = ma_material * mass_concentration
  ct_calc = round(md.mu_to_ct_number(mu_calc,mu_water))

  # Calculate error
  if ct_number:
    mu_meas = md.ct_number_to_mu(ct_number,mu_water)
    mu_delta =  mu_calc - mu_meas
    mu_error = (mu_delta)/mu_meas * 100.0
    ct_delta =  ct_calc - ct_number
  
  # Gather results
  if not quiet:
    ogo.message('Input:')
    ogo.message('  energy     = {:.4f} [keV]'.format(energy))
    ogo.message('  mass_concentration  = {:.4f} [g/cm3]'.format(mass_concentration))
    ogo.message('  material   = {:s}'.format(material))
    ogo.message('Lookup:')
    ogo.message('  mass_dens of {:s} = {:.4f} [g/cm3]'.format(material,mass_dens))
    ogo.message('  mu (water) at {} keV = {:.4f} [/cm]'.format(energy,mu_water))
    ogo.message('  ma ({:s}) at {} keV = {:.4f} [/cm]'.format(material,energy,ma_material))
    ogo.message('Calculated:')
    ogo.message('  mu_calc ({:s})  = {:.4f} [/cm]'.format(material,mu_calc))
    ogo.message('  ct_calc ({:s})  = {:d} [HU]'.format(material,ct_calc))
    if calculate_error:
      ogo.message('Comparison to observed:')
      ogo.message('  ct_number   = {:.1f} [HU]'.format(ct_number))
      ogo.message('  mu_meas   = {:.4f} [/cm]'.format(mu_meas))
      ogo.message('  mu_delta  = {:.4f} [/cm]'.format(mu_delta))
      ogo.message('  ct_delta  = {:.4f} [HU]'.format(ct_delta))
      ogo.message('  mu_error  = {:.4f} %'.format(mu_error))
  
  entry.append( ['energy', '{:.4f}'.format(energy), '[keV]'] )
  entry.append( ['mass_concentration', '{:.4f}'.format(mass_concentration), '[g/cm3]'] )
  entry.append( ['material', '{:s}'.format(material), '[name]'] )
  entry.append( ['mass_dens', '{:.4f}'.format(mass_dens), '[g/cm3]'] )
  entry.append( ['mu_water', '{:.4f}'.format(mu_water), '[/cm]'] )
  entry.append( ['ma_{:s}'.format(material), '{:.4f}'.format(ma_material), '[/cm]'] )
  entry.append( ['mu_calc', '{:.4f}'.format(mu_calc), '[/cm]'] )
  entry.append( ['ct_calc', '{:d}'.format(ct_calc), '[/cm]'] )
  if calculate_error:
    entry.append( ['ct_number', '{:.1f}'.format(ct_number), '[HU]'] )
    entry.append( ['mu_meas', '{:.4f}'.format(mu_meas), '[/cm]'] )
    entry.append( ['mu_delta', '{:.4f}'.format(mu_delta), '[/cm]'] )
    entry.append( ['ct_delta', '{:.4f}'.format(ct_delta), '[H]'] )
    entry.append( ['mu_error', '{:.4f}'.format(mu_error), '[%]'] )  
  
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
  python calculate_linear_attenuation.py --kev 50 --mass_concentration 1.0 --material water
  
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
  parser.add_argument('--ct_number', 
                      type=float, 
                      default=None, 
                      metavar='CTN', 
                      help='Measured CT number in HU')
  parser.add_argument('--header', '-H',
                      action="store_true",
                      help="""Print a header line first (default: %(default)s)""")
  parser.add_argument ("--delimiter", "-d",
                      metavar='VAL', 
                      default = "\t",
                      help="""Delimiter character (default: tab, '\\t').""")
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
