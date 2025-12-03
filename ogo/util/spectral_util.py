# /------------------------------------------------------------------------------+
# | 28-OCT-2025                                                                  |
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

def mass_density():  # mass_density, g/cm3
  rho = {
    'adipose':9.500E-01, #9.3 to 9.7 in ICRU46
    'air':1.205E-03,
    'blood':1.060E+00,
    'bone':1.920E+00,
    'calcium':1.550E+00,
    'cha':3.160E+00,
    'k2hpo4':2.440E+00,
    'muscle':1.050E+00,
    'water':1.000E+00,
    'softtissue':1.000E+00,
    'redmarrow':1.030E+00,
    'yellowmarrow':0.980E+00,
    'spongiosa':1.180E+00,
    'iodine':4.930
    }
  return rho

def ct_number_to_mu(hu=0.0,mu_water=1.0):
  """Convert CT number in HU to linear attenuation coefficient 

  Converts the CT number measured in Hounsfield units, hu, 
  into the linear attenuation coefficient, mu. The linear
  attenuation coefficient of water is required.
  
  mu = mu_water * (hu/1000 + 1)
  """
  return mu_water*hu/1000.0 + mu_water

def mu_to_ct_number(mu=0.0,mu_water=1.0):
  """Convert linear attenuation coefficient to CT number in HU

  Converts the linear attenuation coefficient, mu, into 
  CT number measured in Hounsfield units. The linear
  attenuation coefficient of water is required.
  
  CT_number = (mu - mu_water)/mu_water * 1000
  """
  return (mu - mu_water)/mu_water*1000.0

def mu_from_mass_attenuation_coefficient(mu_rho = 1.0, mass_density = 1.0):
  """Convert the mass attenuation coefficient to mu
  
  Converts the mass attenuation coefficient, mu/rho [cm2/g] to
  linear attenuation coefficient, mu [cm^-1] by dividing by
  mass density, rho [g/cm3]
  """
  return (mu_rho * mass_density)
  
def phantom_linear_fit(x=[1,2,3],y=[6,20,32]):
  """Determine linear fit of density to Hounsfield units

  The dependent variable, x, includes the mass density of
  the phantom rods and the independent variable, y, is the 
  CT number in Hounsfield units of the rods. A minimum of
  two rods is required."""

  if (len(x) != len(y)):
    raise Exception("phantom_linear_fit: Number of dependent and independent variables must be equal.")
    
  if (len(x) <2):
    raise Exception("phantom_linear_fit: Number of rods must be two or more.")
    
  deg = 1
  [m,b] = np.polyfit(x, y, deg)

  return [m,b]
  
def interpolate_mass_attenuation(material="water",keV=90):
  """Interpolate the mass attenuation based on NIST tables
  
  The published NIST mass attenuation [cm2/g] as a function 
  of photon energy [keV] for key materials are interpolated  
  at a given photon energy."""
  
  table_name = material + "_table"
  if table_name not in mass_attenuation_tables.keys():
    raise Exception("interpolate_mass_attenuation: Material {} is not available in current NIST tables.".format(material))
  
  table = mass_attenuation_tables[table_name]
  fp = table['Mass Attenuation [cm2/g]'].to_numpy()
  xp = table['Energy [keV]'].to_numpy()

  if keV < xp.min() or keV > xp.max():
    raise Exception("interpolate_mass_attenuation: Input keV {} is not in range of NIST table.".format(keV))

  mass_attenuation = np.interp(keV,xp,fp)
    
  return mass_attenuation
