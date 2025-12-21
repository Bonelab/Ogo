# /------------------------------------------------------------------------------+
# | 28-OCT-2025                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+

# Imports
import os
import numpy as np
from ogo.dat.MassAttenuationTables import mass_attenuation_tables
from scipy.interpolate import CubicSpline

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
  """Interpolate the mass attenuation based on NIST tables.
  
  Cubic interpolation is implemented. Minimum keV is 30.
  
  The published NIST mass attenuation [cm2/g] as a function 
  of photon energy [keV] for key materials are interpolated  
  at a given photon energy."""
  
  table_name = material + "_table"
  if table_name not in mass_attenuation_tables.keys():
    raise Exception("interpolate_mass_attenuation: Material {} is not available in current NIST tables.".format(material))
  
  table = mass_attenuation_tables[table_name]
  
  # Diagnostic range is > 30 keV 
  table = table[table['Energy [keV]'] >= 30]

  fp = table['Mass Attenuation [cm2/g]'].to_numpy()
  xp = table['Energy [keV]'].to_numpy()

  if keV < xp.min() or keV > xp.max():
    raise Exception("interpolate_mass_attenuation: Input keV {} is not in range of NIST table.".format(keV))

  #mass_attenuation = np.interp(keV,xp,fp) # linear interpolation
  spl = CubicSpline(xp,fp)
  mass_attenuation = spl(keV)
  
  return mass_attenuation

def solve_system_equations(A=np.ones((2,2)), B=np.ones(2), use_numpy=False):
  """Solve a system of equations
  
  This is used for material decomposition in spectral 
  imaging.
  
  It can produce a solution using two possible methods:
  
  1. It uses the Numpy linear algebra solver.
  
  2. Uses Cramer's Rule to solve a system of equations. 
  
  Using Cramer's rule only is only valid for either 2x2 
  or 3x3 systems of equations. Its results are less
  stable if the determinant of A is near zero.
  
  
  
  AX = B, where
  
  ⎡ a b ⎤⎡ x ⎤ = ⎡ e ⎤
  ⎣ c d ⎦⎣ y ⎦   ⎣ f ⎦
  
  ⎡ a b c ⎤⎡ x ⎤   ⎡ j ⎤
  ⎜ d e f ⎟⎜ y ⎟ = ⎜ k ⎟
  ⎣ g h i ⎦⎣ z ⎦   ⎣ l ⎦
  
  """
  
  if use_numpy:
    a = np.array(A)
    b = np.array(B)
    X = np.linalg.solve(a, b)
    return X
  
  size = len(A)

  if size < 2 or size > 3:
    raise ValueError("Matrix size must be 2 or 3.")
    
  if size == 2:
    a = A[0,0]
    b = A[0,1]
    c = A[1,0]
    d = A[1,1]
    e = B[0]
    f = B[1]
    D = a*d - b*c
    Dx = d*e - f*b
    Dy = a*f - c*e
    
    x = Dx / D
    y = Dy / D
    
    return np.array([x,y])
    
  if size == 3:
    a = A[0,0]
    b = A[0,1]
    c = A[0,2]
    d = A[1,0]
    e = A[1,1]
    f = A[1,2]
    g = A[2,0]
    h = A[2,1]
    i = A[2,2]
    j = B[0]
    k = B[1]
    l = B[2]
    D =  a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g)
    #if D==0:
    #  print('ERROR: Determinant is zero.')
    #  return np.zeros(3)
      
#    D =  np.linalg.det(A)
    Dx = j*(e*i - f*h) - b*(k*i - f*l) + c*(k*h - e*l)
    Dy = a*(k*i - f*l) - j*(d*i - f*g) + c*(d*l - g*k)
    Dz = a*(e*l - h*k) - b*(d*l - g*k) + j*(d*h - e*g)

    x = Dx / D
    y = Dy / D
    z = Dz / D
    
    return np.array([x,y,z])

def parse_filename(filename):
  """Parses a filename into parts

  Useful for extracting directory, filename,
  full filename without extension, and extension."""

  basename = os.path.basename(filename) # remove
  dirname = os.path.dirname(filename) # remove
  name, ext = os.path.splitext(filename)
  if 'gz' in ext:
      name = os.path.splitext(name)[0]  # Manages files with double extension
      ext = '.nii' + ext
  
  return basename,dirname,name,ext

def print_matrix(matrix):
  for row in matrix:
      for element in row:
          print(f"{element:4d}", end=" ") # Adjust the '4d' for desired spacing
      print() # Newline after each row

def print_matrix(matrix):
  for row in matrix:
      if len(row)>6:
        row = row[:6]
      for element in row:
          print(f"{element:4.4e}", end=" ") # Adjust the '4d' for desired spacing
      print() # Newline after each row
  