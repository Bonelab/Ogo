# /------------------------------------------------------------------------------+
# | 28-OCT-2025                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+

import os
import numpy as np
from ogo.dat.MassAttenuationTables import mass_attenuation_tables
from scipy.interpolate import CubicSpline
from scipy.interpolate import UnivariateSpline


def mass_density():  # mass_density, g/cm3
  rho = {
    "adipose": 9.500E-01,  # 9.3 to 9.7 in ICRU
    "air": 1.205E-03,
    "blood": 1.060E+00,
    "bone": 1.920E+00,
    "calcium": 1.550E+00,
    "cha": 3.225E+00,
    "k2hpo4": 2.440E+00,
    "muscle": 1.050E+00,
    "water": 1.000E+00,
    "softtissue": 1.000E+00,
    "redmarrow": 1.030E+00,
    "yellowmarrow": 0.980E+00,
    "spongiosa": 1.180E+00,
    "iodine": 4.930,
  }
  return rho


def ct_number_to_mu(hu=0.0, mu_water=1.0):
  """Convert CT number in HU to linear attenuation coefficient."""
  return mu_water * hu / 1000.0 + mu_water


def mu_to_ct_number(mu=0.0, mu_water=1.0):
  """Convert linear attenuation coefficient to CT number in HU."""
  return (mu - mu_water) / mu_water * 1000.0


def mu_from_mass_attenuation_coefficient(mu_rho=1.0, mass_density=1.0):
  """Convert mass attenuation coefficient [cm^2/g] to linear attenuation [cm^-1]."""
  return mu_rho * mass_density


def phantom_linear_fit(x=[1, 2, 3], y=[6, 20, 32]):
  """Determine linear fit of density to Hounsfield units.

  x: densities (dependent), y: CT numbers (independent). Requires at least two points.
  """
  if len(x) != len(y):
    raise ValueError(
      "phantom_linear_fit: Number of dependent and independent variables must be equal."
    )

  if len(x) < 2:
    raise ValueError("phantom_linear_fit: Number of rods must be two or more.")

  deg = 1
  m, b = np.polyfit(x, y, deg)
  return [m, b]


def interpolate_mass_attenuation(material="water", keV=90, method="cubic", smoothing_factor=None):
  """Interpolate mass attenuation [cm2/g] from NIST tables at given energy (keV).

  method: 'linear', 'cubic', or 'log-log' (uses a univariate spline in log-log space).
  smoothing_factor: smoothing parameter 's' for UnivariateSpline when using 'log-log'.
                    If None, s=0 (interpolating spline).

  Default method is 'cubic'.
  """
  method_norm = method.lower().replace("_", "").replace("-", "")
  if method_norm not in {"linear", "cubic", "loglog"}:
    raise ValueError("interpolate_mass_attenuation: method must be 'linear', 'cubic' or 'log-log'.")

  table_name = material + "_table"
  if table_name not in mass_attenuation_tables.keys():
    raise ValueError(
      f"interpolate_mass_attenuation: Material {material} is not available in current NIST tables."
    )

  table = mass_attenuation_tables[table_name]

  # truncate table at material K-edge for selected materials (we need to be conservative to avoid interpolation errors)
  material_key = material.lower()
  k_edges = {
    "iodine": 34.5,  # I K-edge ~33.17 keV
    "calcium": 4.5,   # Ca K-edge ~4.04 keV, but we use a conservative estimate
    "cha": 4.5,       # calcium hydroxyapatite ~ Ca K-edge
    "bone": 4.5,      # bone dominated by Ca K-edge
    "air": 3.5,       # air ~3.2 keV
    "muscle": 4.0,   # muscle - we seem to have an issue with our NIST data not monotonically increasing at values below 4 keV
    "adipose": 3.0  # adipose - we seem to have an issue with our NIST data not monotonically increasing at values below 3 keV
  }

  if material_key in k_edges:
    k = float(k_edges[material_key])
    xp_min = float(table["Energy [keV]"].min())
    xp_max = float(table["Energy [keV]"].max())
    if k > xp_max:
      raise ValueError(
        f"interpolate_mass_attenuation: K-edge {k} keV for material {material} is outside table energy range [{xp_min}, {xp_max}]."
      )
    # only apply truncation if K-edge is above the current minimum
    # note that a problem with this code is that there are two mass attenuations at a given energy, and so we discard both. This affects, for example, iodine at 33 keV.
    if k > xp_min:
      table = table[table["Energy [keV]"] >= k]

  fp = table["Mass Attenuation [cm2/g]"].to_numpy()
  xp = table["Energy [keV]"].to_numpy()

  keV_arr = np.asarray(keV)
  if np.any(keV_arr < xp.min()) or np.any(keV_arr > xp.max()):
    raise ValueError(
      f"interpolate_mass_attenuation: Input keV {keV} is below or close to k-edge and cannot be interpreted from NIST table for material {material}."
    )

  if method_norm == "linear":
    result = np.interp(keV_arr, xp, fp)
    return result.item() if keV_arr.ndim == 0 else result

  if method_norm == "cubic":
    spl = CubicSpline(xp, fp)
    result = spl(keV_arr)
    return result.item() if np.ndim(keV_arr) == 0 else result

  # method_norm == "loglog" -> perform interpolation in log-log space using UnivariateSpline
  if np.any(xp <= 0) or np.any(fp <= 0) or np.any(keV_arr <= 0):
    raise ValueError("interpolate_mass_attenuation: log-log interpolation requires positive energies and values.")


  log_xp = np.log(xp)
  log_fp = np.log(fp)
  log_keV = np.log(keV_arr)

  s = 0 if smoothing_factor is None else float(smoothing_factor)
  spline = UnivariateSpline(log_xp, log_fp, s=s, k=3)
  log_result = spline(log_keV)
  result = np.exp(log_result)
  return result.item() if keV_arr.ndim == 0 else result


def solve_system_equations(A=np.ones((2, 2)), B=np.ones(2), use_numpy=False):
  """Solve AX = B. Supports numpy solver or Cramer's rule for 2x2 and 3x3 systems."""
  if use_numpy:
    a = np.array(A)
    b = np.array(B)
    X = np.linalg.solve(a, b)
    return X

  size = len(A)
  if size < 2 or size > 3:
    raise ValueError("Matrix size must be 2 or 3.")

  if size == 2:
    a = A[0, 0]
    b = A[0, 1]
    c = A[1, 0]
    d = A[1, 1]
    e = B[0]
    f = B[1]
    D = a * d - b * c
    Dx = d * e - f * b
    Dy = a * f - c * e

    x = Dx / D
    y = Dy / D
    return np.array([x, y])

  if size == 3:
    a = A[0, 0]
    b = A[0, 1]
    c = A[0, 2]
    d = A[1, 0]
    e = A[1, 1]
    f = A[1, 2]
    g = A[2, 0]
    h = A[2, 1]
    i = A[2, 2]
    j = B[0]
    k = B[1]
    l = B[2]
    D = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
    Dx = j * (e * i - f * h) - b * (k * i - f * l) + c * (k * h - e * l)
    Dy = a * (k * i - f * l) - j * (d * i - f * g) + c * (d * l - g * k)
    Dz = a * (e * l - h * k) - b * (d * l - g * k) + j * (d * h - e * g)

    x = Dx / D
    y = Dy / D
    z = Dz / D
    return np.array([x, y, z])


def parse_filename(filename):
  """Parse filename into basename, dirname, name (without extension), and extension.
  Handles single extensions (.jpg, .nii, .tar) and common double/extensions (.nii.gz, .tar.gz, .tar.bz2).
  """
  filename = str(filename)
  basename = os.path.basename(filename)
  dirname = os.path.dirname(filename)

  root, ext = os.path.splitext(basename)
  ext = ext.lower()

  # common compressed suffixes that indicate a double extension
  compressed_suffixes = {".gz", ".bz2", ".xz", ".zip"}
  if ext in compressed_suffixes:
    root2, ext2 = os.path.splitext(root)
    if ext2:
      ext = ext2.lower() + ext
      root = root2

  return basename, dirname, root, ext


def print_matrix(matrix, max_cols=6, float_fmt="{:12.4e}"):
  """Print matrix rows (truncate to max_cols). Handles ints and floats uniformly."""
  for row in matrix:
    row_to_print = row[:max_cols] if len(row) > max_cols else row
    for element in row_to_print:
      if isinstance(element, (int, np.integer)):
        print(f"{element:4d}", end=" ")
      else:
        try:
          print(float_fmt.format(float(element)), end=" ")
        except Exception:
          print(str(element), end=" ")
    print()
