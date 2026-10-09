"""Scalar bone laws: density in mg/cm3, modulus and yield stress in MPa.

Density convention is set by the anatomy builder (ash for hip, calibrated QCT
for spine); these functions only convert mg/cm3 to g/cm3.
"""

from typing import Callable

# === Morgan et al. (2003) ===
def morgan_trab_E(density_mgcc: float) -> float:
    """Return Morgan elastic modulus (MPa) from density (mg/cm3)."""
    if density_mgcc < 0:
        return 0.0001
    rho = density_mgcc / 1000.0  # Convert to g/cc
    return 4730 * rho ** 1.56  # MPa

# === Bayraktar et al. (2004) tissue-level constants ===
def bayraktar_trab_E(density_mgcc: float) -> float:
    """Return the constant trabecular tissue modulus (MPa); density is unused."""
    return 18000  # MPa

def bayraktar_trab_yc(density_mgcc: float) -> float:
    """Return the constant trabecular compressive yield stress (MPa)."""
    return 134  # MPa (compressive)

def bayraktar_trab_yt(density_mgcc: float) -> float:
    """Return the constant trabecular tensile yield stress (MPa)."""
    return 83  # MPa (tensile)

# === Kopperdahl et al. (2002) ===
def kopperdahl_trab_E(density_mgcc: float) -> float:
    """Return spine elastic modulus (MPa) from calibrated density (mg/cm3)."""
    if density_mgcc < 0:
        return 0.0001
    rho = density_mgcc / 1000.0
    return 2980 * rho ** 1.05  # MPa

def kopperdahl_trab_yc(density_mgcc: float) -> float:
    """Return compressive yield stress (MPa) from density (mg/cm3)."""
    if density_mgcc < 0:
        return 0.0001
    rho = density_mgcc / 1000.0
    return 37.4 * rho ** 1.39  # MPa

# === Crawford et al. (2003) ===
def crawford_voxel_E(density_mgcc: float) -> float:
    """Return nonnegative Crawford voxel modulus (MPa) from density (mg/cm3)."""
    if density_mgcc < 0:
        return 0.0001
    rho = density_mgcc / 1000.0
    return max(0.0, -34.7 + 3230 * rho)  # MPa


# === Cortical estimates using Bayraktar ===
def bayraktar_cort_E(density_mgcc: float) -> float:
    """Return the constant cortical tissue modulus (MPa); density is unused."""
    return 20000  # MPa

def bayraktar_cort_yc(density_mgcc: float) -> float:
    """Return the constant cortical compressive yield stress (MPa)."""
    return 134  # MPa

def bayraktar_cort_yt(density_mgcc: float) -> float:
    """Return the constant cortical tensile yield stress (MPa)."""
    return 83  # MPa

# === Default  functions (adapted from trabecular with different scaling) ===
def default_E(density_mgcc: float) -> float:
    """Return hip elastic modulus (MPa) from ash density (mg/cm3)."""
    if density_mgcc < 0:
        return 0.0001
    rho = density_mgcc / 1000.0
    return 10500 * rho ** 2.29  # MPa

def default_yc(density_mgcc: float) -> float:
    """Return the default compressive yield stress (MPa)."""
    return 120 

def default_yt(density_mgcc: float) -> float:
    """Return the default tensile yield stress (MPa)."""
    return 100 

# === Define custom functions via command line === #
def constant(value: float) -> Callable[[float], float]:
    """Return a density-independent scalar law."""
    # use as: --elastic_E_func constant_Value
    return lambda density_mgcc: value

def linear(A: float, B: float) -> Callable[[float], float]:
    """Return A + B * density, with density converted to g/cm3."""
    # use as: --elastic_E_func linear_A_b

    return lambda density_mgcc: A + B * (density_mgcc / 1000.0)

def power_law(A: float, b: float) -> Callable[[float], float]:
    """Return A * density**b, with density converted to g/cm3."""
      # use as: --elastic_E_func power_A_b
    return lambda density_mgcc: A * (density_mgcc / 1000.0) ** b
