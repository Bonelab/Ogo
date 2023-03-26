# /------------------------------------------------------------------------------+
# | 22-AUG-2022                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+

"""Calibration of a standard phantom"""

from scipy import stats
import copy
from collections import OrderedDict
from .calibration import Calibration
import numpy as np
import SimpleITK as sitk


class StandardCalibration(Calibration):
    """Calibration for standard calibration phantoms

    This algorithm takes samples of a standard calibration phantom and
    determines the equation of best fit. This is typically done for
    solid phantoms. The calibration equation is straight forward:

    .. math::
        \\rho = m \cdot HU + b

    Given sames of density and Hounsfield Units in the scan field of view,
    an equation of best fit is determined.

    The output densities are in units of whatever calibration phantom was used.
    It is the user's responsibility to know what those units are.
    """  # noqa: W605

    def __init__(self):
        super(StandardCalibration, self).__init__()
        self._slope = 0.0
        self._intercept = 0.0
        self._r_value = 0.0
        self._p_value = 0.0
        self._std_err = 0.0
        self._intercept_std_err = 0.0

        self._hu = None
        self._rho = None

    @property
    def slope(self):
        """Return the calibration slope"""
        return self._slope

    @property
    def intercept(self):
        """Return the calibration intercept"""
        return self._intercept

    @property
    def r_value(self):
        """Return the R value (coefficient of correlation) for the fit"""
        return self._r_value

    @property
    def p_value(self):
        """Return the coefficient of correlation p-value"""
        return self._p_value

    @property
    def std_err(self):
        """Return the standard error around the line of best fit"""
        return self._std_err
    
    @property
    def std_err_intercept(self):
        """Return the standard error for the intercept of the line of best fit"""
        return self._intercept_std_err

    def _fit(self, hounsfield_units, densities):
        """Standard least squares fit"""
 
        results = stats.linregress(hounsfield_units, densities)
        self._slope = results.slope
        self._intercept = results.intercept
        self._r_value = results.rvalue
        self._p_value = results.pvalue
        self._std_err = results.stderr
        self._intercept_std_err = results.intercept_stderr


        # Save for printing
        self._hu = copy.deepcopy(hounsfield_units)
        self._rho = copy.deepcopy(densities)

    def _predict(self, hu):
        """Standard linear equation prediction"""
        return hu * self._slope + self._intercept
    
    def montecarlo_predict(self, hu):
        """Standard linear equation prediction"""
        density_uncertainty = sitk.Sqrt(
            ((hu * self._std_err) ** 2) + (self._intercept_std_err ** 2)
        )
        return density_uncertainty

    def _predict_inverse(self, density):
        """Standard linear equation inverse prediction"""
        return (density - self._intercept) / self._slope

    def get_dict(self):
        return OrderedDict([
            ('Is fit',          self._is_fit),

            ('Slope',           self.slope),
            ('Intercept',       self.intercept),
            ('R value',         self.r_value),
            ('P value',         self.p_value),
            ('Standard Error',  self.std_err),
            ('Intercept Std Err', self.std_err_intercept)

            ('HU',          str(self._hu)),
            ('Densities',   str(self._rho))
        ])
