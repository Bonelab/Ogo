'''Calibration of a standard phantom'''

import numpy as np
from scipy import stats
from .calibration import Calibration


class StandardCalibration(Calibration):
    '''Calibration for standard calibration phantoms

    This algorithm takes samples of a standard calibration phantom and
    determines the equation of best fit. This is typically done for
    solid phantoms. The calibration equation is straight forward:

    .. math::
        \\rho = m * HU + b

    Given sames of density and Hounsfield Units in the scan field of view,
    an equation of best fit is determined.

    The output densities are in units of whatever calibration phantom was used.
    It is the users responsibility to know what those units are.
    '''

    def __init__(self):
        super(StandardCalibration, self).__init__()
        self._slope = 0.0
        self._intercept = 0.0
        self._r_value = 0.0
        self._p_value = 0.0
        self._std_err = 0.0

    @property
    def slope(self):
        '''Return the calibration slope'''
        return self._slope

    @property
    def intercept(self):
        '''Return the calibration intercept'''
        return self._intercept

    @property
    def r_value(self):
        '''Return the R value (coefficient of correlation) for the fit'''
        return self._r_value

    @property
    def p_value(self):
        '''Return the coefficient of correlation p-value'''
        return self._p_value

    @property
    def std_err(self):
        '''Return the standard error around the line of best fit'''
        return self._std_err

    def _fit(self, hounsfield_units, densities):
        '''Standard least squares fit'''
        self._slope, self._intercept, self._r_value, \
            self._p_value, self._std_err = \
            stats.linregress(hounsfield_units, densities)

    def _predict(self, hu):
        '''Standard linear equation prediction'''
        return self._slope * np.array(hu) + self._intercept

    def _predict_inverse(self, density):
        '''Standard linear equation inverse prediction'''
        return (np.array(density) - self._intercept) / self._slope
