'''Calibration of a mindways phantom'''

import numpy as np
from scipy import stats
from .standard_calibration import StandardCalibration


class MindwaysCalibration(StandardCalibration):
    '''Perform calibration on a Mindways phantom QCT PRO.

    There is some nuance to performing density calibration with a Mindways
    phantom. The reason for this is that the material (:math:`K_2HPO_4`) is
    dissolved in water. This creates a slightly different calibration equation
    because the content of water must be controled for.

    The calibration steps are as follows. First, the densities of water and of
    :math:`K_2HPO_4` must be known in the phantom rods. These can be taken from
    the certificate of calibration provided with your phantom and is different
    for each phantom.

    Then, the equation of best fit to the Hounsfield units of each rod:

    .. math::
        \\mu_{ROI} = \\rho_{water} + \\sigma_{ref} * \\rho_{K_2HPO_4} +
            \\beta_{ref}

    The coefficient of correlation return is of this equation.

    Next, there is a conversion from the water-dissolved equation to
    traditional density measures.

    .. math::
        \\sigma_{CT} = \\sigma_{ref} - 0.2174

    .. math::
        \\beta_{CT} = \\beta_{ref} + 999.6

    Finally, Hounsfield units and K2HPO4 equivalent density are related by the
    following equation:

    .. math::
        \\mu_{ROI} = \\sigma_{CT} * \\rho_{K_2HPO_4} + \\beta_{CT}

    In the calibration framework presented in Ogo, we want to solve the
    equation for :math:`\\rho_{K_2HPO_4}` which gives the parameters in
    standard calibration of:

    .. math::
        m = \\frac{1}{\\sigma_{CT}}

    .. math::
        b = \\frac{- \\beta_{CT}}{\\sigma_{CT}}

    [1] QCT PRO User Guide, Mindways Software, Inc. v5.0, rev 20110801
    '''

    def __init__(self):
        super(MindwaysCalibration, self).__init__()

        self._sigma_ref = 0.0
        self._beta_ref = 0.0
        self._sigma_ct = 0.0
        self._beta_ct = 0.0

    @property
    def sigma_ref(self):
        '''Get the computed :math:`\\sigma_{ref}`'''
        return self._sigma_ref

    @property
    def beta_ref(self):
        '''Get the computed :math:`\\beta_{ref}`'''
        return self._beta_ref

    @property
    def sigma_ct(self):
        '''Get the computed :math:`\\sigma_{CT}`'''
        return self._sigma_ct

    @property
    def beta_ct(self):
        '''Get the computed :math:`\\beta_{CT}`'''
        return self._beta_ct

    def fit(self, hounsfield_units, densities, water):
        '''Override Calibration fit method.

        Mindways calibration phantom requires a water density'''
        if len(water) != len(densities) \
                or len(water) != len(hounsfield_units) \
                or len(water) == 0:
            raise RuntimeError('Please provide the water and k2hpo4 values \
              given with your calibration certificate')
        self._fit(hounsfield_units, densities, water)
        self._is_fit = True

    def _fit(self, hounsfield_units, densities, water):
        '''Non-standard Mindways fit'''

        # Fit equation 3
        lhs = np.array(hounsfield_units) - np.array(water)
        self._sigma_ref, self._beta_ref, self._r_value, \
            self._p_value , self._std_err = \
            stats.linregress(densities, lhs)

        # Use equations 4 & 5
        self._sigma_ct = self._sigma_ref - 0.2174
        self._beta_ct = self._beta_ref + 999.6

        # Store slope and intercept
        self._slope = 1.0 / self._sigma_ct
        self._intercept = -1.0 * self._beta_ct / self._sigma_ct
