# /------------------------------------------------------------------------------+
# | 22-AUG-2022                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+

"""Calibration of a mindways phantom"""

import numpy as np
from scipy import stats
import copy
from .standard_calibration import StandardCalibration


class MindwaysCalibration(StandardCalibration):
    """Perform calibration on a Mindways phantom QCT PRO.

    There is some nuance to performing density calibration with a Mindways
    phantom. The reason for this is that the material (|K2HPO4|) is
    dissolved in water. This creates a slightly different calibration equation
    because the content of water must be controled for.

    The calibration steps are as follows. First, the densities of water and of
    |K2HPO4| must be known in the phantom rods. These can be taken from
    the certificate of calibration provided with your phantom and is different
    for each phantom.

    Then, the equation of best fit to the Hounsfield units of each rod:

    .. math::
        \\mu_{ROI} = \\rho_{water} + \\sigma_{ref} \cdot \\rho_{K_2HPO_4} +
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
        \\mu_{ROI} = \\sigma_{CT} \cdot \\rho_{K_2HPO_4} + \\beta_{CT}

    In the calibration framework presented in Ogo, we want to solve the
    equation for :math:`\\rho_{K_2HPO_4}` which gives the parameters in
    standard calibration of:

    .. math::
        m = \\frac{1}{\\sigma_{CT}}

    .. math::
        b = \\frac{- \\beta_{CT}}{\\sigma_{CT}}

    .. |K2HPO4| replace:: K\ :sub:`2`\ HPO\ :sub:`4`

    [1] QCT PRO User Guide, Mindways Software, Inc. v5.0, rev 20110801
    """  # noqa: W605

    def __init__(self):
        super(MindwaysCalibration, self).__init__()

        self._sigma_ref = 0.0
        self._beta_ref = 0.0
        self._sigma_ct = 0.0
        self._beta_ct = 0.0

        self._water = None

    @property
    def sigma_ref(self):
        """Get the computed :math:`\\sigma_{ref}`"""
        return self._sigma_ref

    @property
    def beta_ref(self):
        """Get the computed :math:`\\beta_{ref}`"""
        return self._beta_ref

    @property
    def sigma_ct(self):
        """Get the computed :math:`\\sigma_{CT}`"""
        return self._sigma_ct

    @property
    def beta_ct(self):
        """Get the computed :math:`\\beta_{CT}`"""
        return self._beta_ct

    def fit(self, hounsfield_units, densities, water):
        """Override Calibration fit method.

        Mindways calibration phantom requires a water density"""
        if (
            len(water) != len(densities)
            or len(water) != len(hounsfield_units)
            or len(water) == 0
        ):
            raise RuntimeError('Please provide the water and k2hpo4 values \
              given with your calibration certificate')
        self._fit(hounsfield_units, densities, water)
        self._is_fit = True

        # Save for printing
        self._hu = copy.deepcopy(hounsfield_units)
        self._rho = copy.deepcopy(densities)
        self._water = copy.deepcopy(water)

    def _fit(self, hounsfield_units, densities, water):
        """Non-standard Mindways fit"""

        # Fit equation 3
        lhs = np.array(hounsfield_units) - np.array(water)
        results = stats.linregress(densities, lhs)
        
        self._sigma_ref = results.slope
        self._beta_ref = results.intercept
        self._r_value = results.rvalue
        self._p_value = results.pvalue
        self._std_err = results.stderr
        self._intercept_std_err = results.intercept_stderr
        
        # Use equations 4 & 5
        self._sigma_ct = self._sigma_ref - 0.2174
        self._beta_ct = self._beta_ref + 999.6

        # Store slope and intercept
        self._slope = 1.0 / self._sigma_ct
        self._intercept = -1.0 * self._beta_ct / self._sigma_ct

    def get_dict(self):
        d = super(MindwaysCalibration, self).get_dict()

        d['Water'] = str(self._water)

        d['Sigma ref'] = self.sigma_ref
        d['Beta ref'] = self.beta_ref
        d['Sigma ct'] = self.sigma_ct
        d['Beta ct'] = self.beta_ct

        return d
