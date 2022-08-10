'''Abstract class for calibration'''
from abc import ABCMeta, abstractmethod
ABC = ABCMeta('ABC', (object,), {})  # compatible with Python 2 *and* 3


class Calibration(ABC):
    '''Abstract base class for calibration.

    Classes that overwrite this method should implement the methods:
        - :func:`~Calibration._fit`
        - :func:`~Calibration._predict`
        - :func:`~Calibration._predict_inverse`

    The logic around fitting and then predicting is handled by the base class.

    In general, we seek a mapping from Hounsfield Units to units of density.

    .. math::
        \\phi: HU \\rightarrow \\rho

    Typically, these densities are in a mass-equivalent (such as
    :math:`K_2HPO_4` or hydroxyapatite).
    '''

    def __init__(self):
        '''This class initializes the class to handle '''
        self._is_fit = False

    def fit(self, hounsfield_units, densities):
        '''Given an array of Hounsfield Units and an array of densities, \
        determine the best fit'''
        if len(densities) != len(hounsfield_units) or len(densities) == 0:
            raise RuntimeError('Must have the same number of densities and \
                Hounsfield Units for calibration')
        self._fit(hounsfield_units, densities)
        self._is_fit = True

    def predict(self, hu):
        '''Having fit the calibration, predict the density from the \
        Hounsfield Units in a VOI'''
        if self._is_fit:
            return self._predict(hu)
        else:
            raise RuntimeError('Must fit before predict can be ran')

    def predict_inverse(self, density):
        '''Having fit the calibration, predict the Hounsfield Units from \
        the density in a VOI'''
        if self._is_fit:
            return self._predict_inverse(density)
        else:
            raise RuntimeError('Must fit before predict_inverse can be ran')

    @abstractmethod
    def _fit(self, hounsfield_units, densities):
        '''Internal abstract method for fitting parametric model'''
        raise NotImplementedError('Calibration is an abstract class')

    @abstractmethod
    def _predict(self, hu):
        '''Internal abstract method for making predictions from a model'''
        raise NotImplementedError('Calibration is an abstract class')

    @abstractmethod
    def _predict_inverse(self, density):
        '''Internal abstract method for making inverse predictions from a \
        model'''
        raise NotImplementedError('Calibration is an abstract class')

    @abstractmethod
    def get_dict(self):
        '''Return a dictionary of calibration information'''
        raise NotImplementedError('Calibration is an abstract class')
