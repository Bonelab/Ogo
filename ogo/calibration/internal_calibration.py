'''Internal density calibration'''

from .standard_calibration import StandardCalibration
from .mass_attenuation_tables import mass_attenuation_tables
import pandas as pd
import numpy as np
import scipy.interpolate as interp
from scipy import stats
from collections import OrderedDict
import copy


class InternalCalibration(StandardCalibration):
    '''Perform internal density calibration.

    For details on the method, please see [1]. Mass attenuation coefficients
    as a function of energy come from [2, 3]. Each sample should come from a
    minimum area of 10 mm\ :sup:`2`.

    For material definitions used, see
    :class:`ogo.calibration.mass_attenuation_tables`. Calibrated images are in
    units of mg K\ :sub:`2`\ HPO\ :sub:`4` /cc.

    [1] Michalski, Andrew S., et al. "CT-based internal density calibration for opportunistic skeletal assessment using abdominal CT scans." Medical Engineering & Physics (2020).

    [2] https://www.nist.gov/pml/x-ray-mass-attenuation-coefficients

    [3] White DR, Wilson IJ, Griffith RV. Report 46. J Int Comm Radiat Units Meas 2016 os24:NP-NP. doi: 10.1093/jicru/os24.1.Report46
    '''  # noqa: W605, E501

    def __init__(self,
                 adipose_hu=0, air_hu=0, blood_hu=0, bone_hu=0, muscle_hu=0):
        super(InternalCalibration, self).__init__()

        # User input values
        self._adipose_hu = adipose_hu
        self._air_hu = air_hu
        self._blood_hu = blood_hu
        self._bone_hu = bone_hu
        self._muscle_hu = muscle_hu

        # Computed values
        self._effective_energy = 0
        self._max_r2 = 0

        self._adipose_mass_attenuation = 0
        self._air_mass_attenuation = 0
        self._blood_mass_attenuation = 0
        self._bone_mass_attenuation = 0
        self._muscle_mass_attenuation = 0

        self._K2HPO4_mass_attenuation = 0
        self._CHA_mass_attenuation = 0
        self._triglyceride_mass_attenuation = 0
        self._water_mass_attenuation = 0

        self._hu_to_mass_attenuation_slope = 0
        self._hu_to_mass_attenuation_intercept = 0

        self._hu_to_density_slope = 0
        self._hu_to_density_intercept = 0

        self._voxel_volume = 0

    def predict(self, hu, voxel_volume):
        '''Internal calibration predict method

        It is recommended to cast hu to type float (sitk.sitkFloat64) before
        passing to this function. Operators (*, +, /) are used without
        reference to type to allow this function to work with SimpleITK, numpy,
        and base python. Therefore, no explicite checks for type are performed
        internally.

        :param hu: Input Hounsfield unit
        :param voxel_volume: Voxel volume in  mm\ :sup:`3`
        '''  # noqa: W605
        if self._is_fit:
            return self._predict(hu, voxel_volume)
        else:
            raise RuntimeError('Must fit before predict can be ran')

    def _predict(self, hu, voxel_volume):
        '''Internal calibration _predict method'''

        # Convert voxel volume to cm3
        self._voxel_volume = copy.deepcopy(voxel_volume)
        voxel_volume = voxel_volume / 1000.0

        # Convert HU to mass attenuation coefficient
        u_p = hu * self.hu_to_mass_attenuation_slope + \
            self.hu_to_mass_attenuation_intercept

        # Convert HU to Archimedian density
        arch = hu * self.hu_to_density_slope + \
            self.hu_to_density_intercept

        # Conver Archimedian density to total mass
        mass = arch * voxel_volume

        # Create two component model
        k2hpo4 = mass * (
            (u_p - self.triglyceride_mass_attenuation) /
            (self.K2HPO4_mass_attenuation - self.triglyceride_mass_attenuation)
        )

        # Convert g to mg
        k2hpo4 = k2hpo4 / voxel_volume * 1000.0

        return k2hpo4

    def fit(self):
        '''Override Calibration fit method.

        Internal calibration requires specific sampled values'''
        self._fit()
        self._is_fit = True

    def _fit(self):
        '''Internal calibration fit

        Performs a grid search on all energies to find an energy which
        maximizes the coefficient of determination between sampled
        Hounsfield Units and mass attenuation. The best correlated energy
        is used as the scan effective energy.'''

        self._interpolate_mass_attenuation()
        self._determine_scan_effective_energy()
        self._determine_hu_to_mass_attenuation()
        self._determine_hu_to_density()

    def _interpolate_mass_attenuation(self):
        '''Interpolate the mass attenuation curves'''
        # Constants
        energies = np.arange(1, 200.5, 0.5)

        # Interpolate each table
        self._interpolate_tables = {}
        for table_name in mass_attenuation_tables.keys():
            # Get the table
            this_table = mass_attenuation_tables[table_name]

            # Interpolate to specified energies
            interp_table = interp.griddata(
                this_table['Energy [keV]'],
                this_table['Mass Attenuation [cm2/g]'],
                energies,
                method='linear'
            )

            # Store as the same pandas DataFrame
            self._interpolate_tables[table_name] = pd.DataFrame({
                'Energy [keV]':             energies,
                'Mass Attenuation [cm2/g]': interp_table
            })

        self._interpolate_tables

    def _determine_scan_effective_energy(self):
        '''Determine scan effective energy'''
        def _get_mass_attenuation_at_index(idx):
            '''Return the mass attenuation array at a given index'''
            return [
                self._interpolate_tables['adipose_table'].loc[idx, 'Mass Attenuation [cm2/g]'],  # noqa: E501
                self._interpolate_tables['air_table'].loc[idx, 'Mass Attenuation [cm2/g]'],      # noqa: E501
                self._interpolate_tables['blood_table'].loc[idx, 'Mass Attenuation [cm2/g]'],    # noqa: E501
                self._interpolate_tables['bone_table'].loc[idx, 'Mass Attenuation [cm2/g]'],     # noqa: E501
                self._interpolate_tables['muscle_table'].loc[idx, 'Mass Attenuation [cm2/g]'],   # noqa: E501
            ]

        # Measured HU values
        HU = [
            self.adipose_hu,
            self.air_hu,
            self.blood_hu,
            self.bone_hu,
            self.muscle_hu
        ]

        n = len(self._interpolate_tables['adipose_table'])
        max_index = -1
        max_r2 = -np.Inf

        for i in np.arange(1, n, 1):
            # Get the values at this energy level
            attenuation = _get_mass_attenuation_at_index(i)

            # Lest squares fit
            EE_lr = stats.linregress(HU, attenuation)
            r_squared = EE_lr[2]**2

            # Take best fit
            if r_squared > max_r2:
                max_r2 = r_squared
                max_index = i

        # Set values
        self._effective_energy = self._interpolate_tables['adipose_table'].loc[max_index, 'Energy [keV]']  # noqa: E501
        self._max_r2 = max_r2

        self._adipose_mass_attenuation = self._interpolate_tables['adipose_table'].loc[max_index, 'Mass Attenuation [cm2/g]']  # noqa: E501
        self._air_mass_attenuation = self._interpolate_tables['air_table'].loc[max_index, 'Mass Attenuation [cm2/g]']          # noqa: E501
        self._blood_mass_attenuation = self._interpolate_tables['blood_table'].loc[max_index, 'Mass Attenuation [cm2/g]']      # noqa: E501
        self._bone_mass_attenuation = self._interpolate_tables['bone_table'].loc[max_index, 'Mass Attenuation [cm2/g]']        # noqa: E501
        self._muscle_mass_attenuation = self._interpolate_tables['muscle_table'].loc[max_index, 'Mass Attenuation [cm2/g]']    # noqa: E501

        self._K2HPO4_mass_attenuation = self._interpolate_tables['k2hpo4_table'].loc[max_index, 'Mass Attenuation [cm2/g]']              # noqa: E501
        self._CHA_mass_attenuation = self._interpolate_tables['cha_table'].loc[max_index, 'Mass Attenuation [cm2/g]']                    # noqa: E501
        self._triglyceride_mass_attenuation = self._interpolate_tables['triglyceride_table'].loc[max_index, 'Mass Attenuation [cm2/g]']  # noqa: E501
        self._water_mass_attenuation = self._interpolate_tables['water_table'].loc[max_index, 'Mass Attenuation [cm2/g]']                # noqa: E501

        # Save memory
        del self._interpolate_tables

    def _determine_hu_to_mass_attenuation(self):
        '''Compute HU to mass attenuation relationship'''
        # Get values
        HU = [
            self.adipose_hu,
            self.air_hu,
            self.blood_hu,
            self.bone_hu,
            self.muscle_hu
        ]
        attenuation = [
            self.adipose_mass_attenuation,
            self.air_mass_attenuation,
            self.blood_mass_attenuation,
            self.bone_mass_attenuation,
            self.muscle_mass_attenuation
        ]

        # Perform linear regression
        linreg = stats.linregress(HU, attenuation)

        # Store values
        self._hu_to_mass_attenuation_slope = linreg[0]
        self._hu_to_mass_attenuation_intercept = linreg[1]

    def _determine_hu_to_density(self):
        '''Compute HU to density relationship'''
        # Get HU
        HU = [
            self.adipose_hu,
            self.air_hu,
            self.blood_hu,
            self.bone_hu,
            self.muscle_hu
        ]

        # Convert to densities
        def _convert_value(hu, attenuation):
            water_attenuation = self.water_mass_attenuation
            water_density = 1.0

            return (hu/1000*water_attenuation*water_density + water_attenuation*water_density)/attenuation  # noqa: E501

        densities = [
            _convert_value(self.adipose_hu, self.adipose_mass_attenuation),
            _convert_value(self.air_hu, self.air_mass_attenuation),
            _convert_value(self.blood_hu, self.blood_mass_attenuation),
            _convert_value(self.bone_hu, self.bone_mass_attenuation),
            _convert_value(self.muscle_hu, self.muscle_mass_attenuation)
        ]

        # Perform linear regression
        linreg = stats.linregress(HU, densities)

        # Store values
        self._hu_to_density_slope = linreg[0]
        self._hu_to_density_intercept = linreg[1]

    @property
    def adipose_hu(self):
        '''Sampled adipose hounsfield units'''
        return self._adipose_hu

    @adipose_hu.setter
    def adipose_hu(self, HU):
        '''Sampled adipose hounsfield units'''
        self._adipose_hu = HU

    @property
    def air_hu(self):
        '''Sampled air hounsfield units'''
        return self._air_hu

    @air_hu.setter
    def air_hu(self, HU):
        '''Sampled air hounsfield units'''
        self._air_hu = HU

    @property
    def blood_hu(self):
        '''Sampled blood hounsfield units'''
        return self._blood_hu

    @blood_hu.setter
    def blood_hu(self, HU):
        '''Sampled blood hounsfield units'''
        self._blood_hu = HU

    @property
    def bone_hu(self):
        '''Sampled bone hounsfield units'''
        return self._bone_hu

    @bone_hu.setter
    def bone_hu(self, HU):
        '''Sampled bone hounsfield units'''
        self._bone_hu = HU

    @property
    def muscle_hu(self):
        '''Sampled muscle hounsfield units'''
        return self._muscle_hu

    @muscle_hu.setter
    def muscle_hu(self, HU):
        '''Sampled muscle hounsfield units'''
        self._muscle_hu = HU

    @property
    def effective_energy(self):
        '''Calculated effective energy [keV]'''
        return self._effective_energy

    @property
    def max_r2(self):
        '''Maximum coefficient of determination'''
        return self._max_r2

    @property
    def adipose_mass_attenuation(self):
        '''Adipose mass attenuation coefficient at the effective energy'''
        return self._adipose_mass_attenuation

    @property
    def air_mass_attenuation(self):
        '''Air mass attenuation coefficient at the effective energy'''
        return self._air_mass_attenuation

    @property
    def blood_mass_attenuation(self):
        '''Blood mass attenuation coefficient at the effective energy'''
        return self._blood_mass_attenuation

    @property
    def bone_mass_attenuation(self):
        '''Bone mass attenuation coefficient at the effective energy'''
        return self._bone_mass_attenuation

    @property
    def muscle_mass_attenuation(self):
        '''Muscle mass attenuation coefficient at the effective energy'''
        return self._muscle_mass_attenuation

    @property
    def K2HPO4_mass_attenuation(self):
        '''K2HPO4 mass attenuation coefficient at the effective energy'''
        return self._K2HPO4_mass_attenuation

    @property
    def CHA_mass_attenuation(self):
        '''CHA mass attenuation coefficient at the effective energy'''
        return self._CHA_mass_attenuation

    @property
    def triglyceride_mass_attenuation(self):
        '''Triglyceride mass attenuation coefficient at the effective energy'''
        return self._triglyceride_mass_attenuation

    @property
    def water_mass_attenuation(self):
        '''Water mass attenuation coefficient at the effective energy'''
        return self._water_mass_attenuation

    @property
    def hu_to_mass_attenuation_slope(self):
        '''HU to mass attenuation slope'''
        return self._hu_to_mass_attenuation_slope

    @property
    def hu_to_mass_attenuation_intercept(self):
        '''HU to mass attenuation intercept'''
        return self._hu_to_mass_attenuation_intercept

    @property
    def hu_to_density_slope(self):
        '''HU to density slope'''
        return self._hu_to_density_slope

    @property
    def hu_to_density_intercept(self):
        '''HU to density intercept'''
        return self._hu_to_density_intercept

    def get_dict(self):
        return OrderedDict([
            ('Is fit',       self._is_fit),

            ('Adipose [HU]', self.adipose_hu),
            ('Air [HU]',     self.air_hu),
            ('Blood [HU]',   self.blood_hu),
            ('Bone [HU]',    self.bone_hu),
            ('Muscle [HU]',  self.muscle_hu),

            ('Voxel Volume [mm3]',   self._voxel_volume),

            ('Effective Energy [keV]',   self.effective_energy),
            ('Max R^2',                  self.max_r2),

            ('Adipose u/p [cm2/g]',      self.adipose_mass_attenuation),
            ('Air u/p [cm2/g]',          self.air_mass_attenuation),
            ('Blood u/p [cm2/g]',        self.blood_mass_attenuation),
            ('Bone u/p [cm2/g]',         self.bone_mass_attenuation),
            ('Muscle u/p [cm2/g]',       self.muscle_mass_attenuation),
            ('K2HPO4 u/p [cm2/g]',       self.K2HPO4_mass_attenuation),
            ('CHA u/p [cm2/g]',          self.CHA_mass_attenuation),
            ('Triglyceride u/p [cm2/g]', self.triglyceride_mass_attenuation),
            ('Water u/p [cm2/g]',        self.water_mass_attenuation),

            ('HU-u/p Slope',        self.hu_to_mass_attenuation_slope),
            ('HU-u/p Y-Intercept',  self.hu_to_mass_attenuation_intercept),

            ('HU-Material Density Slope',         self.hu_to_density_slope),
            ('HU-Material Density Y-Intercept',   self.hu_to_density_intercept)
        ])
