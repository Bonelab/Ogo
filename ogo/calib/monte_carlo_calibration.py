# /------------------------------------------------------------------------------+
# | 22-AUG-2022                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+

"""Monte Carlo Internal density calibration"""

from .standard_calibration import StandardCalibration
from ogo.dat.MassAttenuationTables import mass_attenuation_tables
import pandas as pd
import numpy as np
import scipy.interpolate as interp
from scipy import stats
from scipy.optimize import minimize_scalar
from collections import OrderedDict
import SimpleITK as sitk 
import math 
import copy
from tqdm import tqdm



class MonteCarloCalibration(StandardCalibration):
    """Perform internal density calibration.
    For details on the method, please see [1]. Mass attenuation coefficients
    as a function of energy come from [2, 3]. Each sample should come from a
    minimum area of 10 mm\ :sup:`2`.
    For material definitions used, see
    :class:`ogo.calibration.mass_attenuation_tables`. Calibrated images are in
    units of mg K\ :sub:`2`\ HPO\ :sub:`4` /cc.
    [1] Michalski, Andrew S., et al. "CT-based internal density calibration for opportunistic skeletal assessment using abdominal CT scans." Medical Engineering & Physics (2020).
    [2] https://www.nist.gov/pml/x-ray-mass-attenuation-coefficientslabel_list
    [3] White DR, Wilson IJ, Griffith RV. Report 46. J Int Comm Radiat Units Meas 2016 os24:NP-NP. doi: 10.1093/jicru/os24.1.Report46
    """  # noqa: W605, E501

    def __init__(self,
                 iterations=0, adipose_hu=0, air_hu=0, blood_hu=0, bone_hu=0, muscle_hu=0, adipose_std=0, 
                 air_std=0, blood_std=0, bone_std=0, muscle_std=0, label_list=(91, 92, 93, 94, 95)):
        super(MonteCarloCalibration, self).__init__()

        # User input values
        self.iterations = iterations

        self._adipose_hu = adipose_hu
        self._air_hu = air_hu
        self._blood_hu = blood_hu
        self._bone_hu = bone_hu
        self._muscle_hu = muscle_hu

        self.adipose_std = adipose_std
        self.air_std = air_std
        self.blood_std = blood_std
        self.bone_std = bone_std
        self.muscle_std = muscle_std

        self._label_list = list(label_list)

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
        self._hu_to_mass_attenuation_slope_error = 0
        self._hu_to_mass_attenuation_intercept_error = 0

        self._hu_to_density_slope = 0
        self._hu_to_density_intercept = 0
        self._hu_to_density_slope_error = 0
        self._hu_to_density_intercept_error = 0

        self.k2hpo4_stdev = 0
        self.triglyceride_stddev = 0
        self.k2hpo4_standard_error = 0
        self.triglyceride_standard_error = 0

        self._voxel_volume = 0
        
    def _montecarlo_predict(self, hu, voxel_volume):
            """Internal calibration predict method
            It is recommended to cast hu to type float (sitk.sitkFloat64) before
            passing to this function. Operators (*, +, /) are used without
            reference to type to allow this function to work with SimpleITK, numpy,
            and base python. Therefore, no explicite checks for type are performed
            internally.
            :param hu: Input Hounsfield unit
            :param voxel_volume: Voxel volume in  mm\ :sup:`3`
            """  # noqa: W605
            if self._is_fit:
                return self.montecarlo_predict(hu, voxel_volume)
            else:
                raise RuntimeError('Must fit before predict can be run.')

    
    def montecarlo_predict(self, hu, voxel_volume):
        """Internal calibration _predict method"""

        # Convert voxel volume to cm3
        self._voxel_volume = copy.deepcopy(voxel_volume)
        voxel_volume = voxel_volume / 1000.0

        # Convert HU to mass attenuation coefficient
        u_p = (
          hu * self.hu_to_mass_attenuation_slope
          + self.hu_to_mass_attenuation_intercept
        )

        # Convert HU to Archimedian density
        arch = (
            hu * self.hu_to_density_slope
            + self.hu_to_density_intercept
        )

        # Conver Archimedian density to total mass
        mass = arch * voxel_volume

        # Create two component model
        k2hpo4 = mass * (
                (u_p - self.triglyceride_mass_attenuation) /
                (self.K2HPO4_mass_attenuation - self.triglyceride_mass_attenuation)
        )

        # Convert g to mg
        k2hpo4 = k2hpo4 / voxel_volume * 1000.0
        
        error_hu_mass_attenuation = sitk.Sqrt(
            sitk.Square((hu * self._hu_to_mass_attenuation_slope_error))
            + (self._hu_to_mass_attenuation_intercept_error) ** 2
        )

        error_hu_density = sitk.Sqrt(
            sitk.Square((hu * self._hu_to_density_slope_error))
            + (self._hu_to_density_intercept_error) ** 2
        )

        mu_k_minus_mu_t = self.K2HPO4_mass_attenuation - self.triglyceride_mass_attenuation
        sq_mu_k_minus_mu_t = mu_k_minus_mu_t ** 2
        error_mass = sitk.Multiply( error_hu_density, voxel_volume)

        #error propagation 
        uncertainty_k2hpo4 = sitk.Sqrt(
                sitk.Square(sitk.Multiply(sitk.Divide(u_p - self.triglyceride_mass_attenuation, mu_k_minus_mu_t), error_mass))
                + sitk.Square(sitk.Multiply(sitk.Divide(mass, mu_k_minus_mu_t), error_hu_mass_attenuation))
                + 
                sitk.Square(
                    sitk.Multiply(
                        sitk.Divide(
                                sitk.Multiply(sitk.Subtract(self.triglyceride_mass_attenuation, u_p), self.K2HPO4_mass_attenuation),
                                sq_mu_k_minus_mu_t
                            ),
                        self.k2hpo4_standard_error
                    )
                )
                + 
                sitk.Square(
                    sitk.Multiply(
                        sitk.Divide(
                            sitk.Multiply(mass, sitk.Subtract(u_p, self.K2HPO4_mass_attenuation)),
                            sq_mu_k_minus_mu_t
                        ),
                        self.triglyceride_standard_error
                    )
                )
                )
        
        uncertainty_k2hpo4 = uncertainty_k2hpo4 / voxel_volume * 1000.0


        return uncertainty_k2hpo4

    def montecarlofit(self):
        """Override Calibration fit method.
        Internal calibration requires specific sampled values"""
        self.montecarlo_fit()
        self._is_fit = True

    def montecarlo_fit(self):
        """Internal calibration fit
        Performs a grid search on all energies to find an energy which
        maximizes the coefficient of determination between sampled
        Hounsfield Units and mass attenuation. The best correlated energy
        is used as the scan effective energy."""

        self._interpolate_mass_attenuation()
        self.montecarlo_determine_scan_effective_energy()
        self.montecarlo_determine_hu_to_mass_attenuation()
        self.montecarlo_determine_hu_to_density()

    def _subset(self, all_samples):
        """Returns only the samples will be employed in linear regressions.
        """

        sample_subset = []
        if 91 in self._label_list:
            sample_subset.append(all_samples[0])  # adipose
        if 92 in self._label_list:
            sample_subset.append(all_samples[1])  # air
        if 93 in self._label_list:
            sample_subset.append(all_samples[2])  # blood
        if 94 in self._label_list:
            sample_subset.append(all_samples[3])  # bone
        if 95 in self._label_list:
            sample_subset.append(all_samples[4])  # muscle

        return sample_subset

    def _interpolate_mass_attenuation(self):
        """Interpolate the mass attenuation curves"""
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
                'Energy [keV]': energies,
                'Mass Attenuation [cm2/g]': interp_table
            })

    def montecarlo_determine_scan_effective_energy(self):
        """Determine scan effective energy"""

        def _get_mass_attenuation_at_energy(energy):
            """Return the mass attenuation array at a given index"""
            return [self._subset([
                self._interpolate_tables['adipose_table'].loc[energy, 'Mass Attenuation [cm2/g]'],  # noqa: E501
                self._interpolate_tables['air_table'].loc[energy, 'Mass Attenuation [cm2/g]'],  # noqa: E501
                self._interpolate_tables['blood_table'].loc[energy, 'Mass Attenuation [cm2/g]'],  # noqa: E501
                self._interpolate_tables['bone_table'].loc[energy, 'Mass Attenuation [cm2/g]'],  # noqa: E501
                self._interpolate_tables['muscle_table'].loc[energy, 'Mass Attenuation [cm2/g]']  # noqa: E501
            ])
            ]
        
        effective_energy_arr = []
        k2hpo4 = []
        triglyceride = []

        n = len(self._interpolate_tables['adipose_table'])
        max_index = -1
        max_r2 = -np.Inf
    
        for k in tqdm(range(self.iterations)):
            air_hu_rand = np.random.normal(self.air_hu, self.air_std)
            adipose_hu_rand = np.random.normal(self.adipose_hu, self.adipose_std)
            blood_hu_rand = np.random.normal(self.blood_hu, self.blood_std)
            bone_hu_rand = np.random.normal(self.bone_hu, self.bone_std)
            muscle_hu_rand = np.random.normal(self.muscle_hu, self.muscle_std)
        
        
            # Measured HU values
            HU = self._subset([
                adipose_hu_rand,
                air_hu_rand,
                blood_hu_rand,
                bone_hu_rand,
                muscle_hu_rand
            ])

            for i in np.arange(100, n, 1):
                # Get the values at this energy level
                attenuation = _get_mass_attenuation_at_energy(i)

                # Least squares fit
                EE_lr = stats.linregress(HU, attenuation)
                r_squared = EE_lr[2] ** 2

                # Take best fit
                if r_squared > max_r2:
                    max_r2 = r_squared
                    max_index = i

            # Set values
            self._effective_energy = self._interpolate_tables['adipose_table'].loc[max_index, 'Energy [keV]']  # noqa: E501
            self._max_r2 = max_r2
            effective_energy_arr.append(self._effective_energy)
            

            self._adipose_mass_attenuation = self._interpolate_tables['adipose_table'].loc[
                max_index, 'Mass Attenuation [cm2/g]']  # noqa: E501
            self._air_mass_attenuation = self._interpolate_tables['air_table'].loc[
                max_index, 'Mass Attenuation [cm2/g]']  # noqa: E501
            self._blood_mass_attenuation = self._interpolate_tables['blood_table'].loc[
                max_index, 'Mass Attenuation [cm2/g]']  # noqa: E501
            self._bone_mass_attenuation = self._interpolate_tables['bone_table'].loc[
                max_index, 'Mass Attenuation [cm2/g]']  # noqa: E501
            self._muscle_mass_attenuation = self._interpolate_tables['muscle_table'].loc[
                max_index, 'Mass Attenuation [cm2/g]']  # noqa: E501

            self._K2HPO4_mass_attenuation = self._interpolate_tables['k2hpo4_table'].loc[
                max_index, 'Mass Attenuation [cm2/g]']  # noqa: E501
            k2hpo4.append(self._K2HPO4_mass_attenuation)
            self._CHA_mass_attenuation = self._interpolate_tables['cha_table'].loc[
                max_index, 'Mass Attenuation [cm2/g]']  # noqa: E501
            self._triglyceride_mass_attenuation = self._interpolate_tables['triglyceride_table'].loc[
                max_index, 'Mass Attenuation [cm2/g]']  # noqa: E501
            triglyceride.append(self._triglyceride_mass_attenuation)
            self._water_mass_attenuation = self._interpolate_tables['water_table'].loc[
                max_index, 'Mass Attenuation [cm2/g]']  # noqa: E501

        k2hpo4_np = np.array(k2hpo4)
        self._K2HPO4_mass_attenuation = np.mean(k2hpo4_np)
        self.k2hpo4_stdev = np.std(k2hpo4_np)
        self.k2hpo4_standard_error = self.k2hpo4_stdev / np.sqrt(self.iterations)
        triglyceride_np = np.array(triglyceride)
        self._triglyceride_mass_attenuation = np.mean(triglyceride_np)
        self.triglyceride_stddev = np.std(triglyceride_np)
        self.triglyceride_standard_error = self.triglyceride_stddev / np.sqrt(self.iterations)


    def montecarlo_determine_hu_to_mass_attenuation(self):
        """Compute HU to mass attenuation relationship"""
        # Get values
        HU = self._subset([
            self.adipose_hu,
            self.air_hu,
            self.blood_hu,
            self.bone_hu,
            self.muscle_hu
        ])

        attenuation = self._subset([
            self.adipose_mass_attenuation,
            self.air_mass_attenuation,
            self.blood_mass_attenuation,
            self.bone_mass_attenuation,
            self.muscle_mass_attenuation
        ])

        # Perform linear regression
        linreg = stats.linregress(HU, attenuation)

        # Store values
        self._hu_to_mass_attenuation_slope = linreg[0]
        self._hu_to_mass_attenuation_intercept = linreg[1]
        self._hu_to_mass_attenuation_slope_error = linreg.stderr
        self._hu_to_mass_attenuation_intercept_error = linreg.intercept_stderr

    def montecarlo_determine_hu_to_density(self):
        """Compute HU to density relationship"""
        # Get HU
        HU = self._subset([
            self.adipose_hu,
            self.air_hu,
            self.blood_hu,
            self.bone_hu,
            self.muscle_hu
        ])

        # Convert to densities
        def _convert_value(hu, attenuation):
            water_attenuation = self.water_mass_attenuation
            water_density = 1.0

            return (
                           hu / 1000 * water_attenuation * water_density + water_attenuation * water_density) / attenuation  # noqa: E501

        densities = self._subset([
            _convert_value(self.adipose_hu, self.adipose_mass_attenuation),
            _convert_value(self.air_hu, self.air_mass_attenuation),
            _convert_value(self.blood_hu, self.blood_mass_attenuation),
            _convert_value(self.bone_hu, self.bone_mass_attenuation),
            _convert_value(self.muscle_hu, self.muscle_mass_attenuation)
        ])

        # Perform linear regression
        linreg = stats.linregress(HU, densities)

        # Store values
        self._hu_to_density_slope = linreg[0]
        self._hu_to_density_intercept = linreg[1]
        self._hu_to_density_slope_error = linreg.stderr
        self._hu_to_density_intercept_error = linreg.intercept_stderr

    @property
    def adipose_hu(self):
        """Sampled adipose hounsfield units"""
        return self._adipose_hu

    @adipose_hu.setter
    def adipose_hu(self, HU):
        """Sampled adipose hounsfield units"""
        self._adipose_hu = HU

    @property
    def air_hu(self):
        """Sampled air hounsfield units"""
        return self._air_hu

    @air_hu.setter
    def air_hu(self, HU):
        """Sampled air hounsfield units"""
        self._air_hu = HU

    @property
    def blood_hu(self):
        """Sampled blood hounsfield units"""
        return self._blood_hu

    @blood_hu.setter
    def blood_hu(self, HU):
        """Sampled blood hounsfield units"""
        self._blood_hu = HU

    @property
    def bone_hu(self):
        """Sampled bone hounsfield units"""
        return self._bone_hu

    @bone_hu.setter
    def bone_hu(self, HU):
        """Sampled bone hounsfield units"""
        self._bone_hu = HU

    @property
    def muscle_hu(self):
        """Sampled muscle hounsfield units"""
        return self._muscle_hu

    @muscle_hu.setter
    def muscle_hu(self, HU):
        """Sampled muscle hounsfield units"""
        self._muscle_hu = HU

    @property
    def effective_energy(self):
        """Calculated effective energy [keV]"""
        return self._effective_energy

    @property
    def max_r2(self):
        """Maximum coefficient of determination"""
        return self._max_r2

    @property
    def adipose_mass_attenuation(self):
        """Adipose mass attenuation coefficient at the effective energy"""
        return self._adipose_mass_attenuation

    @property
    def air_mass_attenuation(self):
        """Air mass attenuation coefficient at the effective energy"""
        return self._air_mass_attenuation

    @property
    def blood_mass_attenuation(self):
        """Blood mass attenuation coefficient at the effective energy"""
        return self._blood_mass_attenuation

    @property
    def bone_mass_attenuation(self):
        """Bone mass attenuation coefficient at the effective energy"""
        return self._bone_mass_attenuation

    @property
    def muscle_mass_attenuation(self):
        """Muscle mass attenuation coefficient at the effective energy"""
        return self._muscle_mass_attenuation

    @property
    def K2HPO4_mass_attenuation(self):
        """K2HPO4 mass attenuation coefficient at the effective energy"""
        return self._K2HPO4_mass_attenuation

    @property
    def CHA_mass_attenuation(self):
        """CHA mass attenuation coefficient at the effective energy"""
        return self._CHA_mass_attenuation

    @property
    def triglyceride_mass_attenuation(self):
        """Triglyceride mass attenuation coefficient at the effective energy"""
        return self._triglyceride_mass_attenuation

    @property
    def water_mass_attenuation(self):
        """Water mass attenuation coefficient at the effective energy"""
        return self._water_mass_attenuation

    @property
    def hu_to_mass_attenuation_slope(self):
        """HU to mass attenuation slope"""
        return self._hu_to_mass_attenuation_slope

    @property
    def hu_to_mass_attenuation_intercept(self):
        """HU to mass attenuation intercept"""
        return self._hu_to_mass_attenuation_intercept

    @property
    def hu_to_density_slope(self):
        """HU to density slope"""
        return self._hu_to_density_slope

    @property
    def hu_to_density_intercept(self):
        """HU to density intercept"""
        return self._hu_to_density_intercept

    def get_dict(self):
        return OrderedDict([
            ('Is fit', self._is_fit),

            ('Adipose [HU]', self.adipose_hu),
            ('Air [HU]', self.air_hu),
            ('Blood [HU]', self.blood_hu),
            ('Bone [HU]', self.bone_hu),
            ('Muscle [HU]', self.muscle_hu),

            ('Voxel Volume [mm3]', self._voxel_volume),

            ('Effective Energy [keV]', self.effective_energy),
            ('Max R^2', self.max_r2),

            ('Adipose u/p [cm2/g]', self.adipose_mass_attenuation),
            ('Air u/p [cm2/g]', self.air_mass_attenuation),
            ('Blood u/p [cm2/g]', self.blood_mass_attenuation),
            ('Bone u/p [cm2/g]', self.bone_mass_attenuation),
            ('Muscle u/p [cm2/g]', self.muscle_mass_attenuation),
            ('K2HPO4 u/p [cm2/g]', self.K2HPO4_mass_attenuation),
            ('CHA u/p [cm2/g]', self.CHA_mass_attenuation),
            ('Triglyceride u/p [cm2/g]', self.triglyceride_mass_attenuation),
            ('Water u/p [cm2/g]', self.water_mass_attenuation),

            ('HU-u/p Slope', self.hu_to_mass_attenuation_slope),
            ('HU -u/p Slope Standard Error', self._hu_to_mass_attenuation_slope_error)
            ('HU-u/p Y-Intercept', self.hu_to_mass_attenuation_intercept),
            ('HU-u/p Y-Intercept Standard Error', self._hu_to_mass_attenuation_intercept_error)

            ('HU-Material Density Slope', self.hu_to_density_slope),
            ('HU-Material Denstiy Slope Standard Error', self._hu_to_density_slope_error)
            ('HU-Material Density Y-Intercept', self.hu_to_density_intercept)
            ('HU-Material Density Y-Intercept Standard Error', self._hu_to_density_intercept_error)
        ])