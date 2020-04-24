'''Test internal calibration'''

import unittest
from ogo.calibration.internal_calibration import \
    InternalCalibration
import numpy as np
import numpy.testing as npt
from collections import OrderedDict


class TestInternalCalibration(unittest.TestCase):
    '''Test the mass attenuation tables'''

    def test_default_values(self):
        '''Correct default values'''
        m = InternalCalibration()

        self.assertAlmostEqual(0.0, m.adipose_hu)
        self.assertAlmostEqual(0.0, m.air_hu)
        self.assertAlmostEqual(0.0, m.blood_hu)
        self.assertAlmostEqual(0.0, m.bone_hu)
        self.assertAlmostEqual(0.0, m.muscle_hu)

        self.assertAlmostEqual(0.0, m.effective_energy)
        self.assertAlmostEqual(0.0, m.max_r2)

        self.assertAlmostEqual(0.0, m.adipose_mass_attenuation)
        self.assertAlmostEqual(0.0, m.air_mass_attenuation)
        self.assertAlmostEqual(0.0, m.blood_mass_attenuation)
        self.assertAlmostEqual(0.0, m.bone_mass_attenuation)
        self.assertAlmostEqual(0.0, m.muscle_mass_attenuation)
        self.assertAlmostEqual(0.0, m.K2HPO4_mass_attenuation)
        self.assertAlmostEqual(0.0, m.CHA_mass_attenuation)
        self.assertAlmostEqual(0.0, m.triglyceride_mass_attenuation)
        self.assertAlmostEqual(0.0, m.water_mass_attenuation)

        self.assertAlmostEqual(0.0, m.hu_to_mass_attenuation_slope)
        self.assertAlmostEqual(0.0, m.hu_to_mass_attenuation_intercept)

        self.assertAlmostEqual(0.0, m.hu_to_density_slope)
        self.assertAlmostEqual(0.0, m.hu_to_density_intercept)

    def test_constructor_setting(self):
        '''Setting through constructor works'''
        m = InternalCalibration(
            adipose_hu=-100.0,
            air_hu=-10.0,
            blood_hu=10.0,
            bone_hu=100.0,
            muscle_hu=1000.0
        )

        self.assertAlmostEqual(-100.0, m.adipose_hu)
        self.assertAlmostEqual(-10.0, m.air_hu)
        self.assertAlmostEqual(10.0, m.blood_hu)
        self.assertAlmostEqual(100.0, m.bone_hu)
        self.assertAlmostEqual(1000.0, m.muscle_hu)

    def test_set_get_adipose_hu(self):
        '''Set/Get adipose_hu'''
        m = InternalCalibration()
        self.assertAlmostEqual(0.0, m.adipose_hu)
        m.adipose_hu = -10.2
        self.assertAlmostEqual(-10.2, m.adipose_hu)

    def test_set_get_air_hu(self):
        '''Set/Get air_hu'''
        m = InternalCalibration()
        self.assertAlmostEqual(0.0, m.air_hu)
        m.air_hu = -10.2
        self.assertAlmostEqual(-10.2, m.air_hu)

    def test_set_get_blood_hu(self):
        '''Set/Get blood_hu'''
        m = InternalCalibration()
        self.assertAlmostEqual(0.0, m.blood_hu)
        m.blood_hu = -10.2
        self.assertAlmostEqual(-10.2, m.blood_hu)

    def test_set_get_bone_hu(self):
        '''Set/Get bone_hu'''
        m = InternalCalibration()
        self.assertAlmostEqual(0.0, m.bone_hu)
        m.bone_hu = -10.2
        self.assertAlmostEqual(-10.2, m.bone_hu)

    def test_set_get_muscle_hu(self):
        '''Set/Get muscle_hu'''
        m = InternalCalibration()
        self.assertAlmostEqual(0.0, m.muscle_hu)
        m.muscle_hu = -10.2
        self.assertAlmostEqual(-10.2, m.muscle_hu)

    def test_fit_data(self):
        '''Correctly fits example data'''
        # Data
        hu = {
            'adipose_hu':   -97.7360,
            'air_hu':       -988.9582,
            'blood_hu':     29.8856,
            'bone_hu':      1361.9773,
            'muscle_hu':    53.5251
        }

        # Create object
        calib = InternalCalibration(**hu)

        # Calibrate
        calib.fit()

        # Inputs still set correctly
        self.assertAlmostEqual(-97.7360,    calib.adipose_hu)
        self.assertAlmostEqual(-988.9582,   calib.air_hu)
        self.assertAlmostEqual(29.8856,     calib.blood_hu)
        self.assertAlmostEqual(1361.9773,   calib.bone_hu)
        self.assertAlmostEqual(53.5251,     calib.muscle_hu)

        # Test the output values
        self.assertAlmostEqual(96.0,                calib.effective_energy)
        self.assertAlmostEqual(0.9983404998928974,  calib.max_r2)

        self.assertAlmostEqual(0.17104, calib.adipose_mass_attenuation)
        self.assertAlmostEqual(0.15652, calib.air_mass_attenuation)
        self.assertAlmostEqual(0.17214, calib.blood_mass_attenuation)
        self.assertAlmostEqual(0.19298, calib.bone_mass_attenuation)
        self.assertAlmostEqual(0.1719,  calib.muscle_mass_attenuation)

        self.assertAlmostEqual(0.20792, calib.K2HPO4_mass_attenuation)
        self.assertAlmostEqual(0.21330, calib.CHA_mass_attenuation)
        self.assertAlmostEqual(0.17042, calib.triglyceride_mass_attenuation)
        self.assertAlmostEqual(0.17330, calib.water_mass_attenuation)

        self.assertAlmostEqual(1.5474534088249448e-05,  calib.hu_to_mass_attenuation_slope)
        self.assertAlmostEqual(0.17180587621547028,     calib.hu_to_mass_attenuation_intercept)

        self.assertAlmostEqual(0.0008884964716273846,   calib.hu_to_density_slope)
        self.assertAlmostEqual(0.9655496613171292,      calib.hu_to_density_intercept)

    def test_predict_data(self):
        '''Correctly predicts example data'''
        # Data
        hu = {
            'adipose_hu':   -97.7360,
            'air_hu':       -988.9582,
            'blood_hu':     29.8856,
            'bone_hu':      1361.9773,
            'muscle_hu':    53.5251
        }
        voxel_volume_mm = 0.48925677750966123

        # Create object
        calib = InternalCalibration(**hu)

        # Calibrate
        calib.fit()

        # Predict
        expected = np.array([
            2084.294978, 55.369403, -12.496955, 46.245288
        ])
        hu = np.array([
            -3024.0, 44.0, -125.0, 24.0
        ])
        predicted = calib.predict(hu, voxel_volume_mm)

        npt.assert_array_almost_equal(expected, predicted)

    def test_default_get_dict(self):
        '''get_dict returns dictionary with expected keys'''
        c = InternalCalibration()
        d = c.get_dict()

        self.assertTrue(OrderedDict, type(d))

        expected_keys = [
            'Is fit',
            'Adipose [HU]',
            'Air [HU]',
            'Blood [HU]',
            'Bone [HU]',
            'Muscle [HU]',
            'Voxel Volume [mm3]',
            'Effective Energy [keV]',
            'Max R^2',
            'Adipose u/p [cm2/g]',
            'Air u/p [cm2/g]',
            'Blood u/p [cm2/g]',
            'Bone u/p [cm2/g]',
            'Muscle u/p [cm2/g]',
            'K2HPO4 u/p [cm2/g]',
            'CHA u/p [cm2/g]',
            'Triglyceride u/p [cm2/g]',
            'Water u/p [cm2/g]',
            'HU-u/p Slope',
            'HU-u/p Y-Intercept',
            'HU-Material Density Slope',
            'HU-Material Density Y-Intercept'
        ]

        self.assertTrue(len(expected_keys), len(d))
        for a, b in zip(expected_keys, d.keys()):
            self.assertEqual(a, b)

    def test_get_dict_after_predict(self):
        '''get_dict set properly after fit'''
        # Data
        hu = {
            'adipose_hu':   -97.7360,
            'air_hu':       -988.9582,
            'blood_hu':     29.8856,
            'bone_hu':      1361.9773,
            'muscle_hu':    53.5251
        }

        # Create object
        calib = InternalCalibration(**hu)

        # Calibrate
        calib.fit()
        calib.predict(1.0, 0.5)
        d = calib.get_dict()

        self.assertEqual(True, d['Is fit'])

        self.assertAlmostEqual(-97.736,     d['Adipose [HU]'])
        self.assertAlmostEqual(-988.9582,   d['Air [HU]'])
        self.assertAlmostEqual(29.8856,     d['Blood [HU]'])
        self.assertAlmostEqual(1361.9773,   d['Bone [HU]'])
        self.assertAlmostEqual(53.5251,     d['Muscle [HU]'])

        self.assertEqual(0.5,   d['Voxel Volume [mm3]'])

        self.assertAlmostEqual(96.0,                d['Effective Energy [keV]'])
        self.assertAlmostEqual(0.9983405001630987,  d['Max R^2'])

        self.assertAlmostEqual(0.17104,     d['Adipose u/p [cm2/g]'])
        self.assertAlmostEqual(0.15652,     d['Air u/p [cm2/g]'])
        self.assertAlmostEqual(0.17214,     d['Blood u/p [cm2/g]'])
        self.assertAlmostEqual(0.19298,     d['Bone u/p [cm2/g]'])
        self.assertAlmostEqual(0.17190,     d['Muscle u/p [cm2/g]'])

        self.assertAlmostEqual(0.20792,     d['K2HPO4 u/p [cm2/g]'])
        self.assertAlmostEqual(0.21330,     d['CHA u/p [cm2/g]'])
        self.assertAlmostEqual(0.17042,     d['Triglyceride u/p [cm2/g]'])
        self.assertAlmostEqual(0.17330,     d['Water u/p [cm2/g]'])

        self.assertAlmostEqual(1.5474534530872792e-05,  d['HU-u/p Slope'])
        self.assertAlmostEqual(0.17180587608117798,     d['HU-u/p Y-Intercept'])
        self.assertAlmostEqual(0.0008884964685912277,   d['HU-Material Density Slope'])
        self.assertAlmostEqual(0.9655496637163411,      d['HU-Material Density Y-Intercept'])


if __name__ == '__main__':
    unittest.main()
