'''Test the mindways calibration class'''

import unittest
from ogo.calibration.mindways_calibration import MindwaysCalibration
import numpy as np
from collections import OrderedDict


class TestMindwaysCalibration(unittest.TestCase):
    '''Test the mindways calibration class'''

    def test_initial_values(self):
        '''Lock down the initial values to 0'''
        c = MindwaysCalibration()
        self.assertAlmostEqual(c.sigma_ref, 0)
        self.assertAlmostEqual(c.beta_ref, 0)
        self.assertAlmostEqual(c.sigma_ct, 0)
        self.assertAlmostEqual(c.beta_ct, 0)

    def test_predict_without_fit_throws_error(self):
        '''Predicting without fitting should throw an error'''
        c = MindwaysCalibration()

        with self.assertRaises(RuntimeError):
            c.predict(0.0)

        with self.assertRaises(RuntimeError):
            c.predict_inverse(0.0)

    def test_fit_with_different_lengths_throws_exception(self):
        '''If the vectors have different lengths, should throw an exception'''
        c = MindwaysCalibration()

        with self.assertRaises(RuntimeError):
            c.fit([0, 1], [0, 1], [0])

        with self.assertRaises(RuntimeError):
            c.fit([0, 1], [0], [0, 1])

        with self.assertRaises(RuntimeError):
            c.fit([0], [0, 1], [0, 1])

    def test_fit_with_length_zero_is_error(self):
        '''If the vectors have length zero, should throw an exception'''
        c = MindwaysCalibration()

        with self.assertRaises(RuntimeError):
            c.fit([], [], [])

    def test_fit(self):
        '''Correctly fits a basic example'''
        densities = [375.83, 157.05, 58.88, -53.40, -51.83]
        water = [923.20, 1119.52, 1103.57, 1056.95, 1012.25]
        hu = [462.93637565619906, 345.73809215407005, 176.71197512467555, -35.378190313384316, -69.54917011727602]

        c = MindwaysCalibration()
        c.fit(hu, densities, water)

        self.assertAlmostEqual(c.slope,     0.8002826626918712)
        self.assertAlmostEqual(c.intercept, 8.122906532849933)
        self.assertAlmostEqual(c.r_value,   0.9998533422632067)
        self.assertAlmostEqual(c.p_value,   2.1319733486784323e-06)
        self.assertAlmostEqual(c.std_err,   0.014506827937069183)

        self.assertAlmostEqual(c.sigma_ref, 1.46695849554)
        self.assertAlmostEqual(c.beta_ref,  -1009.75004687)
        self.assertAlmostEqual(c.sigma_ct,  1.24955849554)
        self.assertAlmostEqual(c.beta_ct,   -10.1500468666)

    def test_predict(self):
        '''Correctly predicts a basic example'''
        densities = np.array([375.83, 157.05, 58.88, -53.40, -51.83])
        water = np.array([923.20, 1119.52, 1103.57, 1056.95, 1012.25])
        hu = np.array([462.93637565619906, 345.73809215407005, 176.71197512467555, -35.378190313384316, -69.54917011727602])

        c = MindwaysCalibration()
        c.fit(hu, densities, water)

        d_predict = c.predict(hu)

        for i, value in enumerate([378.6028619 , 284.81110752, 149.54243652, -20.18964581, -47.53608852]):
            self.assertAlmostEqual(d_predict[i], value)

    def test_predict_inverse(self):
        '''Correctly predicts the inverse for a basic example'''
        densities = np.array([375.83, 157.05, 58.88, -53.40, -51.83])
        water = np.array([923.20, 1119.52, 1103.57, 1056.95, 1012.25])
        hu = np.array([462.93637565619906, 345.73809215407005, 176.71197512467555, -35.378190313384316, -69.54917011727602])

        c = MindwaysCalibration()
        c.fit(hu, densities, water)

        hu_predict = c.predict_inverse(densities)

        for i, value in enumerate([459.47152251, 186.09311486,  63.42395735, -76.87647053, -74.91466369]):
            self.assertAlmostEqual(hu_predict[i], value)

    def test_default_get_dict(self):
        '''get_dict returns dictionary with expected keys'''
        c = MindwaysCalibration()
        d = c.get_dict()

        self.assertTrue(OrderedDict, type(d))

        expected_keys = [
            'Is fit',
            'Slope',
            'Intercept',
            'R value',
            'P value',
            'Standard Error',
            'HU',
            'Densities',
            'Water',
            'Sigma ref',
            'Beta ref',
            'Sigma ct',
            'Beta ct'
        ]

        self.assertTrue(len(expected_keys), len(d))
        for a, b in zip(expected_keys, d.keys()):
            self.assertEqual(a, b)

    def test_get_dict_after_fit(self):
        '''get_dict set properly after fit'''
        densities = np.array([375.83, 157.05, 58.88, -53.40, -51.83])
        water = np.array([923.20, 1119.52, 1103.57, 1056.95, 1012.25])
        hu = np.array([462.93637565619906, 345.73809215407005, 176.71197512467555, -35.378190313384316, -69.54917011727602])

        c = MindwaysCalibration()
        c.fit(hu, densities, water)
        d = c.get_dict()

        self.assertEqual(True, d['Is fit'])
        self.assertAlmostEqual(0.8002826626887849, d['Slope'])
        self.assertAlmostEqual(8.122906532840492, d['Intercept'])
        self.assertAlmostEqual(0.9998533422632067, d['R value'])
        self.assertAlmostEqual(2.1319733486784323e-06, d['P value'])
        self.assertAlmostEqual(0.014506827937069183, d['Standard Error'])
        self.assertEqual(
           '[462.93637566 345.73809215 176.71197512 -35.37819031 -69.54917012]',
            d['HU']
        )
        self.assertEqual(
            '[375.83 157.05  58.88 -53.4  -51.83]',
            d['Densities']
        )
        self.assertEqual(
            '[ 923.2  1119.52 1103.57 1056.95 1012.25]',
            d['Water']
        )
        self.assertAlmostEqual(1.4669584955448192, d['Sigma ref'])
        self.assertAlmostEqual(-1009.7500468666274, d['Beta ref'])
        self.assertAlmostEqual(1.2495584955448191, d['Sigma ct'])
        self.assertAlmostEqual(-10.15004686662735, d['Beta ct'])


if __name__ == '__main__':
    unittest.main()
