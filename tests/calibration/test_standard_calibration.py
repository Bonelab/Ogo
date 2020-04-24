'''Test the standard calibration class'''

import unittest
from ogo.calibration.standard_calibration import StandardCalibration
import os
from collections import OrderedDict
import numpy as np


class TestStandardCalibration(unittest.TestCase):
    '''Test the standard calibration class'''

    def test_can_instantiate(self):
        '''As an inherited class from an abstract class, make sure we can inherit'''
        StandardCalibration()

    def test_initial_values(self):
        '''Lock down the initial values to 0'''
        c = StandardCalibration()
        self.assertAlmostEqual(c.slope,     0)
        self.assertAlmostEqual(c.intercept, 0)
        self.assertAlmostEqual(c.r_value,   0)
        self.assertAlmostEqual(c.p_value,   0)
        self.assertAlmostEqual(c.std_err,   0)

    def test_predict_without_fit_throws_error(self):
        '''Predicting without fitting should throw an error'''
        c = StandardCalibration()

        with self.assertRaises(RuntimeError):
            c.predict(0.0)

        with self.assertRaises(RuntimeError):
            c.predict_inverse(0.0)

    def test_fit_with_different_lengths_throws_exception(self):
        '''If the vectors have different lengths, should throw an exception'''
        c = StandardCalibration()

        with self.assertRaises(RuntimeError):
            c.fit([0, 1], [0])

    def test_fit_with_length_zero_is_error(self):
        '''If the vectors have length zero, should throw an exception'''
        c = StandardCalibration()

        with self.assertRaises(RuntimeError):
            c.fit([], [])

    def test_fit(self):
        '''Correctly fits a basic example'''
        hu = [490.28414063, 200.31946153,  68.02660165, -77.60183361, -77.74668932]
        d = [375.83, 157.05,  58.88, -53.4 , -51.83]

        c = StandardCalibration()
        c.fit(hu, d)

        self.assertAlmostEqual(c.slope,     0.7539727523297295)
        self.assertAlmostEqual(c.intercept, 6.334410127001945)
        self.assertAlmostEqual(c.r_value,   0.9999866268559224)
        self.assertAlmostEqual(c.p_value,   5.870610556263143e-08)
        self.assertAlmostEqual(c.std_err,   0.0022512884297108417)

    def test_predict(self):
        '''Correctly predicts a basic example'''
        hu = np.array([490.28414063, 200.31946153,  68.02660165, -77.60183361, -77.74668932])
        d = np.array([375.83, 157.05,  58.88, -53.4 , -51.83])

        c = StandardCalibration()
        c.fit(hu, d)

        d_predict = c.predict(hu)

        for i, value in enumerate([375.99529306, 157.36982588,  57.62461421, -52.17525794, -52.28447521]):
            self.assertAlmostEqual(d_predict[i], value)

    def test_predict_inverse(self):
        '''Correctly predicts the inverse for a basic example'''
        hu = np.array([490.28414063, 200.31946153,  68.02660165, -77.60183361, -77.74668932])
        d = np.array([375.83, 157.05,  58.88, -53.4 , -51.83])

        c = StandardCalibration()
        c.fit(hu, d)

        hu_predict = c.predict_inverse(d)

        for i, value in enumerate([490.06491114, 199.89527394,  69.69162972, -79.22621864, -77.14391528]):
            self.assertAlmostEqual(hu_predict[i], value)

    def test_default_get_dict(self):
        '''get_dict returns dictionary with expected keys'''
        c = StandardCalibration()
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
        ]

        self.assertTrue(len(expected_keys), len(d))
        for a, b in zip(expected_keys, d.keys()):
            self.assertEqual(a, b)

    def test_get_dict_after_fit(self):
        '''get_dict set properly after fit'''
        hu = np.array([490.28414063, 200.31946153,  68.02660165, -77.60183361, -77.74668932])
        d = np.array([375.83, 157.05,  58.88, -53.4 , -51.83])

        c = StandardCalibration()
        c.fit(hu, d)
        d = c.get_dict()

        self.assertEqual(True, d['Is fit'])
        self.assertAlmostEqual(0.7539727523331842, d['Slope'])
        self.assertAlmostEqual(6.3344101269433395, d['Intercept'])
        self.assertAlmostEqual(0.9999866268559309, d['R value'])
        self.assertAlmostEqual(5.870610556263143e-08, d['P value'])
        self.assertAlmostEqual(0.0022512884297108417, d['Standard Error'])
        self.assertEqual(
            '[490.28414063 200.31946153  68.02660165 -77.60183361 -77.74668932]',
            d['HU']
        )
        self.assertEqual(
            '[375.83 157.05  58.88 -53.4  -51.83]',
            d['Densities']
        )


if __name__ == '__main__':
    unittest.main()
