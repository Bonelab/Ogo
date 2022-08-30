'''Test the standard calibration class'''

import unittest
from ogo.calibration.standard_calibration import StandardCalibration
import os


class TestStandardCalibration(unittest.TestCase):
    '''Test the standard calibration class'''

    def test_can_instantiate(self):
        '''As an inherited class from an abstract class, make sure we can inherit'''
        c = StandardCalibration()

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
        hu = [490.28414063, 200.31946153,  68.02660165, -77.60183361, -77.74668932]
        d = [375.83, 157.05,  58.88, -53.4 , -51.83]

        c = StandardCalibration()
        c.fit(hu, d)

        d_predict = c.predict(hu)

        for i, value in enumerate([375.99529306, 157.36982588,  57.62461421, -52.17525794, -52.28447521]):
            self.assertAlmostEqual(d_predict[i], value)

    def test_predict_inverse(self):
        '''Correctly predicts the inverse for a basic example'''
        hu = [490.28414063, 200.31946153,  68.02660165, -77.60183361, -77.74668932]
        d = [375.83, 157.05,  58.88, -53.4 , -51.83]

        c = StandardCalibration()
        c.fit(hu, d)

        hu_predict = c.predict_inverse(d)

        for i, value in enumerate([490.06491114, 199.89527394,  69.69162972, -79.22621864, -77.14391528]):
            self.assertAlmostEqual(hu_predict[i], value)

if __name__ == '__main__':
    unittest.main()
