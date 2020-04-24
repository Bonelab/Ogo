'''Test the base abstract class Calibration'''

import unittest
from ogo.calibration.calibration import Calibration


class TestCalibration(unittest.TestCase):
    '''Test the base abstract class Calibration'''

    def test_cannot_instantiate(self):
        with self.assertRaises(TypeError):
            Calibration()


if __name__ == '__main__':
    unittest.main()
