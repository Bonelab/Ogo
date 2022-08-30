import unittest

from ogo.calib.calibration import Calibration


class Test_Calibration(unittest.TestCase):

    def setUp(self) -> None:
        """ Use this method to prep data or files for other tests, if necessary. """
        pass

    def tearDown(self) -> None:
        """ Use this method to clean up after the tests, if necessary."""
        pass

    def test_cannot_instantiate_abc(self) -> None:
        """`Calibration` is an ABC, so we cannot instantiate it"""
        with self.assertRaises(TypeError):
            Calibration()

    @unittest.skip("Placeholder - add more tests here")
    def test_placeholder(self) -> None:
        """ ogo.calib.calibration """
        pass
