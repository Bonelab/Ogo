'''Test the ogo command line interface for calib exists'''

import unittest
import subprocess


class TestCLIForCalib(unittest.TestCase):
    '''Test the ogo command line interface for calib exists'''

    def test_can_call_calib(self):
        '''Can call `ogo calib`'''
        command = ['ogo', 'calib']
        self.assertTrue(subprocess.check_output(command))

    def test_can_call_calib_help(self):
        '''Can call `ogo calib --help`'''
        command = ['ogo', 'calib', '--help']
        self.assertTrue(subprocess.check_output(command))


if __name__ == '__main__':
    unittest.main()
