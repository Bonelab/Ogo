'''Test the ogo command line interface for calib internal'''

import unittest
import subprocess


class TestInternalCLI(unittest.TestCase):
    '''Test the ogo command line interface for calib internal'''

    def test_can_call_internal_help(self):
        '''Can call `ogo calib internal --help`'''
        command = ['ogo', 'calib', 'internal', '--help']
        self.assertTrue(subprocess.check_output(command))

    def test_blank_internal_call_fails(self):
        '''Calling `ogo calib internal` fails'''
        command = ['ogo', 'calib', 'internal']
        with self.assertRaises(subprocess.CalledProcessError):
            subprocess.check_output(command, stderr=subprocess.STDOUT)


if __name__ == '__main__':
    unittest.main()
