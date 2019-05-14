'''Test the ogo command line interface for calib standard'''

import unittest
import subprocess


class TestStandardCLI(unittest.TestCase):
    '''Test the ogo command line interface for calib standard'''

    def test_can_call_standard_help(self):
        '''Can call `ogo calib standard --help`'''
        command = ['ogo', 'calib', 'standard', '--help']
        self.assertTrue(subprocess.check_output(command))

    def test_blank_standard_call_fails(self):
        '''Calling `ogo calib standard` fails'''
        command = ['ogo', 'calib', 'standard']
        with self.assertRaises(subprocess.CalledProcessError):
            subprocess.check_output(command, stderr=subprocess.STDOUT)

if __name__ == '__main__':
    unittest.main()
