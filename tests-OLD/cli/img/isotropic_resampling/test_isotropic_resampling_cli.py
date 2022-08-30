'''Test the command line interface for resample'''

import unittest
import subprocess


class TestIsotropicResamplingCLI(unittest.TestCase):
    '''Test the command line interface for resample'''

    def test_can_call_resample_help(self):
        '''Can call `ogo img resample --help`'''
        command = ['ogo', 'img', 'resample', '--help']
        self.assertTrue(subprocess.check_output(command))

    def test_blank_resample_call_fails(self):
        '''Calling `ogo img resample` fails'''
        command = ['ogo', 'img', 'resample']
        with self.assertRaises(subprocess.CalledProcessError):
            subprocess.check_output(command, stderr=subprocess.STDOUT)

if __name__ == '__main__':
    unittest.main()
