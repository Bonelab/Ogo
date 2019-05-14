'''Test the ogo command line interface for calib midways'''

import unittest
import subprocess


class TestMindwaysCLI(unittest.TestCase):
    '''Test the ogo command line interface for calib midways'''

    def test_can_call_mindways_help(self):
        '''Can call `ogo calib mindways --help`'''
        command = ['ogo', 'calib', 'mindways', '--help']
        self.assertTrue(subprocess.check_output(command))

    def test_blank_mindways_call_fails(self):
        '''Calling `ogo calib mindways` fails'''
        command = ['ogo', 'calib', 'mindways']
        with self.assertRaises(subprocess.CalledProcessError):
            subprocess.check_output(command, stderr=subprocess.STDOUT)

if __name__ == '__main__':
    unittest.main()
