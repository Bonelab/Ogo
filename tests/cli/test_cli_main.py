'''Test the ogo command line interface exists'''

import unittest
import subprocess


class TestCLIExists(unittest.TestCase):
    '''Test the ogo command line interface exists'''

    def test_can_call_ogo(self):
        '''Can call `ogo`'''
        command = ['ogo']
        self.assertTrue(subprocess.check_output(command))

    def test_can_call_ogo_help(self):
        '''Can call `ogo --help`'''
        command = ['ogo', '--help']
        self.assertTrue(subprocess.check_output(command))


if __name__ == '__main__':
    unittest.main()
