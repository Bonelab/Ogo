'''Test the ogo command line interface for img exists'''

import unittest
import subprocess


class TestCLIForImag(unittest.TestCase):
    '''Test the ogo command line interface for img exists'''

    def test_can_call_img(self):
        '''Can call `ogo img`'''
        command = ['ogo', 'img']
        self.assertTrue(subprocess.check_output(command))

    def test_can_call_img_help(self):
        '''Can call `ogo img --help`'''
        command = ['ogo', 'img', '--help']
        self.assertTrue(subprocess.check_output(command))


if __name__ == '__main__':
    unittest.main()
