'''Test the ogoBmdAnalysis command line interface exists'''

import unittest
import subprocess


class TestBMDAnalysisExists(unittest.TestCase):
    '''Test the ogoBmdAnalysis command line interface exists'''

    def test_can_call_ogo_help(self):
        '''Can call `BmdAnalysis --help`'''
        command = ['ogoBmdAnalysis', '--help']
        self.assertTrue(subprocess.check_output(command))


if __name__ == '__main__':
    unittest.main()
