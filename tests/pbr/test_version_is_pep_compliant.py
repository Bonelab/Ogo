'''Test ogo setup version is correct'''

import unittest
from packaging import version
import re
import ogo

class TestVersionIsCompliant(unittest.TestCase):
    '''Test the ogo version is PEP compliant'''

    def test_get_version(self):
        '''Can access __version__'''
        ogo.__version__

    def test_get_version_info(self):
        '''Can access version_info'''
        ogo.version_info

    def test_version_is_pep_compliant(self):
        '''__version__ is PEP 440 compliant'''
        # See https://www.python.org/dev/peps/pep-0440/#id79
        regex = re.compile( 
            r"^\s*" + version.VERSION_PATTERN + r"\s*$",
            re.VERBOSE | re.IGNORECASE,
        )
        self.assertTrue(regex.match(ogo.__version__))

if __name__ == '__main__':
    unittest.main()