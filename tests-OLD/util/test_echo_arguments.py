'''Test ogo setup version is correct'''

import unittest
from ogo.util.echo_arguments import echo_arguments
import os


class TestEchoArguments(unittest.TestCase):
    '''Test the ogo version is PEP compliant'''

    def test_empty_dictionary(self):
        expected = 'title:'+os.linesep
        self.assertEqual(echo_arguments('title', {}), expected)

    def test_simple_dictionary(self):
        expected = 'title:' + os.linesep + '  a  1' + os.linesep
        self.assertEqual(echo_arguments('title', {'a':1}), expected)

if __name__ == '__main__':
    unittest.main()
