'''Test write text functionality'''

import unittest
from ogo.util.write_txt import write_txt
import subprocess
import shutil, tempfile, filecmp
import os


class TestWriteTXT(unittest.TestCase):
    '''Test write text functionality'''

    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.txt_file_name = os.path.join(self.test_dir, 'test.txt')

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_create_file(self):
        # Create test data
        data = {
          'FileName':   'filename.txt',
          'A':          2,
          'B':          'b'
        }

        # Create the expected file
        correct_file_name = os.path.join(self.test_dir, 'test_correct.txt')
        with open(correct_file_name, 'w') as f:
            f.write('FileName,filename.txt' + os.linesep)
            f.write('A,2' + os.linesep)
            f.write('B,b' + os.linesep)

        # Print result
        write_txt(data, self.txt_file_name)
        self.assertTrue(os.path.exists(self.txt_file_name),
            'Could not create file')
        self.assertTrue(filecmp.cmp(correct_file_name, self.txt_file_name),
            'Files are not the same')

    def test_delimiter(self):
        # Create test data
        data = {
          'FileName':   'filename.txt',
          'A':          2,
          'B':          'b'
        }

        # Create the expected file
        correct_file_name = os.path.join(self.test_dir, 'test_correct.csv')
        with open(correct_file_name, 'w') as f:
            f.write('FileName;filename.txt' + os.linesep)
            f.write('A;2' + os.linesep)
            f.write('B;b' + os.linesep)

        # Print result
        write_txt(data, self.txt_file_name, delimiter=';')
        self.assertTrue(os.path.exists(self.txt_file_name),
          'Could not create file')
        self.assertTrue(filecmp.cmp(correct_file_name, self.txt_file_name),
          'Files are not the same')

    def test_append_file(self):
        # Create test data
        data = {
          'FileName':   'filename.txt',
          'A':          2,
          'B':          'b'
        }

        # First write
        write_txt(data, self.txt_file_name)

        self.assertTrue(os.path.exists(self.txt_file_name),
          'Could not create file')

        # Second write
        write_txt(data, self.txt_file_name, mode='a')

        # Create the expected file
        correct_file_name = os.path.join(self.test_dir, 'test_correct.csv')
        with open(correct_file_name, 'w') as f:
            f.write('FileName,filename.txt' + os.linesep)
            f.write('A,2' + os.linesep)
            f.write('B,b' + os.linesep)
            f.write('FileName,filename.txt' + os.linesep)
            f.write('A,2' + os.linesep)
            f.write('B,b' + os.linesep)

        self.assertTrue(filecmp.cmp(correct_file_name, self.txt_file_name),
          'Files are not the same')


if __name__ == '__main__':
    unittest.main()
