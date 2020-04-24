'''Test write csv functionality'''

import unittest
from ogo.util.write_csv import write_csv
import subprocess
import shutil, tempfile, filecmp
import os


class TestWriteCSV(unittest.TestCase):
    '''Test write csv functionality'''

    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.csv_file_name = os.path.join(self.test_dir, 'test.csv')

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
        correct_file_name = os.path.join(self.test_dir, 'test_correct.csv')
        with open(correct_file_name, 'w') as f:
          f.write('FileName,A,B' + os.linesep)
          f.write('filename.txt,2,b' + os.linesep)

        # Print result
        write_csv(data, self.csv_file_name)
        self.assertTrue(os.path.exists(self.csv_file_name),
          'Could not create file')
        self.assertTrue(filecmp.cmp(correct_file_name, self.csv_file_name),
          'Files are not the sample')

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
          f.write('FileName;A;B' + os.linesep)
          f.write('filename.txt;2;b' + os.linesep)

        # Print result
        write_csv(data, self.csv_file_name, delimiter=';')
        self.assertTrue(os.path.exists(self.csv_file_name),
          'Could not create file')
        self.assertTrue(filecmp.cmp(correct_file_name, self.csv_file_name),
          'Files are not the same')

    def test_append_file(self):
        # Create test data
        data = {
          'FileName':   'filename.txt',
          'A':          2,
          'B':          'b'
        }

        # First write
        write_csv(data, self.csv_file_name)

        self.assertTrue(os.path.exists(self.csv_file_name),
          'Could not create file')

        # Second write
        write_csv(data, self.csv_file_name)

        # Create the expected file
        correct_file_name = os.path.join(self.test_dir, 'test_correct.csv')
        with open(correct_file_name, 'w') as f:
          f.write('FileName,A,B' + os.linesep)
          f.write('filename.txt,2,b' + os.linesep)
          f.write('filename.txt,2,b' + os.linesep)

        self.assertTrue(filecmp.cmp(correct_file_name, self.csv_file_name),
          'Files are not the same')


if __name__ == '__main__':
    unittest.main()
