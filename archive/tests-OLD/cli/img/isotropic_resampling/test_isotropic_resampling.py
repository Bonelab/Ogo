'''Test the logic for the function `ogo img resample`'''

import unittest
import subprocess
import shutil, tempfile
import SimpleITK as sitk
import os


class TestIsotropicResampling(unittest.TestCase):
    '''Test the logic for the function `ogo img resample`'''

    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.input_file = os.path.join(self.test_dir, 'input.nii')
        self.output_file = os.path.join(self.test_dir, 'output.nii')
    
        # Create image
        image = sitk.GaussianSource(
                    sitk.sitkUInt8,
                    [10, 10, 10],       # size
                    [1.0, 1.0, 1.0],    # sigma
                    [2.5, 5.0, 3.75],   # mean
                    255,                # scale
                    [0.0, 0.0, 0.0],    # origin
                    [0.5, 1.0, 0.75]    # spacing
                )
        sitk.WriteImage(image, self.input_file)

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_input_image_exists(self):
        self.assertTrue(os.path.isfile(self.input_file))

    def test_can_call_resample_help(self):
        '''Can call `ogo img resample --help`'''
        command = ['ogo', 'img', 'resample', self.input_file, self.output_file]
        self.assertTrue(subprocess.check_output(command))
        self.assertTrue(os.path.isfile(self.output_file))

if __name__ == '__main__':
    unittest.main()
