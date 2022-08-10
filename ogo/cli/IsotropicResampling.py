# /------------------------------------------------------------------------------+
# | 19-JUL-2022                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+

# Imports
import argparse
import os
import SimpleITK as sitk
import numpy as np
from ogo.util.echo_arguments import echo_arguments
import ogo.util.Helper as ogo

def IsotropicResampling(input_filename, output_filename, iso_resolution, overwrite=False):

    # Check if output exists and should overwrite
    if os.path.isfile(output_filename) and not overwrite:
        result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(output_filename))
        if result.lower() not in ['y', 'yes']:
            print('Not overwriting. Exiting...')
            os.sys.exit()
    
    # Read input
    if not os.path.isfile(input_filename):
        os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(input_filename))
    image = sitk.ReadImage(input_filename)
    
    # Determine sigma
    ogo.message('Resampling to a resolution of {}'.format(iso_resolution))
    sigma = iso_resolution * np.sqrt(2 * np.log(2)) / np.pi
    ogo.message('  For a half-power cut-off frequency of {:0.6f}, '
               'a standard deviation of {:0.6f} is being used.'.format(1/(2.0*iso_resolution), sigma))
    
    # Filter each sampling directions as needed
    for i, spacing in enumerate(image.GetSpacing()):
        if spacing > iso_resolution:
            ogo.message('  No antialiasing in direction {}'.format(i))
            continue

        ogo.message('  Antialiasing in direction {}'.format(i))
        image = sitk.RecursiveGaussian(
          image,
          sigma, False,
          sitk.RecursiveGaussianImageFilter.ZeroOrder,
          i
        )
        
    # Determine the output size
    resolution = [iso_resolution for i in range(image.GetDimension())]
    size = [int(np.ceil(s * i / o)) for o, i, s in
            zip(resolution, image.GetSpacing(), image.GetSize())]
    ogo.message('  Input Size:  {}'.format(image.GetSize()))
    ogo.message('  Output Size: {}'.format(size))

    ogo.message('Finding minimum in image')
    stats = sitk.StatisticsImageFilter()
    stats.Execute(image)
    vox_min = stats.GetMinimum()
    ogo.message('  Minimum intensity: {}'.format(vox_min))

    ogo.message('Resampling...')
    transform = sitk.Euler3DTransform()
    transform.SetIdentity()

    ogo.message('  Using BSpline interpolation')
    
    output = sitk.Resample(
      image,
      size,
      transform,
      sitk.sitkBSpline,
      image.GetOrigin(),
      resolution,
      image.GetDirection(),
      vox_min,
      image.GetPixelID()
    )

    ogo.message('Writing output to file {}'.format(output_filename))
    sitk.WriteImage(output, output_filename)    
    
def main():
    # Setup description
    description='''
Utility to resample a CT image to an isotropic resolution. It resamples using 
BSpline interpolation.
'''
    epilog='''
Example call: 
ogoIsotropicResampling input.nii output.nii --iso_resolution 1.0
'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoIsotropicResampling",
        description=description,
        epilog=epilog
    )
    parser.add_argument('input_filename', help='Input image file (*.nii, *.nii.gz)')
    parser.add_argument('output_filename', help='Output image file (*.nii, *.nii.gz)')
    parser.add_argument('--iso_resolution', type=float, default=1.0, metavar='DIM', help='Isotropic resolution (default: %(default)s)')
    parser.add_argument('-o', '--overwrite', action='store_true', help='Overwrite output without asking')

    print()

    # Parse and display
    args = parser.parse_args()
    print(echo_arguments('IsotropicResampling', vars(args)))

    # Run program
    IsotropicResampling(**vars(args))

if __name__ == '__main__':
    main()
