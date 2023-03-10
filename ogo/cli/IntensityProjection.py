# /------------------------------------------------------------------------------+
# | 09-MAR-2023                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+

script_version = 1.0

# Imports
import argparse
import os
import sys
import vtk
import math
import numpy as np
import SimpleITK as sitk
from scipy.spatial import procrustes
from datetime import date
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from ogo.util.echo_arguments import echo_arguments
import ogo.util.Helper as ogo
import ogo.dat.OgoMasterLabels as lb

def get_labels(ct):
    filt = sitk.LabelShapeStatisticsImageFilter()
    filt.Execute(ct)
    labels = filt.GetLabels()
    return labels

# +------------------------------------------------------------------------------+
def IntensityProjection(input_image, mask, output_image, projection_labels, projection_type, rotation, overwrite):

    # Check if output exists and should overwrite
    if output_image:
        if os.path.isfile(output_image) and not overwrite:
            result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(output_image))
            if result.lower() not in ['y', 'yes']:
                ogo.message('Not overwriting. Exiting...')
                os.sys.exit()

    #rotation_center = (0, 0, 0)
    #axis = (0,0,1)
    #angle = np.pi * rotation / 180.
    #translation = (0,0,0)
    #scale_factor = 1.0
    #similarity = sitk.Similarity3DTransform(scale_factor, axis, angle, translation, rotation_center)
    #
    #transform = sitk.AffineTransform(3)
    #transform.SetMatrix(similarity.GetMatrix())
    #transform.SetTranslation(similarity.GetTranslation())
    #transform.SetCenter(similarity.GetCenter())
    #
    ##transform.Rotate(0,0,radians)
    #print('dimension: {}'.format(transform.GetDimension()))
    #print('angle (deg): {}'.format(rotation))
    #print('angle (rad): {}'.format(angle))
    #print('transform: \n{}'.format(transform))

    # Read input
    if not os.path.isfile(input_image):
        os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(input_image))

    if not (input_image.lower().endswith('.nii') or input_image.lower().endswith('.nii.gz')):
        os.sys.exit('[ERROR] Input must be type NIFTI file: \"{}\"'.format(input_image))
    
    ogo.message('Reading image: ')
    ogo.message('\"{}\"'.format(input_image))
    ct = sitk.ReadImage(input_image, sitk.sitkUInt8)
    
    # Check if mask is provided and set things up
    if mask:
        # Read mask
        if not os.path.isfile(mask):
            os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(mask))
        if not (mask.lower().endswith('.nii') or mask.lower().endswith('.nii.gz')):
            os.sys.exit('[ERROR] Mask must be type NIFTI file: \"{}\"'.format(mask))
        ogo.message('Reading mask: ')
        ogo.message('\"{}\"'.format(mask))
        mask = sitk.ReadImage(mask, sitk.sitkUInt8)
        labels = get_labels(mask)
        ogo.message('Mask image contains the following labels:')
        ogo.message('  [' + ', '.join('{:d}'.format(i) for i in labels) + ']')
        if not projection_labels:
            projection_labels = labels
        else:
            for label in projection_labels:
                if label not in labels:
                    os.sys.exit('[ERROR] Label {} is not in input mask.'.format(label))
        ogo.message('Projections will be calculated for the following labels:')
        ogo.message('  [' + ', '.join('{:d}'.format(i) for i in projection_labels) + ']')
        
    else:
        ogo.message('No mask for quantitative analysis provided.')

    #BinaryProjectionImageFilter

    # Set the morphological operation
    ogo.message('Setting projection type to: [{}]'.format(projection_type))
    if projection_type == 'median':
        projectionFilt = sitk.MedianProjectionImageFilter()
    elif projection_type == 'mean':
        projectionFilt = sitk.MeanProjectionImageFilter()
    elif projection_type == 'minimum':
        projectionFilt = sitk.MinimumProjectionImageFilter()
    elif projection_type == 'stdev':
        projectionFilt = sitk.StandardDeviationProjectionImageFilter()
    elif projection_type == 'maximum':
        projectionFilt = sitk.MaximumProjectionImageFilter()
    elif projection_type == 'sum':
        projectionFilt = sitk.SumProjectionImageFilter()
    else:
        os.sys.exit('[ERROR] Unknown projection type: {}'.format(projection_type))
    
    # Set up transform
    rotation_center = (0, 0, 0)
    axis = (0,0,1) # Z-axis
    angle = np.pi * rotation / 180.
    translation = (0,0,0)
    scale_factor = 1.0
    similarity = sitk.Similarity3DTransform(scale_factor, axis, angle, translation, rotation_center)

    transform = sitk.AffineTransform(3)
    transform.SetMatrix(similarity.GetMatrix())
    transform.SetTranslation(similarity.GetTranslation())
    transform.SetCenter(similarity.GetCenter())
    
    # interpolation options are sitk.sitkBSpline, sitk.sitkNearestNeighbor, sitk.sitkLinear
    ogo.message('Applying transform of Z-axis rotation of {:.2f} degrees.'.format(rotation))
    ct_transformed = sitk.Resample(ct, transform, sitk.sitkBSpline)
    
    # Perform projection
    ogo.message('Performing intensity projection')
    projection = projectionFilt.Execute(ct_transformed)
    
    # Apply mask if suppled
    if mask:
        mask_transformed = sitk.Resample(mask, transform, sitk.sitkNearestNeighbor)
        
        stats = sitk.LabelIntensityStatisticsImageFilter()
        
        for label in projection_labels:
            try:
                desc = lb.labels_dict[label]['LABEL']
            except KeyError:
                desc = 'unknown label'
            ogo.message('  processing label {:>2d} ({})'.format(label,desc))
            
            stats.Execute(mask_transformed,ct_transformed)
            print('mean = {}'.format(stats.GetMean(label)))
    
    
    if output_image:
        ogo.message('Writing projection file:')
        ogo.message('  {}'.format(output_image))
        sitk.WriteImage(projection, output_image)
    else:
        ogo.message('No output projection image written to file.')            
    
    ogo.message('Done.')
    
    
def main():
    # Setup description
    description = '''
Generate an intensity projection of a 3D dataset. 

User has the option to select the rotation in degrees around the Z-axis 
(longitudinal axis) and to define a mask within which to generate an output.
'''

    epilog = '''
Example call: 
     
ogoIntensityProjection image.nii.gz --mask mask.nii.gz --rotation 10

ogoIntensityProjection /Users/skboyd/Desktop/ML/test/retro.nii.gz --mask /Users/skboyd/Desktop/ML/test/retro_mask.nii.gz --labels 7 8 9 10 --overwrite

'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoIntensityProjection",
        description=description,
        epilog=epilog
    )

    parser.add_argument('input_image', help='Input image file (*.nii, *.nii.gz)')
    parser.add_argument('--mask', default=None, metavar='MASK', help='Labels for quantitative output (default: %(default)s)')
    parser.add_argument('--output_image', default=None, metavar='OUTPUT', help='Output image file (default: %(default)s)')
    parser.add_argument('--projection_labels', type=int, nargs='*', default=[], metavar='LABEL', help='List of labels in mask to use (default: all)')
    parser.add_argument('--projection_type', default='sum', choices=['median', 'mean', 'maximum', 'minimum', 'stdev', 'sum'],
                         help='Select type of projection (default: %(default)s)')
    parser.add_argument('--rotation', type=float, default=90.0, metavar='ROTZ',help='Rotation in degrees about Z axis (default: %(default)s deg)')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite output file without asking')

    # Parse and display
    args = parser.parse_args()
    print(echo_arguments('IntensityProjection', vars(args)))

    # Run program
    IntensityProjection(**vars(args))


if __name__ == '__main__':
    main()
