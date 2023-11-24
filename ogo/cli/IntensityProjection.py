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
import yaml
import SimpleITK as sitk
from scipy.spatial import procrustes
from datetime import date
from datetime import datetime
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from ogo.util.echo_arguments import echo_arguments
import ogo.util.Helper as ogo
import ogo.dat.OgoMasterLabels as lb

def get_labels(ct):
    filt = sitk.LabelShapeStatisticsImageFilter()
    filt.Execute(ct)
    labels = filt.GetLabels()
    return labels

# This function is from https://github.com/rock-learning/pytransform3d/blob/7589e083a50597a75b12d745ebacaa7cc056cfbd/pytransform3d/rotations.py#L302
def matrix_from_axis_angle(a):
    """ Compute rotation matrix from axis-angle.
    This is called exponential map or Rodrigues' formula.
    Parameters
    ----------
    a : array-like, shape (4,)
        Axis of rotation and rotation angle: (x, y, z, angle)
    Returns
    -------
    R : array-like, shape (3, 3)
        Rotation matrix
    """
    ux, uy, uz, theta = a
    c = np.cos(theta)
    s = np.sin(theta)
    ci = 1.0 - c
    R = np.array([[ci * ux * ux + c,
                   ci * ux * uy - uz * s,
                   ci * ux * uz + uy * s],
                  [ci * uy * ux + uz * s,
                   ci * uy * uy + c,
                   ci * uy * uz - ux * s],
                  [ci * uz * ux - uy * s,
                   ci * uz * uy + ux * s,
                   ci * uz * uz + c],
                  ])

    # This is equivalent to
    # R = (np.eye(3) * np.cos(theta) +
    #      (1.0 - np.cos(theta)) * a[:3, np.newaxis].dot(a[np.newaxis, :3]) +
    #      cross_product_matrix(a[:3]) * np.sin(theta))

    return R


def resample(image, transform, interpolator):
    """
    This function resamples (updates) an image using a specified transform
    :param image: The sitk image we are trying to transform
    :param transform: An sitk transform (ex. resizing, rotation, etc.
    :return: The transformed sitk image
    """
    reference_image = image
    default_value = 0
    return sitk.Resample(image, reference_image, transform,
                         interpolator, default_value)


def get_center(img):
    """
    This function returns the physical center point of a 3d sitk image
    :param img: The sitk image we are trying to find the center of
    :return: The physical center point of the image
    """
    width, height, depth = img.GetSize()
    return img.TransformIndexToPhysicalPoint((int(np.ceil(width/2)),
                                              int(np.ceil(height/2)),
                                              int(np.ceil(depth/2))))


def rotation3d(image, theta_z, interpolator, show=False):
    """
    This function rotates an image across each of the x, y, z axes by theta_x, theta_y, and theta_z degrees
    respectively
    :param image: An sitk MRI image
    :param theta_x: The amount of degrees the user wants the image rotated around the x axis
    :param theta_y: The amount of degrees the user wants the image rotated around the y axis
    :param theta_z: The amount of degrees the user wants the image rotated around the z axis
    :param interpolator: Type of interpolation (sitk.sitkNearestNeighbor, sitk.sitkLinear, sitk.sitkBSpline)
    :param show: Boolean, whether or not the user wants to see the result of the rotation
    :return: The rotated image
    """
    theta_z = np.deg2rad(theta_z)
    euler_transform = sitk.Euler3DTransform()
    image_center = get_center(image)
    euler_transform.SetCenter(image_center)

    direction = image.GetDirection()
    axis_angle = (direction[2], direction[5], direction[8], theta_z)
    np_rot_mat = matrix_from_axis_angle(axis_angle)
    euler_transform.SetMatrix(np_rot_mat.flatten().tolist())
    resampled_image = resample(image, euler_transform, interpolator)
    if show:
        print(euler_transform.GetMatrix())
        slice_num = int(input("Enter the index of the slice you would like to see"))
        plt.imshow(sitk.GetArrayFromImage(resampled_image)[slice_num])
        plt.show()
    return resampled_image
    
# +------------------------------------------------------------------------------+
def IntensityProjection(input_image, output_image, mask_image, mask_image_2d, yaml_file, selection, projection_type, rotation, overwrite):

    # Check if output image exists and should overwrite
    if os.path.isfile(output_image) and not overwrite:
        result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(output_image))
        if result.lower() not in ['y', 'yes']:
            ogo.message('Not overwriting. Exiting...')
            os.sys.exit()

    # Check if output report exists and should overwrite
    if not yaml_file is None:
        if os.path.isfile(yaml_file) and not overwrite:
            result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(yaml_file))
            if result.lower() not in ['y', 'yes']:
                ogo.message('Not overwriting. Exiting...')
                os.sys.exit()
    
    # Read input
    if not os.path.isfile(input_image):
        os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(input_image))

    if not (input_image.lower().endswith('.nii') or input_image.lower().endswith('.nii.gz')):
        os.sys.exit('[ERROR] Input must be type NIfTI file: \"{}\"'.format(input_image))
    
    ogo.message('Reading input image: ')
    ogo.message('  {}'.format(input_image))
    ct = sitk.ReadImage(input_image, sitk.sitkInt16) # Read in raw CT data as 16-bit signed integer
    ogo.aix_nifti(input_image,ct)
    ogo.message('[WARNING] By forcing input image to Int16 some round-off error occurs.')
    ogo.message('          We could read in the CT without casting and then use')
    ogo.message('          ct.GetPixelID() to use the correct image type throughout.')
    
    # Set the type of projection
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
    
    # Remove negative values (e.g. air) from CT image
    ogo.message('Removing negative values in input CT image.')
    ct_no_zeros = sitk.Cast(ct>0, sitk.sitkInt16)
    ct = ct * ct_no_zeros
    
    interpolator = sitk.sitkLinear # options are sitk.sitkBSpline, sitk.sitkLinear
    
    # Transform the image and then create the projection
    ogo.message('Applying transform of Z-axis rotation of {:.2f} degrees.'.format(rotation)) 
    ct_transformed = rotation3d(ct, rotation, interpolator, False) 

    ogo.message('Performing intensity projection.')
    ct_projection = projectionFilt.Execute(ct_transformed)
    
    voxel_volume = ct.GetSpacing()[0] * ct.GetSpacing()[1] * ct.GetSpacing()[2] # mm3
    projected_voxel_area = ct_projection.GetSpacing()[1] * ct_projection.GetSpacing()[2] # area in mm of projected voxels
    projected_voxel_height = ct_projection.GetSpacing()[0]
    projected_voxel_volume = ct_projection.GetSpacing()[0] * ct_projection.GetSpacing()[1] * ct_projection.GetSpacing()[2]
    
    # Collection information into dictionary for YAML file
    info_dict = {}
    info_dict['runtime']={}
    info_dict['runtime']['time']=datetime.now().strftime("%H:%M:%S")
    info_dict['runtime']['script']=os.path.splitext(os.path.basename(sys.argv[0]))[0]
    info_dict['runtime']['version']=script_version
    info_dict['input_parameters']={}
    info_dict['input_parameters']['input_image']=input_image
    info_dict['input_parameters']['mask_image']=mask_image
    info_dict['input_parameters']['output_image']=output_image
    info_dict['input_parameters']['mask_image_2d']=mask_image_2d
    info_dict['input_parameters']['yaml_file']=yaml_file
    info_dict['input_parameters']['selection']=list(selection)
    info_dict['input_parameters']['projection_type']=projection_type
    info_dict['input_parameters']['rotation']=rotation
    info_dict['input_parameters']['overwrite']=overwrite
    info_dict['input_image']={}
    info_dict['input_image']['dim']=list(ct.GetSize())
    info_dict['input_image']['el_size_mm']=list(ct.GetSpacing())
    info_dict['input_image']['pixel_type']=ct.GetPixelIDTypeAsString()
    info_dict['output_image']={}
    info_dict['output_image']['dim']=list(ct_projection.GetSize())
    info_dict['output_image']['el_size_mm']=list(ct_projection.GetSpacing())
    info_dict['output_image']['pixel_type']=ct_projection.GetPixelIDTypeAsString()
    info_dict['output_image']['interpolator']=interpolator
            
    # If a mask is provided perform an analysis
    if mask_image:
     
        # Read mask
        if not os.path.isfile(mask_image):
            os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(mask_image))
        if not (mask_image.lower().endswith('.nii') or mask_image.lower().endswith('.nii.gz')):
            os.sys.exit('[ERROR] Mask must be type NIfTI file: \"{}\"'.format(mask_image))
        ogo.message('Reading mask: ')
        ogo.message('  {}'.format(mask_image))
        ct_mask = sitk.ReadImage(mask_image, sitk.sitkUInt8)
        ogo.aix_nifti(mask_image,ct_mask)
        
        labels = get_labels(ct_mask)
        ogo.message('Mask image contains the following labels:')
        ogo.message('  [' + ', '.join('{:d}'.format(i) for i in labels) + ']')
        
        if not selection:
            selection = labels
        else:
            for label in selection:
                if label not in labels:
                    os.sys.exit('[ERROR] Label {} is not in input mask.'.format(label))
        ogo.message('Projections will be calculated for the following labels:')
        ogo.message('  [' + ', '.join('{:d}'.format(i) for i in selection) + ']')
        
        # Check if mask_image_2d image exists and should overwrite
        if not mask_image_2d is None:
            if os.path.isfile(mask_image_2d) and not overwrite:
                result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(mask_image_2d))
                if result.lower() not in ['y', 'yes']:
                    ogo.message('Not overwriting. Exiting...')
                    os.sys.exit()
            if not (mask_image_2d.lower().endswith('.nii') or mask_image_2d.lower().endswith('.nii.gz')):
                os.sys.exit('[ERROR] Mask 2D image must be type NIfTI file: \"{}\"'.format(mask_image_2d))
            
        mask_transformed = rotation3d(ct_mask, rotation, sitk.sitkNearestNeighbor, False)
        
        maskProjectionFilt = sitk.BinaryProjectionImageFilter() # Special projection for labels
        
        stats = sitk.LabelIntensityStatisticsImageFilter()
        
        info_dict['analysis']={}
        info_dict['analysis']['labels_available']=list(labels)
        info_dict['analysis']['labels_used']=list(selection)
                
        seg_2d = maskProjectionFilt.Execute(mask_transformed)
        seg_2d = seg_2d<0 # reset to all zeros
        seg_3d = mask_transformed<0 # reset to all zeros
        
        for label in selection:
            try:
                desc = lb.labels_dict[label]['LABEL']
            except KeyError:
                desc = 'unknown label'
            ogo.message('  processing label {:>2d} ({})'.format(label,desc))
            
            # Area measures
            mask = mask_transformed==label
            mask_projection = maskProjectionFilt.Execute(mask)
            mask_projection = label * mask_projection
            stats.Execute(mask_projection,ct_projection)
            area_2d = stats.GetPhysicalSize(label) / projected_voxel_height / 100.0 # cm2 --> from mm3
            bmc_2d = stats.GetSum(label) / 1000.0 / 1000.0 * voxel_volume # g --> from mg/cm3
            bmd_2d = bmc_2d / area_2d # g/cm2
                        
            # Volume measures
            stats.Execute(label*mask, ct_transformed)
            volume_3d = stats.GetPhysicalSize(label) / 1000.0 # cm3 --> from mm3
            bmd_3d = stats.GetMean(label) # mg/cm3
            bmc_3d = bmd_3d * volume_3d / 1000.0 # g --> from mg
            
            # Capture total 2D projection 
            bin_seg_2d = seg_2d>0
            bin_part_2d = mask_projection>0
            overlap_2d = (bin_seg_2d + bin_part_2d)==2
            mask_2d = 1 - overlap_2d
            this_label_2d = sitk.Mask(mask_projection, mask_2d)
            seg_2d = seg_2d + label*(this_label_2d>0)
            
            # Capture total 3D projection 
            bin_seg_3d = seg_3d>0
            bin_part_3d = mask>0
            overlap_3d = (bin_seg_3d + bin_part_3d)==2
            mask_3d = 1 - overlap_3d
            this_label_3d = sitk.Mask(mask, mask_3d)
            seg_3d = seg_3d + label*(this_label_3d>0)
            
            info_dict['analysis'][label]={}
            info_dict['analysis'][label]['desc'] = desc
            info_dict['analysis'][label]['label'] = label
            info_dict['analysis'][label]['bmd_2d'] = bmd_2d
            info_dict['analysis'][label]['bmd_3d'] = bmd_3d
            info_dict['analysis'][label]['bmc_3d'] = bmc_3d
            info_dict['analysis'][label]['area_2d'] = area_2d
            info_dict['analysis'][label]['volume_3d'] = volume_3d
                     
        # Write projection of all bones
        if not mask_image_2d is None:
            #temp_name = temp_base_name + '_2D_' + '{:.0f}'.format(rotation) + '.nii.gz'
            sitk.WriteImage(seg_2d, mask_image_2d)
        
        # Integral calculations
        seg_2d = seg_2d>0 # Set all labels to 1
        seg_3d = seg_3d>0 # Set all labels to 1
        
        label=1
        stats.Execute(seg_2d,ct_projection)
        area_2d = stats.GetPhysicalSize(label) / projected_voxel_height / 100.0 # cm2 --> from mm3
        bmc_2d = stats.GetSum(label) / 1000.0 / 1000.0 * voxel_volume # g --> from mg/cm3
        bmd_2d = bmc_2d / area_2d # g/cm2
        
        stats.Execute(seg_3d, ct_transformed)
        volume_3d = stats.GetPhysicalSize(label) / 1000.0 # cm3 --> from mm3
        bmd_3d = stats.GetMean(label) # mg/cm3
        bmc_3d = bmd_3d * volume_3d / 1000.0 # g --> from mg
                
        # Finalize report on mask
        info_dict['analysis']['integral']={}
        info_dict['analysis']['integral']['desc'] = 'Integral'
        info_dict['analysis']['integral']['label'] = 0
        info_dict['analysis']['integral']['bmd_2d'] = bmd_2d
        info_dict['analysis']['integral']['bmd_3d'] = bmd_3d
        info_dict['analysis']['integral']['bmc_3d'] = bmc_3d
        info_dict['analysis']['integral']['area_2d'] = area_2d
        info_dict['analysis']['integral']['volume_3d'] = volume_3d
             
    else:
     ogo.message('No mask for quantitative analysis provided.')
    
    if output_image:
        ogo.message('Writing projection file:')
        ogo.message('  {}'.format(output_image))
        sitk.WriteImage(ct_projection, output_image)
    else:
        ogo.message('No output projection image written to file.')            
    
    info_dict['message']={}
    info_dict['message']['warning']='The quantitative outputs have not been fully checked. There may be errors in density, mass and volume/area.'
    
    if yaml_file:
        ogo.message('Saving report to file:')
        ogo.message('  {}'.format(yaml_file))
        with open(yaml_file, 'w') as file:
            documents = yaml.dump(info_dict, file, sort_keys=False)
        
    ogo.message('Done.')
    
    
def main():
    # Setup description
    description = '''
Generate an intensity projection of a 3D dataset. User can select the rotation
in degrees around the Z-axis. A typical A-P scan uses a rotation angle of 90 
degrees.

If a label mask is supplied that identifies the bones, then aBMD and vBMD and
associated measures of area, volume and bone mass are calculated. It is 
optional to output the projected 2D mask.

Negative values in the CT scan (e.g. air) are removed before calculating the
projected BMC and aBMD measurements.
'''

    epilog = '''
Example call: 
     
ogoIntensityProjection image.nii.gz image_2d.nii.gz --mask_image mask.nii.gz --rotation 45

'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoIntensityProjection",
        description=description,
        epilog=epilog
    )

    parser.add_argument('input_image',  metavar='NIfTI', 
                                        help='Input image file (*.nii, *.nii.gz)')
    parser.add_argument('output_image', metavar='NIfTI', 
                                        help='Output image file (*.nii, *.nii.gz)')
    parser.add_argument('--mask_image', default=None, metavar='NIfTI', 
                                        help='Bone labels (*.nii.gz, default: %(default)s)')
    parser.add_argument('--mask_image_2d', default=None, metavar='NIfTI', 
                                        help='Print the 2d projected mask (*.nii.gz, default: %(default)s)')
    parser.add_argument('--yaml_file', default=None, metavar='FILE', 
                                        help='Write report to file (*.yaml, default: %(default)s)')
    parser.add_argument('--selection', type=int, nargs='*', default=[], metavar='LABEL', 
                                        help='List of mask labels to use (default: all)')
    parser.add_argument('--projection_type', default='sum', 
                                        choices=['median', 'mean', 'maximum', 'minimum', 'stdev', 'sum'],
                                        help='Select type of projection (default: %(default)s)')
    parser.add_argument('--rotation', type=float, default=90.0, metavar='ANGLE',
                                        help='Rotation in degrees about Z axis (default: %(default)s deg)')
    parser.add_argument('--overwrite', action='store_true', 
                                        help='Overwrite output file without asking')

    # Parse and display
    args = parser.parse_args()
    print(echo_arguments('IntensityProjection', vars(args)))

    # Run program
    IntensityProjection(**vars(args))


if __name__ == '__main__':
    main()
