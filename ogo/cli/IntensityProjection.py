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
def IntensityProjection(input_image, mask_image, output_image, selection, projection_type, rotation, print_projected_masks, overwrite):

    # Check if output exists and should overwrite
    if output_image:
        if os.path.isfile(output_image) and not overwrite:
            result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(output_image))
            if result.lower() not in ['y', 'yes']:
                ogo.message('Not overwriting. Exiting...')
                os.sys.exit()

    # Read input
    if not os.path.isfile(input_image):
        os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(input_image))

    if not (input_image.lower().endswith('.nii') or input_image.lower().endswith('.nii.gz')):
        os.sys.exit('[ERROR] Input must be type NIFTI file: \"{}\"'.format(input_image))
    
    ogo.message('Reading image: ')
    ogo.message('\"{}\"'.format(input_image))
    ct = sitk.ReadImage(input_image, sitk.sitkInt16) # Read in raw CT data as 16-bit signed integer
    
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
    
    # Remove negative values (e.g. air) from CT image
    ogo.message('Removing negative values in input CT image.')
    ct_no_zeros = sitk.Cast(ct>0, sitk.sitkInt16)
    ct = ct * ct_no_zeros
    
    # Transform the image and then create the projection
    ogo.message('Applying transform of Z-axis rotation of {:.2f} degrees.'.format(rotation)) 
    ct_transformed = sitk.Resample(ct, transform, sitk.sitkBSpline) # sitk.sitkBSpline, sitk.sitkNearestNeighbor, sitk.sitkLinear
    ogo.message('Performing intensity projection.')
    ct_projection = projectionFilt.Execute(ct_transformed)
    
    voxel_volume = ct.GetSpacing()[0] * ct.GetSpacing()[1] * ct.GetSpacing()[2] # mm3
    projected_voxel_area = ct_projection.GetSpacing()[1] * ct_projection.GetSpacing()[2] # area in mm of projected voxels
    projected_voxel_height = ct_projection.GetSpacing()[0]
    projected_voxel_volume = ct_projection.GetSpacing()[0] * ct_projection.GetSpacing()[1] * ct_projection.GetSpacing()[2]
    
    # Start the report
    report = ''
    report += '  {:>20s}\n'.format('_______________________________________________________________________Report')
    report += '  {:>20s} {:s}\n'.format('input image',input_image)
    report += '  {:>20s} {:s}\n'.format('output image',output_image if output_image else 'None')
    report += '  {:>20s} {:s}\n'.format('mask image',mask_image if mask_image else 'None')
    report += '\n'
    report += '  {:>20s}\n'.format('input image')
    report += '  {:>20s} '.format('dim:') + ' '.join('{:12d}'.format(i) for i in ct.GetSize()) + '\n'
    report += '  {:>20s} '.format('el_size_mm:') + ' '.join('{:12.4f}'.format(i) for i in ct.GetSpacing()) + '\n'
    report += '  {:>20s} {:12.4f} [mm3]\n'.format('voxel volume:',voxel_volume)
    report += '  {:>20s} {:>36s}\n'.format('pixel type:',ct.GetPixelIDTypeAsString())
    
    if output_image:
        report += '\n'
        report += '  {:>20s}\n'.format('output image')
        report += '  {:>20s} '.format('dim:') + ' '.join('{:12d}'.format(i) for i in ct_projection.GetSize()) + '\n'
        report += '  {:>20s} '.format('el_size_mm:') + ' '.join('{:12.4f}'.format(i) for i in ct_projection.GetSpacing()) + '\n'
        report += '  {:>20s} {:12.4f} [mm2]\n'.format('voxel area:',projected_voxel_area)
        report += '  {:>20s} {:12.4f} [mm3]\n'.format('voxel volume:',projected_voxel_volume)
    report += '\n'
    report += '  {:>20s} {:s}\n'.format('projection type:',str(projection_type))
    report += '  {:>20s} {:.1f} [deg]\n'.format('Z-axis rotation:',rotation)
    
    # If a mask is provided perform an analysis
    if mask_image:
        
        # Read mask
        if not os.path.isfile(mask_image):
            os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(mask_image))
        if not (mask_image.lower().endswith('.nii') or mask_image.lower().endswith('.nii.gz')):
            os.sys.exit('[ERROR] Mask must be type NIFTI file: \"{}\"'.format(mask_image))
        ogo.message('Reading mask: ')
        ogo.message('\"{}\"'.format(mask_image))
        ct_mask = sitk.ReadImage(mask_image, sitk.sitkUInt8)
        
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
        
        # Establish base name for temporary mask files
        temp_base_name, ext = os.path.splitext(mask_image)
        if 'gz' in ext:
            temp_base_name = os.path.splitext(temp_base_name)[0]  # Manages files with double extension
            ext = '.nii' + ext
        
        mask_transformed = sitk.Resample(ct_mask, transform, sitk.sitkNearestNeighbor)
        
        maskProjectionFilt = sitk.BinaryProjectionImageFilter() # Special projection for labels
        
        stats = sitk.LabelIntensityStatisticsImageFilter()
        
        # Add to the report if using a mask
        report += '\n'
        report += '  {:>20s}'.format('labels available: ') + ','.join('{:d}'.format(i) for i in labels) + '\n'
        report += '  {:>20s}'.format('labels used: ') + ','.join('{:d}'.format(i) for i in selection) + '\n'
        report += '\n'
        report += '  {:>20s} {:>10s} {:>10s} {:>10s} {:>10s} {:>10s}\n'.format('ROI','aBMD','vBMD','MASS','AREA','VOLUME')
        report += '  {:>20s} {:>10s} {:>10s} {:>10s} {:>10s} {:>10s}\n'.format('--','[g/cm2]','[mg/cm3]','[grams]','[cm2]','[cm3]')
        
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
            
            report += '  {:>15s}{:>5s} {:10.3f} {:10.1f} {:10.2f} {:10.2f} {:10.1f}\n'.format('('+desc+')',str(label),bmd_2d,bmd_3d,bmc_3d,area_2d,volume_3d)
            
            if print_projected_masks:
                temp_name = temp_base_name + '_TEMP_2D_' + '{:.0f}_'.format(rotation) + desc.replace(' ','') + '.nii.gz'
                sitk.WriteImage(mask_projection, temp_name)
            
        # Write projection of all bones
        if print_projected_masks:
            temp_name = temp_base_name + '_TEMP_2D_' + '{:.0f}_'.format(rotation) + 'all' + '.nii.gz'
            sitk.WriteImage(seg_2d, temp_name)
        
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
        report += '  {:>15s}{:>5s} {:10.3f} {:10.1f} {:10.2f} {:10.2f} {:10.1f}\n'.format('(Integral)','--',bmd_2d,bmd_3d,bmc_3d,area_2d,volume_3d)
        
    else:
        ogo.message('No mask for quantitative analysis provided.')
    
    if output_image:
        ogo.message('Writing projection file:')
        ogo.message('  {}'.format(output_image))
        sitk.WriteImage(ct_projection, output_image)
    else:
        ogo.message('No output projection image written to file.')            
    
    report += '\n'
    report += 'WARNING! The quantitative outputs have not been fully checked.\n'
    report += '         There may be errors in density, mass and volume/area.\n'
    
    print(report)
    ogo.message('Done.')
    
    
def main():
    # Setup description
    description = '''
Generate an intensity projection of a 3D dataset. User can select the rotation
in degrees around the Z-axis. A typical A-P scan uses a rotation angle of 90 
degrees.

If a label mask is supplied that identifies the bones, then aBMD and vBMD and
associated measures of area, volume and bone mass are calculated.

Negative values in the CT scan (e.g. air) are removed before calculating the
projected BMC and aBMD measurements.
'''

    epilog = '''
Example call: 
     
ogoIntensityProjection image.nii.gz --mask_image mask.nii.gz --rotation 92

'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoIntensityProjection",
        description=description,
        epilog=epilog
    )

    parser.add_argument('input_image', 
                        help='Input image file (*.nii, *.nii.gz)')
    parser.add_argument('--output_image', default=None, metavar='NIfTI', 
                        help='Output image file (*.nii.gz, default: %(default)s)')
    parser.add_argument('--mask_image', default=None, metavar='NIfTI', 
                        help='Bone labels (*.nii.gz, default: %(default)s)')
    parser.add_argument('--selection', type=int, nargs='*', default=[], metavar='LABEL', 
                        help='List of mask labels to use (default: all)')
    parser.add_argument('--projection_type', default='sum', 
                        choices=['median', 'mean', 'maximum', 'minimum', 'stdev', 'sum'],
                        help='Select type of projection (default: %(default)s)')
    parser.add_argument('--rotation', type=float, default=90.0, metavar='ANGLE',
                        help='Rotation in degrees about Z axis (default: %(default)s deg)')
    parser.add_argument('--print_projected_masks', action='store_true', 
                        help='Print the 2d projected masks as TEMP files (*.nii.gz)')
    parser.add_argument('--overwrite', action='store_true', 
                        help='Overwrite output file without asking')

    # Parse and display
    args = parser.parse_args()
    print(echo_arguments('IntensityProjection', vars(args)))

    # Run program
    IntensityProjection(**vars(args))


if __name__ == '__main__':
    main()
