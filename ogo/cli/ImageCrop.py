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
import math
import numpy as np
from ogo.util.echo_arguments import echo_arguments
import ogo.util.Helper as ogo
import ogo.dat.OgoMasterLabels as lb
from ogo.util.cli_tools import (
    check_overwrite,
    check_image_dimensions_match,
    validate_integer_range,
    validate_input_file,
    validate_output_file,
    print_image_info
)

def ImageCrop(input_file, output_file, input_label_file, output_label_file, voi, labels, slab, offset, overwrite, quiet):
    
    if not quiet:
        ogo.message('Output image will be:')
        ogo.message('      "{}"'.format(output_file))

    if output_label_file is not None and not quiet:
        ogo.message('Output label will be:')
        ogo.message('      "{}"'.format(output_label_file))

    # Check if outputs exist and should overwrite
    # Collect all output files that need to be checked
    output_files = [output_file]
    if output_label_file is not None:
        output_files.append(output_label_file)
    
    # Check all output files together - only prompt once if needed
    if not overwrite:
        existing_files = [f for f in output_files if os.path.isfile(f)]
        if existing_files:
            files_list = '\n      '.join(['"{}"'.format(f) for f in existing_files])
            result = input(f'The following file(s) already exist:\n      {files_list}\n    Overwrite? [y/n]: ')
            if result.lower() not in ['y', 'yes']:
                print('Not overwriting. Exiting...')
                return 1

    # Read input
    if not quiet:
        ogo.message('Reading input CT image to be cropped:')
        ogo.message('      "{}"'.format(input_file))
    ct = sitk.ReadImage(input_file)
    if not quiet:
        print_image_info(input_file, ct)
    
    # If using labels to define bounds
    if input_label_file is not None:
        if not quiet:
            ogo.message('Reading label image')
            ogo.message('      "{}"'.format(input_label_file))
        if os.path.isfile(input_label_file):
            ct_labels = sitk.ReadImage(input_label_file, sitk.sitkUInt8)
            if not quiet:
                print_image_info(input_label_file, ct_labels)
            
        # check dimensions of label and ct images are the same
        if not check_image_dimensions_match(ct, ct_labels, "input_file", "input_label_file"):
            return 1

        if not labels:
            ogo.message('[ERROR] Labels must be defined to crop using labels.')
            return 1
            
        # create mask of labels
        ct_mask = ct_labels==-1 # we want these to be all zeros
        filt = sitk.LabelShapeStatisticsImageFilter()
        filt.Execute(ct_labels)
        
        for label in filt.GetLabels():
            try:
                desc = lb.labels_dict[label]['LABEL']
            except KeyError:
                desc = 'unknown label'
            
            bounding_box = filt.GetBoundingBox(label)
            if not quiet:
                ogo.message('  label {:2d} ({:s}):'.format(label,desc))
            
            if label in labels:
                if not quiet:
                    ogo.message('    [' + ', '.join('{:3d}'.format(i) for i in bounding_box) + '] - used to define VOI')
                bin_part = ct_labels==label
                mask = 1 - bin_part
                ct_mask = sitk.Mask(ct_mask, mask)
                ct_mask = ct_mask + bin_part
            
            else:
                if not quiet:
                    ogo.message('    [' + ', '.join('{:3d}'.format(i) for i in bounding_box) + ']')
                
        filt.Execute(ct_mask)
        voi = list(filt.GetBoundingBox(1))  # Convert tuple to list
        if not quiet:
            ogo.message('')
            ogo.message('  VOI based on labels:')
            ogo.message('    [' + ', '.join('{:3d}'.format(i) for i in voi) + ']')
        
        # Apply offset if specified
        if offset is not None:
            dim = ct.GetSize()
            # Calculate new VOI bounds with offset, clamped to image dimensions
            requested_start_x = voi[0] - offset[0]
            requested_end_x = voi[0] + voi[3] + offset[0]
            requested_start_y = voi[1] - offset[1]
            requested_end_y = voi[1] + voi[4] + offset[1]
            requested_start_z = voi[2] - offset[2]
            requested_end_z = voi[2] + voi[5] + offset[2]
            
            new_start_x = max(0, requested_start_x)
            new_end_x = min(dim[0], requested_end_x)
            new_start_y = max(0, requested_start_y)
            new_end_y = min(dim[1], requested_end_y)
            new_start_z = max(0, requested_start_z)
            new_end_z = min(dim[2], requested_end_z)
            
            # Check if clamping occurred and warn the user
            clamped_dims = []
            if requested_start_x < 0 or requested_end_x > dim[0]:
                clamped_dims.append('X')
            if requested_start_y < 0 or requested_end_y > dim[1]:
                clamped_dims.append('Y')
            if requested_start_z < 0 or requested_end_z > dim[2]:
                clamped_dims.append('Z')
            
            if clamped_dims and not quiet:
                ogo.message('[WARNING] Offset expansion exceeded image boundaries in dimension(s): {}'.format(', '.join(clamped_dims)))
                ogo.message('[WARNING] VOI has been clamped to fit within image bounds [0-{}, 0-{}, 0-{}]'.format(dim[0]-1, dim[1]-1, dim[2]-1))
            
            # Update VOI with new values
            voi[0] = new_start_x
            voi[1] = new_start_y
            voi[2] = new_start_z
            voi[3] = new_end_x - new_start_x
            voi[4] = new_end_y - new_start_y
            voi[5] = new_end_z - new_start_z
            
            if not quiet:
                ogo.message('')
                ogo.message('  VOI after applying offset [{}, {}, {}]:'.format(*offset))
                ogo.message('    [' + ', '.join('{:3d}'.format(i) for i in voi) + ']')

    # Reset VOI if slab
    dim = ct.GetSize()
    if slab:
        if not quiet:
            ogo.message('Using Z bounds only.')
        voi = [0, 0, voi[2], dim[0], dim[1], voi[5]]
        if not quiet:
            ogo.message('')
            ogo.message('  Final VOI:')
            ogo.message('    [' + ', '.join('{:3d}'.format(i) for i in voi) + ']')
    else:
        if not quiet:
            ogo.message('')
            ogo.message('  Final VOI:')
            ogo.message('    [' + ', '.join('{:3d}'.format(i) for i in voi) + ']')
        
    # Take the subvolume
    if voi[0]<0 or voi[1]<0 or voi[2]<0 or voi[3]>dim[0] or voi[4]>dim[1] or voi[5]>dim[2]:
        ogo.message('[ERROR] Defined VOI is out of range: {} {} {} {} {} {}'.format(*voi))
        return 1
    
    ct_out = ct[voi[0]:(voi[0]+voi[3]), voi[1]:(voi[1]+voi[4]), voi[2]:(voi[2]+voi[5])]
    if not quiet:
        ogo.message('Writing output image to file {}'.format(output_file))
    sitk.WriteImage(ct_out, output_file)
    if not quiet:
        print_image_info(output_file, ct_out)
    
    if input_label_file is not None and output_label_file is not None:
        ct_label_out = ct_labels[voi[0]:(voi[0]+voi[3]), voi[1]:(voi[1]+voi[4]), voi[2]:(voi[2]+voi[5])]
        if not quiet:
            ogo.message('Writing output labels to file {}'.format(output_label_file))
        sitk.WriteImage(ct_label_out, output_label_file)
        if not quiet:
            print_image_info(output_label_file, ct_label_out)
    
    if not quiet:
        ogo.message('Done ogoImageCrop!')
    
    return 0


def main():
    # Setup description
    description = '''
Utility to crop a NIfTI file to a subvolume.

There are two main methods of cropping. One is to define the volume of interest (VOI)
and the other is to use the labels in a mask of the image to be cropped.

Define the VOI using six integers:origX, origY, origZ, dimX, dimY, dimZ
For example: --voi 50 100 150 200 250 300
This crops from voxel (50,100,150) with dimensions (200,250,300).

Define an input file with labels (usually from ML or similar) and specify which labels
to use as a basis for cropping. It will automatically find the VOI that includes all
the specified labels.

For example: --input_label_file labels.nii.gz --labels 5 10 15
This finds the smallest VOI containing labels 5, 10, and 15, then crops.

Optional flags:

--slab:   When using label-based cropping, this option keeps the full X and Y
          dimensions and only crops in the Z direction based on the label bounds.
          Useful for extracting vertebral slabs or similar structures.

--offset: When using label-based cropping, this option adds a symmetric offset to the
          VOI defined by the labels. For example, --offset 10 10 5 adds 10 voxels on 
          each side in X (20 total), 10 on each side in Y (20 total), and 5 on each 
          side in Z (10 total).

Outputs are the cropped image and optionally the cropped label file if an input label 
file was used. The output label file will contain the same labels as the input label 
file, but only within the cropped region.          
'''
    epilog = '''
Example calls: 
ogoImageCrop input.nii.gz output.nii.gz --voi 0 512 0 512 137 448 
ogoImageCrop input.nii.gz output.nii.gz --input_label_file labels.nii.gz -l 8 9 10 --slab
ogoImageCrop input.nii.gz output_cropped.nii.gz --input_label_file labels.nii.gz -l 5 10 15 --offset 10 10 5 -ow -q
'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoImageCrop",
        description=description,
        epilog=epilog
    )
    parser.add_argument('input_file', type=validate_input_file(['.nii', '.nii.gz']), help='Input image file (*.nii, *.nii.gz)')
    parser.add_argument('output_file', type=validate_output_file(['.nii', '.nii.gz']), help='Output CT image file (*.nii, *.nii.gz)')
    parser.add_argument('--input_label_file', type=validate_input_file(['.nii', '.nii.gz']), metavar='FILE', default=None, help='Label image used to define crop (*.nii, *.nii.gz)')
    parser.add_argument('--output_label_file', type=validate_output_file(['.nii', '.nii.gz']), metavar='FILE', default=None, help='Output label file (*.nii, *.nii.gz)')
    parser.add_argument('--voi', type=int, nargs=6, default=[0,1,0,1,0,1], metavar='0', help='[origX, origY, origZ, dimX, dimY, dimZ')
    parser.add_argument('--labels', '-l', type=validate_integer_range(0, 255), nargs='*', default=[], metavar='LABEL', help='List of labels to define crop region (default: %(default)s)')
    parser.add_argument('--offset', type=int, nargs=3, default=None, metavar=('X', 'Y', 'Z'), help='Symmetric offset to expand VOI in label-based cropping [offsetX, offsetY, offsetZ]')
    parser.add_argument('--slab', action='store_true', help='Crop a Z slab (overrides X and Y VOI)')
    parser.add_argument('--overwrite', '-ow', action='store_true', help='Overwrite output without asking')
    parser.add_argument('--quiet', '-q', action='store_true', help='Suppress informational output to stdout')

    # Parse and display
    args = parser.parse_args()
    if not args.quiet:
        ogo.message(echo_arguments('ImageCrop', vars(args)))

    # Run program
    result = ImageCrop(**vars(args))
    return result


if __name__ == '__main__':
    main()
