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
    print_image_info
)

def ImageCrop(input_file, output_file, output_label_file, input_label_file, voi, labels, slab, overwrite, quiet):
    
    # Define output_file if not set explicitly
    if output_file is None:
        basename = os.path.basename(input_file)
        name, ext = os.path.splitext(input_file)
        if 'gz' in ext:
            name = os.path.splitext(name)[0]  # Manages files with double extension
            ext = '.nii' + ext
        output_file = name + '_CROP.nii.gz'
    if not quiet:
        ogo.message('Output image will be:')
        ogo.message('      "{}"'.format(output_file))

    # Define output_label_file if not set explicitly
    if output_label_file is None and input_label_file is not None:
        basename = os.path.basename(input_label_file)
        name, ext = os.path.splitext(input_label_file)
        if 'gz' in ext:
            name = os.path.splitext(name)[0]  # Manages files with double extension
            ext = '.nii' + ext
        output_label_file = name + '_CROP.nii.gz'
    if input_label_file is not None and not quiet:
        ogo.message('Output label will be:')
        ogo.message('      "{}"'.format(output_label_file))

    # Check if output exists and should overwrite
    if not check_overwrite(output_file, overwrite):
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
        voi = filt.GetBoundingBox(1)
        if not quiet:
            ogo.message('')
            ogo.message('  VOI based on labels:')
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
    
    if input_label_file is not None:
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
Utility to crop a NIFTI file to a subvolume. 

Cropping can be performed based on the labels within an image. A list
of labels can be defined.
'''
    epilog = '''
Example calls: 
ogoImageCrop input.nii.gz --output_file output.nii.gz --voi 0 512 0 512 137 448 
ogoImageCrop input.nii.gz --output_file output.nii.gz --input_label_file labels.nii.gz -l 8 9 10 --slab
ogoImageCrop input.nii.gz --input_label_file labels.nii.gz -l 5 10 15 -ow -q
'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoImageCrop",
        description=description,
        epilog=epilog
    )
    parser.add_argument('input_file', type=validate_input_file(['.nii', '.nii.gz']), help='Input image file (*.nii, *.nii.gz)')
    parser.add_argument('--output_file', metavar='FILE', default=None, help='Output CT image file (*.nii, *.nii.gz)')
    parser.add_argument('--output_label_file', metavar='FILE', default=None, help='Output label file (*.nii, *.nii.gz)')
    parser.add_argument('--input_label_file', metavar='FILE', default=None, help='Label image used to define crop (*.nii, *.nii.gz)')
    parser.add_argument('--voi', type=int, nargs=6, default=[0,1,0,1,0,1], metavar='0', help='[origX, origY, origZ, dimX, dimY, dimZ')
    parser.add_argument('--labels', '-l', type=validate_integer_range(0, 255), nargs='*', default=[], metavar='LABEL', help='List of labels to define crop region (default: %(default)s)')
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
