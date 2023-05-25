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

def printImageInfo(im):
    dim = im.GetSize()
    spacing = im.GetSpacing()
    origin = im.GetOrigin()
    phys_dim = [x * y for x, y in zip(dim, spacing)]
    position = [math.floor(x / y) for x, y in zip(origin, spacing)]
    
    guard = '!-------------------------------------------------------------------------------'
    print(guard)
    print('!> dim                            {:>8}  {:>8}  {:>8}'.format(*dim))
    print('!> off                            {:>8}  {:>8}  {:>8}'.format('-', '-', '-'))
    print('!> pos                            {:>8}  {:>8}  {:>8}'.format(*position))
    print('!> element size in mm             {:>8.4f}  {:>8.4f}  {:>8.4f}'.format(*spacing))
    print('!> phys dim in mm                 {:>8.4f}  {:>8.4f}  {:>8.4f}'.format(*phys_dim))
    print(guard)
    
def ImageCrop(input_image, output_image, output_label, label_image, voi, labels, slab, overwrite):
    
    # Define output_image if not set explicitly
    if output_image is None:
        basename = os.path.basename(input_image)
        name, ext = os.path.splitext(input_image)
        if 'gz' in ext:
            name = os.path.splitext(name)[0]  # Manages files with double extension
            ext = '.nii' + ext
        output_image = name + '_CROP.nii.gz'
    ogo.message('Output image will be:')
    ogo.message('      \"{}\"'.format(output_image))

    # Define output_label if not set explicitly
    if output_label is None and label_image is not None:
        basename = os.path.basename(label_image)
        name, ext = os.path.splitext(label_image)
        if 'gz' in ext:
            name = os.path.splitext(name)[0]  # Manages files with double extension
            ext = '.nii' + ext
        output_label = name + '_CROP.nii.gz'
    if label_image is not None:
        ogo.message('Output label will be:')
        ogo.message('      \"{}\"'.format(output_label))

    # Check if output exists and should overwrite
    if os.path.isfile(output_image) and not overwrite:
        result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(output_image))
        if result.lower() not in ['y', 'yes']:
            ogo.message('Not overwriting. Exiting...')
            os.sys.exit()

    # Read input
    if not os.path.isfile(input_image):
        os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(input_image))

    ogo.message('Reading input CT image to be cropped:')
    ogo.message('      \"{}\"'.format(input_image))
    ct = sitk.ReadImage(input_image)
    printImageInfo(ct)
    
    # If using labels to define bounds
    if label_image is not None:
        ogo.message('Reading label image')
        ogo.message('      \"{}\"'.format(label_image))
        if os.path.isfile(label_image):
            ct_labels = sitk.ReadImage(label_image, sitk.sitkUInt8)
            printImageInfo(ct_labels)
            
        # check dimensions of label and ct images are the same
        if not ((ct.GetSize()[0] == ct_labels.GetSize()[0]) and \
                (ct.GetSize()[1] == ct_labels.GetSize()[1]) and \
                (ct.GetSize()[2] == ct_labels.GetSize()[2])):
            ogo.message('[ERROR] Label image must be same dimensions as input image.')
            os.sys.exit()

        if not labels:
            ogo.message('Labels must be defined to crop using labels.')
            os.sys.exit()
            
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
            ogo.message('  label {:2d} ({:s}):'.format(label,desc))
            
            if label in labels:
                ogo.message('    [' + ', '.join('{:3d}'.format(i) for i in bounding_box) + '] - used to define VOI')
                bin_part = ct_labels==label
                mask = 1 - bin_part
                ct_mask = sitk.Mask(ct_mask, mask)
                ct_mask = ct_mask + bin_part
            
            else:
                ogo.message('    [' + ', '.join('{:3d}'.format(i) for i in bounding_box) + ']')
                
        filt.Execute(ct_mask)
        voi = filt.GetBoundingBox(1)
        ogo.message('')
        ogo.message('  VOI based on labels:')
        ogo.message('    [' + ', '.join('{:3d}'.format(i) for i in voi) + ']')

    # Reset VOI if slab
    dim = ct.GetSize()
    if slab:
        ogo.message('Using Z bounds only.')
        voi = [0, 0, voi[2], dim[0], dim[1], voi[5]]
        ogo.message('')
        ogo.message('  Final VOI:')
        ogo.message('    [' + ', '.join('{:3d}'.format(i) for i in voi) + ']')
    else:
        ogo.message('')
        ogo.message('  Final VOI:')
        ogo.message('    [' + ', '.join('{:3d}'.format(i) for i in voi) + ']')
        
    # Take the subvolume
    if voi[0]<0 or voi[1]<0 or voi[2]<0 or voi[3]>dim[0] or voi[4]>dim[1] or voi[5]>dim[2]:
        os.sys.exit('[ERROR] Defined VOI is out of range: {} {} {} {} {} {}'.format(*voi))
    
    ct_out = ct[voi[0]:(voi[0]+voi[3]), voi[1]:(voi[1]+voi[4]), voi[2]:(voi[2]+voi[5])]
    ogo.message('Writing output image to file {}'.format(output_image))
    sitk.WriteImage(ct_out, output_image)
    printImageInfo(ct_out)
    
    if label_image is not None:
        ct_label_out = ct_labels[voi[0]:(voi[0]+voi[3]), voi[1]:(voi[1]+voi[4]), voi[2]:(voi[2]+voi[5])]
        ogo.message('Writing output labels to file {}'.format(output_label))
        sitk.WriteImage(ct_label_out, output_label)
        printImageInfo(ct_label_out)
        
    ogo.message('Done ogoImageCrop!')


def main():
    # Setup description
    description = '''
Utility to crop a NIFTI file to a subvolume. 

Cropping can be performed based on the labels within an image. A list
of labels can be defined.
'''
    epilog = '''
Example calls: 
ogoImageCrop input.nii.gz output.nii.gz --voi 0 512 0 512 137 448 
ogoImageCrop input.nii.gz output.nii.gz --label_image labels.nii.gz --labels 8 9 10 --slab
'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoImageCrop",
        description=description,
        epilog=epilog
    )
    parser.add_argument('input_image', help='Input image file (*.nii, *.nii.gz)')
    parser.add_argument('--output_image', metavar='FILE', default=None, help='Output CT image file (*.nii, *.nii.gz)')
    parser.add_argument('--output_label', metavar='FILE', default=None, help='Output label file (*.nii, *.nii.gz)')
    parser.add_argument('--label_image', metavar='FILE', default=None, help='Label image used to define crop (*.nii, *.nii.gz)')
    parser.add_argument('--voi', type=int, nargs=6, default=[0,1,0,1,0,1], metavar='0', help='[origX, origY, origZ, dimX, dimY, dimZ')
    parser.add_argument('--labels', type=int, nargs='*', default=[8], metavar='LABEL', help='List of labels to define crop region (default: %(default)s)')
    parser.add_argument('--slab', action='store_true', help='Crop a Z slab (overrides X and Y VOI)')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite output without asking')

    # Parse and display
    args = parser.parse_args()
    print(echo_arguments('ImageCrop', vars(args)))

    # Run program
    ImageCrop(**vars(args))


if __name__ == '__main__':
    main()
