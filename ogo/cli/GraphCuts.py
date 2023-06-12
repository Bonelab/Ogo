# /------------------------------------------------------------------------------+
# | 18-JAN-2023                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+

script_version = 1.0

# Imports
import argparse
import os
import sys
import subprocess
import vtk
import glob
import numpy as np
import SimpleITK as sitk
from datetime import date
from ogo.util.echo_arguments import echo_arguments
import ogo.util.Helper as ogo
import ogo.dat.OgoMasterLabels as lb
        
# +------------------------------------------------------------------------------+
def sheetness(input_image, sheet_image, path_binaries, enhance, \
              num_sigma, min_sigma, max_sigma, air_thres, metal_thres, trace_weight, \
              overwrite, func):

    # Check that the binaries are available to run the graphcuts algorithm
    sheetness_binary_path = ogo.find_executable('Sheetness2',path_binaries)
    if sheetness_binary_path is None:
        ogo.message('[ERROR] You must build the binaries for graph cuts for this to work.')
        ogo.message('        See instructions for installation here:')
        ogo.message('{}'.format(ogo.find_executable('OGO_GRAPHCUTS_INSTALL.txt',path_binaries)))
        os.sys.exit()

    # Check if input image exists
    if not os.path.isfile(input_image):
        ogo.message('[ERROR] Cannot find input file.')
        ogo.message('  {}'.format(input_image))
        os.sys.exit()
    if not (input_image.lower().endswith('.nii') or input_image.lower().endswith('.nii.gz')):
        ogo.message('[ERROR] Input must be type NIFTI file.')
        ogo.message('  {}'.format(input_image))
        os.sys.exit()
    
    # Define sheet_image if not set explicitly
    if sheet_image is None:
        basename = os.path.basename(input_image)
        name, ext = os.path.splitext(input_image)
        if 'gz' in ext:
            name = os.path.splitext(name)[0]  # Manages files with double extension
            ext = '.nii' + ext
        sheet_image = name + '_SHEET.nii.gz'
    
    # Define skin_image
    name, ext = os.path.splitext(sheet_image)
    if 'gz' in ext:
        name = os.path.splitext(name)[0]  # Manages files with double extension
        ext = '.nii' + ext
    skin_image = name.replace('_SHEET','') + '_SKIN.nii.gz'
    
    # Check if output sheet_image exists and should overwrite
    if os.path.isfile(sheet_image) and not overwrite:
        result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(sheet_image))
        if result.lower() not in ['y', 'yes']:
            ogo.message('Not overwriting. Exiting...')
            os.sys.exit()
    if not (sheet_image.lower().endswith('.nii') or sheet_image.lower().endswith('.nii.gz')):
        os.sys.exit('[ERROR] Output sheet image must be type NIFTI file: \"{}\"'.format(sheet_image))
    
    ogo.message('{}'.format('-------------- Files'))
    ogo.message('{:>7s}: {}'.format('input',input_image))
    ogo.message('{:>7s}: {}'.format('skin',skin_image))
    ogo.message('{:>7s}: {}'.format('sheet',sheet_image))
    ogo.message('{}'.format('--------- Parameters'))
    ogo.message('{:>15s}: {:>12s}'.format('enhance',str(enhance)))
    ogo.message('{:>15s}: {:12d}'.format('num_sigma',num_sigma))
    ogo.message('{:>15s}: {:12.2f}'.format('min_sigma',min_sigma))
    ogo.message('{:>15s}: {:12.2f}'.format('max_sigma',max_sigma))
    ogo.message('{:>15s}: {:12.2f}'.format('air_thres',air_thres))
    ogo.message('{:>15s}: {:12.2f}'.format('metal_thres',metal_thres))
    ogo.message('{:>15s}: {:12.2f}'.format('trace_weight',trace_weight))
    ogo.message('{}'.format('------------- Binary'))
    ogo.message('  {}'.format(sheetness_binary_path))
    ogo.message('Starting analysis...')

    if (enhance=='bright'):
        enhance_mode=True
    else:
        enhance_mode=False
            
    # Assemble the command for executing the Sheetness2 function
    cmd = [
      sheetness_binary_path, input_image, skin_image, sheet_image, enhance_mode,
      num_sigma, min_sigma, max_sigma, air_thres, metal_thres, trace_weight
    ]
    
    cmd = [str(x) for x in cmd]
    ogo.message('  CMD: {}'.format(cmd))
    exit()
    
    res = subprocess.check_output(cmd)
    
    #ogo.message('Results: {}'.format(res))
    ogo.message('{}'.format('------------ Outputs'))
    ogo.message('{:>15s}: {}'.format('sheet_image',sheet_image))
    ogo.message('{:>15s}: {}'.format('skin_image',skin_image))
    
    ogo.message('Done.')

# +------------------------------------------------------------------------------+
def periosteal(mark_image, sheet_image, peri_image, path_binaries, \
               gc_lambda, sigma, conn_filter, labels, cleanup, overwrite, func):

    # Check that the binaries are available to run the graphcuts algorithm
    periosteal_binary_path = ogo.find_executable('PeriostealSegmentation',path_binaries)
    if periosteal_binary_path is None:
        ogo.message('[ERROR] You must build the binaries for graph cuts for this to work.')
        ogo.message('        See instructions for installation here:')
        ogo.message('{}'.format(ogo.find_executable('OGO_GRAPHCUTS_INSTALL.txt',path_binaries)))
        os.sys.exit()
        
    # Check if input mark image exists
    if not os.path.isfile(mark_image):
        ogo.message('[ERROR] Cannot find input file.')
        ogo.message('  {}'.format(mark_image))
        os.sys.exit()
    if not (mark_image.lower().endswith('.nii') or mark_image.lower().endswith('.nii.gz')):
        ogo.message('[ERROR] Input must be type NIFTI file.')
        ogo.message('  {}'.format(mark_image))
        os.sys.exit()
    
    # Define sheet_image if not set explicitly
    if sheet_image is None:
        basename = os.path.basename(mark_image)
        name, ext = os.path.splitext(mark_image)
        if 'gz' in ext:
            name = os.path.splitext(name)[0]  # Manages files with double extension
            ext = '.nii' + ext
        sheet_image = name.replace('_MARK','') + '_SHEET.nii.gz'

    # Check if input sheet image exists
    if not os.path.isfile(sheet_image):
        ogo.message('[ERROR] Cannot find input sheet file.')
        ogo.message('  {}'.format(sheet_image))
        os.sys.exit()
    if not (sheet_image.lower().endswith('.nii') or sheet_image.lower().endswith('.nii.gz')):
        ogo.message('[ERROR] Input sheet image must be type NIFTI file.')
        ogo.message('  {}'.format(sheet_image))
        os.sys.exit()

    # Define peri_image if not set explicitly
    if peri_image is None:
        basename = os.path.basename(mark_image)
        name, ext = os.path.splitext(mark_image)
        if 'gz' in ext:
            name = os.path.splitext(name)[0]  # Manages files with double extension
            ext = '.nii' + ext
        peri_image = name.replace('_MARK','') + '_PERI.nii.gz'

    # Check if output peri exists and should overwrite
    if os.path.isfile(peri_image) and not overwrite:
        result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(peri_image))
        if result.lower() not in ['y', 'yes']:
            ogo.message('Not overwriting. Exiting...')
            os.sys.exit()
    
    # Read the input mark image so that we can check it is the correct type of image
    # (We expect either 8bit or 16bit unsigned images for the labels)
    ct_mark = sitk.ReadImage(mark_image, sitk.sitkUInt16)
    pixel_id = ct_mark.GetPixelID()
    if ((pixel_id is sitk.sitkUInt8) or \
        (pixel_id is sitk.sitkUInt16)):
        ogo.message('Input mark image type is {}'.format(ct_mark.GetPixelIDTypeAsString()))
    else:
        ogo.message('[ERROR] Input mark image with labels is usually unsigned 8 or 16 bit.')
        ogo.message('        Input mark image type is {}'.format(ct_mark.GetPixelIDTypeAsString()))
        ogo.message('  {}'.format(mark_image))
        os.sys.exit()

    # Establish the parameters
    ogo.message('{}'.format('-------------- Files'))
    ogo.message('{:>7s}: {}'.format('mark',mark_image))
    ogo.message('{:>7s}: {}'.format('sheet',sheet_image))
    ogo.message('{:>7s}: {}'.format('peri',peri_image))
    ogo.message('{}'.format('--------- Parameters'))
    ogo.message('{:>15s}: {:12.2f}'.format('gc_lambda',gc_lambda))
    ogo.message('{:>15s}: {:12.2f}'.format('sigma',sigma))
    ogo.message('{:>15s}: {:12d}'.format('conn_filter',conn_filter))
    ogo.message('{:>15s}: '.format('labels')+' '.join('{:2d}'.format(c_label) for c_label in labels))
    ogo.message('{}'.format('------------- Binary'))
    ogo.message('  {}'.format(periosteal_binary_path))
    
    ogo.message('Starting analysis...')
    ogo.message('  writing to working directory {}'.format(os.path.dirname(peri_image)))
    ogo.message('')

    # Gather all labels in image
    filt = sitk.LabelShapeStatisticsImageFilter()
    filt.Execute(ct_mark)
    n_labels = filt.GetNumberOfLabels()
    img_labels = filt.GetLabels()

    # Establish base name for temporary files
    temp_base_name, ext = os.path.splitext(peri_image)
    if 'gz' in ext:
        temp_base_name = os.path.splitext(temp_base_name)[0]  # Manages files with double extension
        ext = '.nii' + ext

    for label in labels:

        try:
            label_name = lb.labels_dict[label]['LABEL']
        except KeyError:
            label_name = 'no name'

        if label not in img_labels:
            ogo.message('[WARNING] Label {:d} ({}) not in marked image. No analysis possible.'.format(label,label_name))
            labels.remove(label)
    
        else:
            ogo.message('  processing label {:>2d} ({})'.format(label,label_name))

            temp_name = temp_base_name + '_TEMP_' + label_name.replace(' ','') + '.nii.gz'
            ogo.message('    {}'.format(os.path.basename(temp_name)))
            
            # Assemble the command for executing the Sheetness2 function
            cmd = [
              periosteal_binary_path, sheet_image, mark_image, temp_name,
              gc_lambda, sigma, label, conn_filter
            ]
            
            cmd = [str(x) for x in cmd]
            #ogo.message('  CMD: {}'.format(cmd))
            
            res = subprocess.check_output(cmd)
    
            #ogo.message('Results: {}'.format(res))
        ogo.message('')

    # Combine the images
    ogo.message('Combine the images...')
    ogo.message('  reading from working directory {}'.format(os.path.dirname(peri_image)))

    first_label = labels[0]
    label_name = lb.labels_dict[first_label]['LABEL']
    temp_name = temp_base_name + '_TEMP_' + label_name.replace(' ','') + '.nii.gz'
    ogo.message('    {}'.format(os.path.basename(temp_name)))
    
    seg = sitk.ReadImage(str(temp_name), sitk.sitkUInt8)
    seg = first_label*(seg>0)
    if cleanup:
        os.remove(temp_name)

    for label in labels:
        if label == first_label:
            continue
        
        label_name = lb.labels_dict[label]['LABEL']
        temp_name = temp_base_name + '_TEMP_' + label_name.replace(' ','') + '.nii.gz'
        ogo.message('    {}'.format(os.path.basename(temp_name)))
        
        this_label = sitk.ReadImage(str(temp_name), sitk.sitkUInt8)
        
        bin_seg = seg>0
        bin_this = this_label>0
        overlap = (bin_seg + bin_this)==2
        mask = 1 - overlap
        this_label = sitk.Mask(this_label, mask)
        
        seg = seg + label*(this_label>0)
        
        if cleanup:
            os.remove(temp_name)
    
    ogo.message('')
    ogo.message('Writing output file.')
    ogo.message('    {}'.format(peri_image))
    sitk.WriteImage(seg, peri_image)
    
    ogo.message('Done.')
    
    
def main():
    # Setup description
    description = '''
The graph cuts algorithm is used to segment objects in CT scans and it
involves two steps. First, a \'sheetness\' algorithm pre-processes the CT 
image. Second, a labelled image is fed into the algorithm to define the
periosteal surfaces for each labelled bone. You must create the sheetness
image before you can generate the periosteal surfaces.

This python program is only a convenience tool for accessing the graph cuts
software available here:
https://gridcut.com/downloads.php

It is necessary to install and compile the binaries for graph cuts before 
this script will work. To learn more, read OGO_GRAPHCUTS_INSTALL.txt located
in Ogo/ogo/util/graphcuts . If you try and execute ogoGraphCuts without 
installing the binaries first you will be redirected OGO_GRAPHCUTS_INSTALL.txt

sheetness – Computes an image enhancement using eigenvalues of the local 
            Hessian matrix over many scales.
'''

    epilog = '''
Cite:

Boykov Y, Funka-Lea G, 2006. Graph cuts and efficient N-D image segmentation. 
Int J Comput Vision 70, 109-131. doi = 10.1007/s11263-006-7934-5 

Example calls: 
ogoGraphCuts sheetness image.nii.gz
ogoGraphCuts periosteal image_MARK.nii.gz
'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoGraphCuts",
        description=description,
        epilog=epilog
    )
    subparsers = parser.add_subparsers()

    # Sheetness
    parser_sheetness = subparsers.add_parser('sheetness')
    parser_sheetness.add_argument('input_image', metavar='FILE IN', help='Input raw CT image file (*.nii, *.nii.gz)')
    parser_sheetness.add_argument('--sheet_image', metavar='FILE', help='Output sheetness image (*.nii, *.nii.gz) (default: input_image_SHEET.nii.gz)')
    parser_sheetness.add_argument('--path_binaries', metavar='PATH', default='/Users/', help='Start search path for graphcut binary (default: %(default)s)')
    parser_sheetness.add_argument('--enhance', default='bright', choices=['bright', 'dark'], help='Select enhancement mode (default: %(default)s)')
    parser_sheetness.add_argument('--num_sigma', metavar='INT', type=int, default=2, help='Number of sigma steps to use for enhancing in Hessian filter (default: %(default)s)')
    parser_sheetness.add_argument('--min_sigma', metavar='FLOAT', type=float, default=0.5, help='Minimum sigma (default: %(default)s)')
    parser_sheetness.add_argument('--max_sigma', metavar='FLOAT', type=float, default=1.0, help='Maximum sigma (default: %(default)s)')
    parser_sheetness.add_argument('--air_thres', metavar='FLOAT', type=float, default=-400.0, help='Threshold for determining air (default: %(default)s)')
    parser_sheetness.add_argument('--metal_thres', metavar='FLOAT', type=float, default=1200.0, help='Threshold for determining metal (default: %(default)s)')
    parser_sheetness.add_argument('--trace_weight', metavar='FLOAT', type=float, default=0.05, help='Frobenius norm weight for reducing noise (default: %(default)s)')
    parser_sheetness.add_argument('--overwrite', action='store_true', help='Overwrite output image without asking (default: %(default)s)')
    parser_sheetness.set_defaults(func=sheetness)

    # Periosteal
    parser_periosteal = subparsers.add_parser('periosteal')
    parser_periosteal.add_argument('mark_image', metavar='FILE IN', help='Input marked image mask file (*.nii, *.nii.gz) (typically input_image_MARK.nii.gz)')
    parser_periosteal.add_argument('--sheet_image', metavar='FILE', help='Input sheetness CT image file (*.nii, *.nii.gz)  (default: input_image_SHEET.nii.gz)')
    parser_periosteal.add_argument('--peri_image', metavar='FILE', help='Output segmented image file (*.nii, *.nii.gz)  (default: input_image_PERI.nii.gz)')
    parser_periosteal.add_argument('--path_binaries', metavar='PATH', default='/Users/', help='Start search path for graphcut binary (default: %(default)s)')
    parser_periosteal.add_argument('--gc_lambda', metavar='FLOAT', type=float, default=50.0, help='Smoothness term (default: %(default)s)')
    parser_periosteal.add_argument('--sigma', metavar='FLOAT', type=float, default=0.25, help='Boundry term noise (default: %(default)s)')
    parser_periosteal.add_argument('--conn_filter', metavar='INT', type=int, default=1, help='Number of connected bones to detect (default: %(default)s)')
    parser_periosteal.add_argument('--labels', type=int, nargs='*', default=[1,2,3,4,5,6,7,8,9,10], metavar='LABEL', help='List of labels for analysis (default: %(default)s)')
    parser_periosteal.add_argument('--cleanup', action='store_true', help='Deletes individual segmentations (default: %(default)s)')
    parser_periosteal.add_argument('--overwrite', action='store_true', help='Overwrite output image without asking (default: %(default)s)')
    parser_periosteal.set_defaults(func=periosteal)

    # Parse and display
    args = parser.parse_args()
    #print(echo_arguments('GraphCuts', vars(args)))

    # Run program
    args.func(**vars(args))


if __name__ == '__main__':
    main()
