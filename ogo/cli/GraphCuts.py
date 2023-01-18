# /------------------------------------------------------------------------------+
# | 16-JAN-2023                                                                  |
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
def sheetness(input_image, sheet_image, path_binaries, enhance_bright, \
              num_sigma, min_sigma, max_sigma, air_thres, metal_thres, trace_weight, \
              overwrite, func):

    # Check that the binaries are available to run the graphcuts algorithm
    sheetness_binary_path = ogo.find_executable('Sheetness2',path_binaries)
    if sheetness_binary_path is None:
        ogo.message('[ERROR] You must have built the binaries for graph cuts for this to work.')
        ogo.message('        See instructions for installation in Ogo/ogo/util/graphcuts/.')
    
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
    ogo.message('{:>15s}: {:>12s}'.format('enhance_bright',str(enhance_bright)))
    ogo.message('{:>15s}: {:12d}'.format('num_sigma',num_sigma))
    ogo.message('{:>15s}: {:12.2f}'.format('min_sigma',min_sigma))
    ogo.message('{:>15s}: {:12.2f}'.format('max_sigma',max_sigma))
    ogo.message('{:>15s}: {:12.2f}'.format('air_thres',air_thres))
    ogo.message('{:>15s}: {:12.2f}'.format('metal_thres',metal_thres))
    ogo.message('{:>15s}: {:12.2f}'.format('trace_weight',trace_weight))
    ogo.message('{}'.format('------------- Binary'))
    ogo.message('  {}'.format(sheetness_binary_path))
    ogo.message('Starting analysis...')
    
    # Assemble the command for executing the Sheetness2 function
    cmd = [
      sheetness_binary_path, input_image, skin_image, sheet_image, enhance_bright,
      num_sigma, min_sigma, max_sigma, air_thres, metal_thres, trace_weight
    ]
    
    cmd = [str(x) for x in cmd]
    #ogo.message('  CMD: {}'.format(cmd))

    res = subprocess.check_output(cmd)
    
    #ogo.message('Results: {}'.format(res))
    ogo.message('{}'.format('------------ Outputs'))
    ogo.message('{:>15s}: {}'.format('sheet_image',sheet_image))
    ogo.message('{:>15s}: {}'.format('skin_image',skin_image))
    
    ogo.message('Done.')

# +------------------------------------------------------------------------------+
def periosteal(mark_image, sheet_image, output_image, overwrite, func):

    # Check if output exists and should overwrite
    if os.path.isfile(output_image) and not overwrite:
        result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(output_image))
        if result.lower() not in ['y', 'yes']:
            ogo.message('Not overwriting. Exiting...')
            os.sys.exit()
    
    # Read input mark image
    if not os.path.isfile(mark_image):
        os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(mark_image))

    if not (mark_image.lower().endswith('.nii') or mark_image.lower().endswith('.nii.gz')):
        os.sys.exit('[ERROR] Input mark image must be type NIFTI file: \"{}\"'.format(mark_image))
    
    ogo.message('Reading image: ')
    ogo.message('\"{}\"'.format(mark_image))
    ct_mark = sitk.ReadImage(mark_image, sitk.sitkUInt8)
    
    # Check if sheetness image exists
    if os.path.isfile(sheet_image):
        ogo.message('Found image: ')
        ogo.message('\"{}\"'.format(sheet_image))
        
    exit()
    
    ogo.message('Done.')
    
    
def main():
    # Setup description
    description = '''
The graph cuts algorithm is used to segment objects in CT scans and it
involves two steps. First, a \'sheetness\' algorithm pre-processes the CT 
image. Second, a roughly labelled image is fed into the algorithm and is the 
basis to complete the segmentation.

This python program is only a convenience tool for accessing the graph cuts
software available here:
https://gridcut.com/downloads.php

It is necessary to install this software before executing graph cuts. To
learn more run with --installation_notes.
'''

    epilog = '''
Cite:

Boykov Y, Funka-Lea G, 2006. Graph cuts and efficient N-D image segmentation. 
Int J Comput Vision 70, 109-131. doi = 10.1007/s11263-006-7934-5 

Example calls: 
ogoGraphCuts sheetness image.nii.gz
ogoGraphCuts periosteal image.nii.gz image_sheet.nii.gz

ogoGraphCuts sheetness /Users/skboyd/Desktop/ML/test/kub.nii.gz
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
    parser_sheetness.add_argument('input_image', help='Input raw CT image file (*.nii, *.nii.gz)')
    parser_sheetness.add_argument('--sheet_image', metavar='FILE', help='Output sheetness image (*.nii, *.nii.gz) (default: input_image_SHEET.nii.gz)')
    parser_sheetness.add_argument('--path_binaries', metavar='PATH', default='/Users/', help='Start search path for graphcut binary (default: %(default)s)')
    parser_sheetness.add_argument('--enhance_bright', action='store_false', help='Enhance bright objects (default: %(default)s)')
    parser_sheetness.add_argument('--num_sigma', metavar='INT', type=int, default=2, help='How many sigmas to use for enhancing (default: %(default)s)')
    parser_sheetness.add_argument('--min_sigma', metavar='FLOAT', type=float, default=0.5, help='How many sigmas to use for enhancing (default: %(default)s)')
    parser_sheetness.add_argument('--max_sigma', metavar='FLOAT', type=float, default=1.0, help='How many sigmas to use for enhancing (default: %(default)s)')
    parser_sheetness.add_argument('--air_thres', metavar='FLOAT', type=float, default=-400.0, help='Threshold for determining air (default: %(default)s)')
    parser_sheetness.add_argument('--metal_thres', metavar='FLOAT', type=float, default=1200.0, help='Threshold for determining metal (default: %(default)s)')
    parser_sheetness.add_argument('--trace_weight', metavar='FLOAT', type=float, default=0.05, help='Weight for reducing noise (default: %(default)s)')
    parser_sheetness.add_argument('--overwrite', action='store_true', help='Overwrite output image without asking')
    parser_sheetness.set_defaults(func=sheetness)

    # Periosteal
    parser_periosteal = subparsers.add_parser('periosteal')
    parser_periosteal.add_argument('mark_image', help='Input marked image mask file (*.nii, *.nii.gz) (typically _MARK.nii.gz)')
    parser_periosteal.add_argument('sheet_image', help='Input sheetness CT image file (*.nii, *.nii.gz)')
    parser_periosteal.add_argument('output_image', help='Output segmented image file (*.nii, *.nii.gz)')
    parser_periosteal.add_argument('--overwrite', action='store_true', help='Overwrite output image without asking')
    parser_periosteal.set_defaults(func=periosteal)

    # Parse and display
    args = parser.parse_args()
    #print(echo_arguments('GraphCuts', vars(args)))

    # Run program
    args.func(**vars(args))


if __name__ == '__main__':
    main()
