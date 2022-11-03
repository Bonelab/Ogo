# /------------------------------------------------------------------------------+
# | 24-OCT-2022                                                                  |
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
import numpy as np
import SimpleITK as sitk
from datetime import date
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from ogo.util.echo_arguments import echo_arguments
import ogo.util.Helper as ogo
import ogo.dat.OgoMasterLabels as lb

# +------------------------------------------------------------------------------+
def validate(input_image, report_file, print_parts, overwrite, func):
    
    # Check if output exists and should overwrite
    if not report_file is None:
        if os.path.isfile(report_file) and not overwrite:
            result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(outfile))
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
    ct = sitk.ReadImage(input_image)
    
    # Create base filename (and some extra information)
    basename = os.path.basename(input_image)
    name, ext = os.path.splitext(input_image)
    if 'gz' in ext:
        name = os.path.splitext(name)[0]  # Manages files with double extension
        ext = '.nii' + ext
    dirname = os.path.dirname(input_image)
    
    # Start the report
    report = ''
    report += '  {:>27s}\n'.format('_________________Validation')
    report += '  {:>27s} {:s}\n'.format('ID:', os.path.basename(input_image))
    report += '  {:>27s} {:s}\n'.format('python script:', os.path.splitext(os.path.basename(sys.argv[0]))[0])
    report += '  {:>27s} {:.2f}\n'.format('version:', script_version)
    report += '  {:>27s} {:s}\n'.format('creation date:', str(date.today()))
    report += '\n'
    
    # Command line
    cmd_line = ''
    cmd_line += '  {:>27s}\n'.format('_________________Command line')
    cmd_line += '\n'
    cmd_line += '  {}\n'.format('Edit the command below. Replace NL with the new labels.')
    cmd_line += '\n'
    #cmd_line += '  {}\n'.format('ogoValidate repair \\')
    cmd_line += '  {}\n'.format('python ../../code/Ogo/ogo/cli/Validate.py repair \\')
    cmd_line += '    {} \\\n'.format(input_image)    
    cmd_line += '    {}_REPAIR{} \\\n'.format(name,ext)
    cmd_line += '    {}\n'.format('--relabel_parts \\')
        
    # Gather all labels in image and determine number of parts each label is broken into
    ogo.message('Gather list labels in image and determine number of parts.')
    filt = sitk.LabelShapeStatisticsImageFilter() # Used to get labels in image
    filt.Execute(ct)
    n_labels = filt.GetNumberOfLabels()
    labels = filt.GetLabels()
    
    thres = sitk.BinaryThresholdImageFilter() # Used to isolate labels and find # of parts
    thres.SetInsideValue(1)
    thres.SetOutsideValue(0)
    conn = sitk.ConnectedComponentImageFilter()
    conn.SetFullyConnected(True)
    stats = sitk.LabelIntensityStatisticsImageFilter()
    
    report += '  {:>27s}\n'.format('___________________Report')
    report += '  {:>27s} {:s}\n'.format('number of labels:',str(n_labels))
    report += '  {:>27s} {:>8s} {:>6s} {:>10s} {:>22s}\n'.format('LABEL','SIZE','PART','PART_Sz','CENTROID')
        
    for idx,label in enumerate(labels):
        ogo.message('  processing label {} ({})'.format(label,lb.labels_dict[label]['LABEL']))
        desc = lb.labels_dict[label]['LABEL']
        #centroid = filt.GetCentroid(label)
        size = filt.GetPhysicalSize(label)
        
        thres.SetUpperThreshold(label)
        thres.SetLowerThreshold(label)
        ct_thres = thres.Execute(ct)
        
        ct_conn = conn.Execute(ct,ct_thres)
        
        ct_conn_sorted = sitk.RelabelComponent(ct_conn, sortByObjectSize=True) # could use minimumObjectSize
        
        stats.Execute(ct_conn_sorted,ct)
        n_parts = stats.GetNumberOfLabels()
        
        if (print_parts):
            label_fname = '{}_lab{:02}{}'.format(name,label,ext)
            sitk.WriteImage(ct_conn_sorted, label_fname)
        
        report += '  {:>22s}{:>5s} {:8.0f} '.format('('+desc+')',str(label),size)        
        # Generate report data
        for part in stats.GetLabels(): # for each part of a given label
            c = stats.GetCentroid(part)
            if part>1:
                report += '  {:>22s}{:>5s} {:8s} {:6d} {:10.1f} ({:6.1f},{:6.1f},{:6.1f})\n'.format('','','',part,stats.GetPhysicalSize(part),c[0],c[1],c[2])
            else:
                report += '{:6d} {:10.1f} ({:6.1f},{:6.1f},{:6.1f})\n'.format(part,stats.GetPhysicalSize(part),c[0],c[1],c[2])
        
        # Generate command line
        if n_parts>1 and label <80:
            cmd_line += '    {:d} \\\n'.format(label)                  # current label
            for part in stats.GetLabels(): 
                cmd_line += '      {:d} {:s} \\\n'.format(part,'NL')      # part number and new label

    # largest_component_binary_image = sorted_component_image == 1   # gets the component 1 from the image

    
    # # Erode
    # erode = sitk.BinaryErodeImageFilter()
    # erode.SetKernelRadius(4)
    # erode.SetForegroundValue(1)
    # ct_erode = erode.Execute(ct)

    ogo.message('Printing report')
    print('\n')
    print(report)
    print(cmd_line)
    exit()


    # EXAMPLE 1
    print('EXAMPLE 1')
    binary_image = sitk.Image([256,256], sitk.sitkUInt8)
    binary_image[20:50, 25:125] = 1
    binary_image[70:95, 40:60] = 1
    binary_image[70:95, 80:140] = 1

    # 1. Convert binary image into a connected component image, each component has an integer label.
    # 2. Relabel components so that they are sorted according to size (there is an
    #    optional minimumObjectSize parameter to get rid of small components).
    # 3. Get largest connected componet, label==1 in sorted component image.
    component_image = sitk.ConnectedComponent(binary_image)
    sorted_component_image = sitk.RelabelComponent(component_image, sortByObjectSize=True)
    largest_component_binary_image = sorted_component_image == 1

    # Report results
    print(report)                                                                
    if report_file:
        ogo.message('Saving validate report to file:')
        ogo.message('      \"{}\"'.format(report_file))
        txt_file = open(report_file, "w")
        txt_file.write(report)
        txt_file.close()
    
    ogo.message('Done.')

# +------------------------------------------------------------------------------+
def repair(input_image, output_image, relabel_parts, overwrite, func):

    # Check if output exists and should overwrite
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
    ct = sitk.ReadImage(input_image)
    
    #ogo.message('Done.')
    
    
def main():
    # Setup description
    description = '''
Validates segmented output of machine learning tools of the skeleton
and provides some basic tools for repairing of errors. The 'validate'
and 'repair' options are separated.

The validation performs the following checks:
  Number of connected components per bone label (should be 1)
  
  Each bone labelled with just one label
  Each bone labelled the correct label
  Number of bones matches expectation
  Bones are abnormally small or large

The repair options include:
  Relabel all voxels of a specified bone (connected component)
  Remove fragments (by size or by label)

'''

    epilog = '''
Example calls: 
ogoValidate check image.nii.gz

ogoValidate repair image.nii.gz --component 13 --label 34
ogoValidate repair image.nii.gz --component 13 --label 0

ogoValidate validate /Users/skboyd/Desktop/ML/test/robust/perfect/RETRO_TRAIN_CORR.nii.gz
ogoValidate validate /Users/skboyd/Desktop/ML/test/robust/frag/RETRO_01638_NNUNET.nii.gz
ogoValidate validate /Users/skboyd/Desktop/ML/test/robust/frag/RETRO_02317_NNUNET.nii.gz
ogoValidate validate /Users/skboyd/Desktop/ML/test/robust/mixed/RETRO_01964_NNUNET.nii.gz

'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoValidate",
        description=description,
        epilog=epilog
    )
    subparsers = parser.add_subparsers()

    # Validate
    parser_validate = subparsers.add_parser('validate')
    parser_validate.add_argument('input_image', help='Input image file (*.nii, *.nii.gz)')
    parser_validate.add_argument('--report_file', metavar='FILE', help='Write validation report to file (*.txt)')
    parser_validate.add_argument('--print_parts', action='store_true', help='Writes N label output images showing component parts')
    parser_validate.add_argument('--overwrite', action='store_true', help='Overwrite validation report without asking')
    parser_validate.set_defaults(func=validate)

    # Repair
    parser_repair = subparsers.add_parser('repair')
    parser_repair.add_argument('input_image', help='Input image file (*.nii, *.nii.gz)')
    parser_repair.add_argument('output_image', help='Output image file (*.nii, *.nii.gz)')
    parser_repair.add_argument('--relabel_parts', type=int, nargs='*', default=[0], metavar='LABEL PART NEWLABEL', help='List of labels, parts, and new labels; space separated. Example shows relabelling L4 part 2 to  L3: (e.g. 7 2 8)')
    #parser_repair.add_argument('--component', type=int, default=0, metavar='#',help='Connected component number (default: %(default)s)')
    #parser_repair.add_argument('--label', type=int, default=0, metavar='#',help='Label of bone (default: %(default)s)')
    parser_repair.add_argument('--overwrite', action='store_true', help='Overwrite output image without asking')
    parser_repair.set_defaults(func=repair)

    # Parse and display
    args = parser.parse_args()
    print(echo_arguments('Validate', vars(args)))

    # Run program
    args.func(**vars(args))


if __name__ == '__main__':
    main()
