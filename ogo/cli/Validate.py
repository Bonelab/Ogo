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

def get_labels(ct):
    filt = sitk.LabelShapeStatisticsImageFilter()
    filt.Execute(ct)
    labels = filt.GetLabels()
    return labels

# +------------------------------------------------------------------------------+
def validate(input_image, report_file, print_parts, expected_labels, overwrite, func):
    
    passed_all_tests = True
    validate_expected_labels_in_image = True
    validate_singulary_connected_bone = True
    #validate_positions_of_bones = True   # implement later?
    #validate_bone_size_within_range = True             # implement later?
    
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
            
    # Gather all labels in image and determine number of parts each label is broken into
    filt = sitk.LabelShapeStatisticsImageFilter() # Used to get labels in image
    filt.Execute(ct)
    n_labels = filt.GetNumberOfLabels()
    labels = filt.GetLabels()
    label_repair_list = [] # list of (label, part, new_label)
    
    thres = sitk.BinaryThresholdImageFilter() # Used to isolate labels and find # of parts
    thres.SetInsideValue(1)
    thres.SetOutsideValue(0)
    conn = sitk.ConnectedComponentImageFilter()
    conn.SetFullyConnected(True)
    stats = sitk.LabelIntensityStatisticsImageFilter()
    
    # Check if expected labels are in image
    test_labels = []
    for label in expected_labels:
        if label in labels:
            test_labels.append(label)
        else:
            try:
                label_name = lb.labels_dict[label]['LABEL']
            except KeyError:
                label_name = 'unknown label'
            ogo.message('[WARNING] Expected label {} ({}) was not found in image.'.format(label,label_name))
            validate_expected_labels_in_image = False
    if test_labels == []:
        os.sys.exit('[ERROR] None of the expected labels were found in image.')
    labels = test_labels
    
    ogo.message('Gather each label in the image and determine number of parts.')
    
    report += '  {:>27s}\n'.format('_____________________Report')
    report += '  {:>27s} {:s}\n'.format('number of labels:',str(n_labels))
    report += '  {:>27s} {:>8s} {:>6s} {:>10s} {:>22s} {:>6s}\n'.format('LABEL','SIZE','PART','VOL','CENTROID','OFFSET')
    report += '  {:>27s} {:>8s} {:>6s} {:>10s} {:>22s} {:>6s}\n'.format('#','voxels','#','mm3','(mm,mm,mm)','mm')
        
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
            sitk.WriteImage(ct_conn_sorted, label_fname) # ERROR: Should write image with original label (replacelabel)
        
        report += '  {:>22s}{:>5s} {:8.0f} '.format('('+desc+')',str(label),size)        
        # Generate report data
        for part in stats.GetLabels(): # for each part of a given label
            c = stats.GetCentroid(part)
            if part==1:
                ref_c = c
                report += '{:6d} {:10.1f} ({:6.1f},{:6.1f},{:6.1f})\n'.format(part,stats.GetPhysicalSize(part),c[0],c[1],c[2])
            else:
                validate_singulary_connected_bone = False
                distance_between_centroids = np.sqrt((ref_c[0] - c[0])**2 + (ref_c[1] - c[1])**2 + (ref_c[2] - c[2])**2)
                report += '  {:>22s}{:>5s} {:8s} {:6d} {:10.1f} ({:6.1f},{:6.1f},{:6.1f}) {:6.1f}\n'\
                          .format('','','',part,stats.GetPhysicalSize(part),c[0],c[1],c[2],distance_between_centroids)
        
        # Generate suggestion of new label for command line suggestion
        if n_parts>1:
            for part in stats.GetLabels():
                if part>1:
                    swap = [label,part,label]
                    label_repair_list.append(swap)
    
    # Final report
    report += '\n'
    report += '  {:>27s}\n'.format('_______________Final report')
    if validate_expected_labels_in_image:
        report += '  {:>27s} {:s}\n'.format('Each expected label found:','PASSED')
    else:
        report += '  {:>27s} {:s}\n'.format('Each expected label found:','FAILED')
    if validate_singulary_connected_bone:
        report += '  {:>27s} {:s}\n'.format('Each label a single part:','PASSED')
    else:
        report += '  {:>27s} {:s}\n'.format('Each label a single part:','FAILED')
    report += '\n'
    if validate_singulary_connected_bone and validate_expected_labels_in_image:
        passed_all_tests = True
    else:
        passed_all_tests = False
    if passed_all_tests:
        report += '  {:>27s} {:s}\n'.format('ALL:','PASSED')
    else:
        report += '  {:>27s} {:s}\n'.format('ALL:','FAILED')
    
    # Command line
    cmd_line = ''
    cmd_line += '  {:>27s}\n'.format('_______________Command line')
    cmd_line += '\n'
    cmd_line += '  {}\n'.format('Edit the command below. Should replaced labels be zero?')
    cmd_line += '\n'
    cmd_line += '  {}\n'.format('ogoValidate repair \\')
    cmd_line += '    {} \\\n'.format(input_image)    
    cmd_line += '    {}_REPAIR{} \\\n'.format(name,ext)
    cmd_line += '    {}\n'.format('--relabel_parts \\')
    for idx,swap in enumerate(label_repair_list):
        cmd_line += '    '+' '.join('{:d}'.format(i) for i in swap)
        if idx<len(label_repair_list)-1:
            cmd_line += ' \\\n'
        else:
            cmd_line += '\n'
             
    # # Erode
    # erode = sitk.BinaryErodeImageFilter()
    # erode.SetKernelRadius(4)
    # erode.SetForegroundValue(1)
    # ct_erode = erode.Execute(ct)

    ogo.message('Printing report')
    print('\n')
    print(report)
    if report_file:
        ogo.message('Saving report to file:')
        ogo.message('      \"{}\"'.format(report_file))
        txt_file = open(report_file, "w")
        txt_file.write(report)
        txt_file.close()

    if label_repair_list != []:
        print(cmd_line)
    
    ogo.message('Done.')

# +------------------------------------------------------------------------------+
def repair(input_image, output_image, relabel_parts, remove_by_volume, overwrite, func):

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
    
    # Get all the available labels in the input image
    labels = get_labels(ct)
    n_labels = len(labels)
    ogo.message('Input image contains the following labels:')
    ogo.message('  [' + ', '.join('{:d}'.format(i) for i in labels) + ']')

    # Prepare to go through labels and parts
    filt = sitk.LabelShapeStatisticsImageFilter()
    filt.Execute(ct)
    n_labels = filt.GetNumberOfLabels()
    labels = filt.GetLabels()
    
    labelmapfilt = sitk.LabelImageToLabelMapFilter()
    labelimagefilt = sitk.LabelMapToLabelImageFilter()
    mergefilt = sitk.MergeLabelMapFilter()
    mergefilt.Aggregate
    #mergefilt.Strict
    #mergefilt.Keep
    #mergefilt.Pack
    
    #sitk.ChangeLabelImageFilter() - used in ogoMerge
    #sitk.ChangeLabelLabelMapFilter() – not used. why?
    
    thres_label = sitk.BinaryThresholdImageFilter() # Used to isolate labels and find # of parts
    thres_label.SetInsideValue(1)
    thres_label.SetOutsideValue(0)
    thres_part = sitk.BinaryThresholdImageFilter() # Used to isolate labels and find # of parts
    conn = sitk.ConnectedComponentImageFilter()
    conn.SetFullyConnected(True)
    stats = sitk.LabelIntensityStatisticsImageFilter()    
    
    if relabel_parts == [] and remove_by_volume == 0:
        ogo.message('No operations requested. Suggest setting one of the follwoing:')
        ogo.message('  --relabel_parts')
        ogo.message('  --remove_by_volume')
        
    if (relabel_parts):
        ct_base = sitk.Image(ct) # Do we need to make a deep copy
        ct_base = labelmapfilt.Execute(ct_base)
        
        if np.remainder(len(relabel_parts),3) != 0:
            os.sys.exit('[ERROR] Incorrect argument for --relabel_parts. Tuples of 3 expected.')
        n_relabel_parts = int(len(relabel_parts) / 3)
        relabel_parts = np.reshape(relabel_parts, (n_relabel_parts,3))
        ogo.message('Relabelling {:d} parts in image:'.format(n_relabel_parts))
        for swap in relabel_parts:
            label = int(swap[0])
            part = int(swap[1])
            new_label = int(swap[2])
            
            thres_label.SetUpperThreshold(label)
            thres_label.SetLowerThreshold(label)
            ct_thres = thres_label.Execute(ct)
            ct_conn = conn.Execute(ct,ct_thres)
            ct_conn_sorted = sitk.RelabelComponent(ct_conn, sortByObjectSize=True) # could use minimumObjectSize
            
            filt.Execute(ct_conn_sorted)
            size = filt.GetPhysicalSize(part)
            
            thres_part.SetUpperThreshold(part)
            thres_part.SetLowerThreshold(part)
            thres_part.SetInsideValue(new_label)
            thres_part.SetOutsideValue(0)
            ct_part = thres_part.Execute(ct_conn_sorted)

            ogo.message('  label {:d} ({:s}), part {} ({:.1f} mm3) --> {} ({})'\
                        .format(label,lb.labels_dict[label]['LABEL'],part,size,new_label,lb.labels_dict[new_label]['LABEL']))
            
            ct_part = labelmapfilt.Execute(ct_part)
            ct_base = mergefilt.Execute(ct_part,ct_base)

        ct_final = labelimagefilt.Execute(ct_base)
        
    final_labels = get_labels(ct_final)
    ogo.message('Final image contains the following labels:')
    ogo.message('  [' + ', '.join('{:d}'.format(i) for i in final_labels) + ']')
    
    ogo.message('')
    ogo.message('Writing merged output image to file:')
    ogo.message('  {}'.format(output_image))
    sitk.WriteImage(ct_final, output_image)            

    ogo.message('Done.')
    
    
def main():
    # Setup description
    description = '''
Validates segmented output of machine learning tools of the skeleton
and provides some basic tools for repairing of errors. The 'validate'
and 'repair' options are separated.

The validation performs the following checks:
  
  Each label expected in the image is found
  Each label is a single part
  Each bone size is within expected range (not yet implemented)
  Each bone position is as expected (not yet implemented)

The repair options include:
  Relabel all voxels of a specified bone and part
  Remove fragments by size

'''

    epilog = '''
Example calls: 
ogoValidate check image.nii.gz

ogoValidate repair image.nii.gz --component 13 --label 34
ogoValidate repair image.nii.gz --component 13 --label 0

ogoValidate validate /Users/skboyd/Desktop/ML/test/robust/frag/RETRO_01638_NNUNET.nii.gz
ogoValidate validate /Users/skboyd/Desktop/ML/test/robust/frag/RETRO_02317_NNUNET.nii.gz
ogoValidate validate /Users/skboyd/Desktop/ML/test/robust/mixed/RETRO_01964_NNUNET.nii.gz

ogoValidate validate /Users/skboyd/Desktop/ML/test/robust/perfect/RETRO_TRAIN_CORR.nii.gz
ogoValidate repair \
    /Users/skboyd/Desktop/ML/test/robust/perfect/RETRO_TRAIN_CORR.nii.gz \
    /Users/skboyd/Desktop/ML/test/robust/perfect/RETRO_TRAIN_CORR_REPAIR.nii.gz \
    --relabel_parts \
    3 2 3 \
    5 2 5 \
    6 2 6 \
    9 2 9

ogoValidate validate /Users/skboyd/Desktop/ML/test/robust/frag/RETRO_02317_NNUNET.nii.gz
ogoValidate repair \
    /Users/skboyd/Desktop/ML/test/robust/frag/RETRO_02317_NNUNET.nii.gz \
    /Users/skboyd/Desktop/ML/test/robust/frag/RETRO_02317_NNUNET_REPAIR.nii.gz \
    --relabel_parts \
    3 2 3 \
    3 3 3 \
    3 4 3 \
    4 2 4 \
    5 2 5 \
    5 3 5 \
    5 4 5 \
    5 5 5 \
    5 6 5 \
    5 7 5 \
    5 8 5 \
    5 9 5 \
    5 10 5 \
    6 2 6 \
    6 3 6
ogoValidate repair \
    /Users/skboyd/Desktop/ML/test/robust/frag/RETRO_02317_NNUNET.nii.gz \
    /Users/skboyd/Desktop/ML/test/robust/frag/RETRO_02317_NNUNET_REPAIR.nii.gz \
    --relabel_parts \
    5 2 5 \
    5 3 5 

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
    parser_validate.add_argument('--expected_labels', type=int, nargs='*', default=[1,2,3,4,5,6,7,8,9,10], metavar='LABEL', help='List of labels expected in image (default: %(default)s)')
    parser_validate.add_argument('--overwrite', action='store_true', help='Overwrite validation report without asking')
    parser_validate.set_defaults(func=validate)

    # Repair
    parser_repair = subparsers.add_parser('repair')
    parser_repair.add_argument('input_image', help='Input image file (*.nii, *.nii.gz)')
    parser_repair.add_argument('output_image', help='Output image file (*.nii, *.nii.gz)')
    parser_repair.add_argument('--relabel_parts', type=int, nargs='*', default=[], metavar='LABEL PART NEWLABEL', help='List of labels, parts, and new labels; space separated. Example shows relabelling L4 part 2 to  L3: (e.g. 7 2 8)')
    parser_repair.add_argument('--remove_by_volume', type=float, nargs=1, default=0, metavar='VOL', help='Removes all parts by volumens threshold (default: %(default)s)')
    parser_repair.add_argument('--overwrite', action='store_true', help='Overwrite output image without asking')
    parser_repair.set_defaults(func=repair)

    # Parse and display
    args = parser.parse_args()
    print(echo_arguments('Validate', vars(args)))

    # Run program
    args.func(**vars(args))


if __name__ == '__main__':
    main()
