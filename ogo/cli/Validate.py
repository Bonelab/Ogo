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
import math
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

def is_smaller(filt,percent_threshold,label1,label2):
    
    if label1 not in filt.GetLabels():
        ogo.message('[WARNING] Cannot check if smaller because label {} is not in image.'.format(label1))
        return False
    if label2 not in filt.GetLabels():
        ogo.message('[WARNING] Cannot check if smaller because label {} is not in image.'.format(label2))
        return False
    
    size1 = filt.GetPhysicalSize(label1)
    size2 = filt.GetPhysicalSize(label2)
    average_size = (size1+size2)/2
    thres = average_size * percent_threshold / 100.0
    
    #print('size1 = {:10.3f}'.format(size1))
    #print('size2 = {:10.3f}'.format(size2))
    #print('thres = {:10.3f}'.format(thres))
    #print('diff = {:10.3f}'.format((size1 - size2)))
    
    if (size2 - size1) > thres:
        return True
    else:
        return False
        
# +------------------------------------------------------------------------------+
def validate(input_image, report_file, print_parts, expected_labels, overwrite, func):
    
    # Validation is based on 'innocent until proven guilty' (i.e. default is True)
    QA_labels = {}
    QualityAssurance = {
      "All Pass":  True,
      "Expected Labels": True,
      "Labels":   QA_labels
    }

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
    ct = sitk.ReadImage(input_image, sitk.sitkUInt8)
    
    # Create base filename (and some extra information)
    basename = os.path.basename(input_image)
    name, ext = os.path.splitext(input_image)
    if 'gz' in ext:
        name = os.path.splitext(name)[0]  # Manages files with double extension
        ext = '.nii' + ext
    dirname = os.path.dirname(input_image)
    
    # Start the report
    report = ''
    report += '  {:>27s}\n'.format('_____________________________________________________________________________Validation')
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
    
    conn = sitk.ConnectedComponentImageFilter()
    conn.SetFullyConnected(True)
    stats = sitk.LabelIntensityStatisticsImageFilter()
    
    # Check if expected labels are in image
    test_labels = []
    QualityAssurance["Expected Labels"] = True
    for label in expected_labels:
        if label in labels:
            test_labels.append(label)
            QA_labels[label] = True
        else:
            try:
                label_name = lb.labels_dict[label]['LABEL']
            except KeyError:
                label_name = 'unknown label'
            ogo.message('[WARNING] Expected label {} ({}) was not found in image.'.format(label,label_name))
            QualityAssurance["Expected Labels"] = False
            QualityAssurance["All Pass"] = False
    if test_labels == []:
        os.sys.exit('[ERROR] None of the expected labels were found in image.')
    labels = test_labels
    
    ogo.message('Gather each label in the image and determine number of parts.')
    
    report += '  {:>27s}\n'.format('_________________________________________________________________________________Report')
    report += '  {:>27s} {:s}\n'.format('number of labels:',str(n_labels))
    report += '  {:>27s} {:>6s} {:>10s} {:>10s} {:>19s} {:>6s}\n'.format('LABEL','PART','VOL','VOX','CENTROID','OFFSET')
    report += '  {:>27s} {:>6s} {:>10s} {:>10s} {:>19s} {:>6s}\n'.format('#','#','mm3','#','(X,Y,Z)','mm')

    # Quality check by examining femur symmetry and pelvis symmetry
    percent_femur = 3.0
    if is_smaller(filt,percent_femur,1,2): # Is the right femur smaller than left femur?
        QA_labels[1] = False
        QualityAssurance["All Pass"] = False
        ogo.message('[WARNING] Right femur is {}% smaller than left femur: FAIL'.format(percent_femur))
    if is_smaller(filt,3.0,2,1):
        QA_labels[2] = False
        QualityAssurance["All Pass"] = False
        ogo.message('[WARNING] Left femur is {}% smaller than right femur: FAIL'.format(percent_femur))
    percent_pelvis = 3.0
    if is_smaller(filt,percent_pelvis,3,4): # Is the right pelvis smaller than left pelvis?
        QA_labels[3] = False
        QualityAssurance["All Pass"] = False
        ogo.message('[WARNING] Right pelvis is {}% smaller than left pelvis: FAIL'.format(percent_pelvis))
    if is_smaller(filt,percent_pelvis,4,3):
        QA_labels[4] = False
        QualityAssurance["All Pass"] = False
        ogo.message('[WARNING] Left pelvis is {}% smaller than right pelvis: FAIL'.format(percent_pelvis))

    for idx,label in enumerate(labels):
        ogo.message('  processing label {:>2d} ({})'.format(label,lb.labels_dict[label]['LABEL']))
        desc = lb.labels_dict[label]['LABEL']
        #centroid = filt.GetCentroid(label)
        #size = filt.GetPhysicalSize(label) # volume in mm3
        
        ct_thres = ct==label
        
        ct_conn = conn.Execute(ct,ct_thres)
        ct_conn_sorted = sitk.RelabelComponent(ct_conn, sortByObjectSize=True) # could use minimumObjectSize
        
        stats.Execute(ct_conn_sorted,ct)
        n_parts = stats.GetNumberOfLabels()
        
        if (print_parts):
            label_fname = '{}_lab{:02}{}'.format(name,label,ext)
            sitk.WriteImage(ct_conn_sorted, label_fname) # ERROR: Should write image with original label (replacelabel)
        
        report += '  {:>22s}{:>5s} '.format('('+desc+')',str(label))        
        # Generate report data
        for part in stats.GetLabels(): # for each part of a given label
            centroid = stats.GetCentroid(part)
            bounding_box = stats.GetBoundingBox(part)
            bb=[0]*3
            bb[0] = bounding_box[0] + int(math.ceil(bounding_box[3]/2))
            bb[1] = bounding_box[1] + int(math.ceil(bounding_box[4]/2))
            bb[2] = bounding_box[2] + int(math.ceil(bounding_box[5]/2))

            if part==1:
                ref_centroid = centroid
                report += '{:6d} {:10.1f} {:10d} ({:5d},{:5d},{:5d})\n'.format(part,stats.GetPhysicalSize(part),stats.GetNumberOfPixels(part),bb[0],bb[1],bb[2])
            else:
                QA_labels[label] = False # any label with more than one part cannot be valid (unless it's a broken bone!)
                QualityAssurance["All Pass"] = False
                distance_between_centroids = np.sqrt((ref_centroid[0] - centroid[0])**2 + (ref_centroid[1] - centroid[1])**2 + (ref_centroid[2] - centroid[2])**2)
                report += '  {:>22s}{:>5s} {:6d} {:10.1f} {:10d} ({:5d},{:5d},{:5d}) {:6.1f}\n'\
                          .format('','',part,stats.GetPhysicalSize(part),stats.GetNumberOfPixels(part),bb[0],bb[1],bb[2],distance_between_centroids)
        
        # Generate suggestion of new label for command line suggestion
        if n_parts>1:
            for part in stats.GetLabels():
                if part>1:
                    new_label = label
                    if label == 1:
                        new_label = 2
                    if label == 2:
                        new_label = 1
                    if label == 3:
                        new_label = 4
                    if label == 4:
                        new_label = 3
                    swap = [label,part,new_label]
                    label_repair_list.append(swap)
    
    
    # Final report
    report += '\n'
    report += '  {:>27s}\n'.format('______________________________________________________________________Quality assurance')

    report += '  {:>27s}: '.format('Label')+' '.join('{:5d}'.format(label) for label,passed in QualityAssurance["Labels"].items())
    report += '\n'
    report += '  {:>27s}: '.format('Test result')+' '.join('{:>5s}'.format("PASS" if passed else "FAIL") for label,passed in QualityAssurance["Labels"].items())
    report += '\n'
    
    report += '\n'
    report += '  {:>27s}: {:>5s}\n'.format('All expected labels found',"PASS" if QualityAssurance["Expected Labels"] else "FAIL")

    report += '\n'
    report += '  {:>27s}: {:>5s}\n'.format('Final assessment',"PASS" if QualityAssurance["All Pass"] else "FAIL")

    # Command line
    cmd_line = ''
    cmd_line += '  {:>27s}\n'.format('___________________________________________________________________________Command line')
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
             
    ogo.message('Printing report')
    print('\n')
    print(report)
    report += cmd_line
    
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
    ct = sitk.ReadImage(input_image, sitk.sitkUInt8)
    
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
    
    conn = sitk.ConnectedComponentImageFilter()
    conn.SetFullyConnected(True)
    
    if relabel_parts == [] and remove_by_volume == 0:
        ogo.message('No operations requested. Suggest setting one of the follwoing:')
        ogo.message('  --relabel_parts')
        ogo.message('  --remove_by_volume')
        
    if (relabel_parts):
        ct_base = sitk.Image(ct) # Do we need to make a deep copy
        
        if np.remainder(len(relabel_parts),3) != 0:
            os.sys.exit('[ERROR] Incorrect definition for --relabel_parts. Tuples of 3 are expected.')

        n_relabel_parts = int(len(relabel_parts) / 3)
        relabel_parts = np.reshape(relabel_parts, (n_relabel_parts,3))
        
        ogo.message('Relabelling {:d} parts in image:'.format(n_relabel_parts))
        for swap in relabel_parts:
            label = int(swap[0])
            part = int(swap[1])
            new_label = int(swap[2])
                        
            ct_thres = ct==label
            ct_conn = conn.Execute(ct,ct_thres)
            ct_conn_sorted = sitk.RelabelComponent(ct_conn, sortByObjectSize=True) # could use minimumObjectSize
            
            filt.Execute(ct_conn_sorted)
            size = filt.GetPhysicalSize(part)
            
            ct_part = new_label*(ct_conn_sorted==part)
            
            ogo.message('  label {:d} ({:s}), part {} ({:.1f} mm3) --> {} ({})'\
                        .format(label,lb.labels_dict[label]['LABEL'],part,size,new_label,lb.labels_dict[new_label]['LABEL']))
            
            bin_part = ct_part>0
            mask = 1 - bin_part
            ct_base = sitk.Mask(ct_base, mask)
            ct_base = ct_base + new_label*bin_part
        
        ct_final = ct_base
        
    final_labels = get_labels(ct_final)
    ogo.message('Final image contains the following labels:')
    ogo.message('  [' + ', '.join('{:d}'.format(i) for i in final_labels) + ']')
    
    ogo.message('')
    ogo.message('Writing merged output image to file:')
    ogo.message('  {}'.format(output_image))
    sitk.WriteImage(ct_final, output_image)            
    
    # Command line
    cmd_line = ''
    cmd_line += '  {:>27s}\n'.format('___________________________________________________________________________Command line')
    cmd_line += '\n'
    cmd_line += '  {}\n'.format('ogoValidate validate \\')
    cmd_line += '    {} \n'.format(output_image)    
    
    print(cmd_line)
    
    ogo.message('Done.')
    
    
def main():
    # Setup description
    description = '''
Validates segmented output of machine learning tools of the skeleton and 
provides tools for repairing errors. The 'validate' and 'repair' options are 
separated.

Labels are used to define bones or other tissues. If a label defines a bone it 
is expected that it would be one connected component. The rationale is that all 
bones are a single connected component (unless it is a broken bone). If a label
is not single connected then it will have \'parts\'.

Validation performs the following checks:
  
  Each label is found in the image (list may be defined by user)
  Each label is a single connected part
  Each label femurs or pelvis are the same volume within a tolerance
  Each label for a bone are within an expected range (not yet implemented)
  Each bone position is as expected (not yet implemented)

The repair options include:
  Relabel voxels defined by its \'label\' and \'part\'.
  Remove any \'parts\' of a given size (not implemented)
  
NOTES:
The validation process cycles through the list of expected labels. For each 
label (e.g. Left Femur) it determined how many connected components there are 
and defines these as \'parts\'. If a bone is correctly labelled there should 
only be one part. However, it is possible that only one part is defined, but 
that the rest of the bone is labelled with an incorrect label (e.g. a portion 
of Left Femur is labelled Right Femur). Without visualizing the results of 
segmentation, this is difficult to detect. 

A bone is likely properly labelled if (a) it has only one \'part\', and (b) if 
its volume is similar to opposite bone. This works for femur and pelvis.

In cases of no symmetry (e.g. lumbar spine), it is likely properly labelled if 
(a) it has only one \'part\', and (b) the volume is within a reasonable range.
'''

    epilog = '''
Example calls: 
ogoValidate validate image.nii.gz

# This relabels label 3, part 2 to label 5
ogoValidate repair image.nii.gz image_repaired.nii.gz --relabel_parts 3 2 5

# This relabels label 3, part 2 to label 5 as well as label 3, part 3 to label 5
ogoValidate repair image.nii.gz image_repaired.nii.gz \
    --relabel_parts \
    3 2 5 \
    3 3 5


ogoValidate repair image.nii.gz --component 13 --label 0

ogoValidate validate /Users/skboyd/Desktop/ML/test/robust/frag/RETRO_01638_NNUNET.nii.gz
ogoValidate validate /Users/skboyd/Desktop/ML/test/robust/frag/RETRO_02317_NNUNET.nii.gz
ogoValidate validate /Users/skboyd/Desktop/ML/test/robust/mixed/RETRO_01964_NNUNET.nii.gz

  ogoValidate validate /Users/skboyd/Desktop/ML/data/ARTININ/seg/nnUNet/ARTININ_0002.nii.gz
  
  ogoValidate repair \
    /Users/skboyd/Desktop/ML/data/ARTININ/seg/nnUNet/ARTININ_0002.nii.gz \
    /Users/skboyd/Desktop/ML/data/ARTININ/seg/nnUNet/ARTININ_0002_REPAIR.nii.gz \
    --overwrite \
    --relabel_parts \
    1 2 3 \
    2 2 1 \
    4 2 3

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
