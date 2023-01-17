# /------------------------------------------------------------------------------+
# | 27-NOV-2022                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+

#Imports
import argparse
import os
import sys
import time
import math
import vtk
import vtkbone
import numpy as np
import SimpleITK as sitk
import ogo.dat.OgoMasterLabels as lb

from ogo.util.echo_arguments import echo_arguments
import ogo.util.Helper as ogo
from vtk.util.numpy_support import vtk_to_numpy
    
def resolve_totalsegmentator_label_and_id(filename):
    labels = [v['LABEL'] for k, v in lb.labels_dict.items() if ((v['TOTSEG']+'.nii.gz') in filename)]
    if labels:
        keys = [k for k, v in lb.labels_dict.items() if v['LABEL'] == labels[0]]
        if keys:
            return labels[0],keys[0]
    
    return None,None # return None,None if cannot be resolved
    
def get_labels(ct):
    filt = sitk.LabelShapeStatisticsImageFilter()
    filt.Execute(ct)
    labels = filt.GetLabels()
    return labels
    
def merge_labels(input_filenames, output_filename, add_multilabel, merge_method, collection, overwrite=False):
        
    # Check if output exists and should overwrite
    if os.path.isfile(output_filename) and not overwrite:
        result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(output_filename))
        if result.lower() not in ['y', 'yes']:
            ogo.message('Not overwriting. Exiting...')
            os.sys.exit()
            
    # Check if output file is in correct format
    if not (output_filename.lower().endswith('.nii') or output_filename.lower().endswith('.nii.gz')):
        os.sys.exit('[ERROR] Output must be type NIFTI file: \"{}\"'.format(output_filename))
    
    # Check if input files exist and are correct format
    ogo.message('The following {} images are available for merging:'.format(len(input_filenames)))
    max_number_input_files_to_print = 15
    for idx,input_filename in enumerate(input_filenames):
        if not os.path.isfile(input_filename):
            os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(input_filename))
        if not (input_filename.lower().endswith('.nii') or input_filename.lower().endswith('.nii.gz')):
            os.sys.exit('[ERROR] Input must be type NIFTI file: \"{}\"'.format(input_filename))
        ogo.message('  {}'.format(input_filename))
        if idx>max_number_input_files_to_print:
            ogo.message('  ...not showing remaining {} input files...'.format(len(input_filenames)-idx-1))
            ogo.message('')
            break
        
    if len(input_filenames)<2:
        os.sys.exit('[ERROR] At least two input files are required.')
        
    # Establish a valid list of images for collection
    if collection == 'all':
        valid_list = [v['LABEL'] for k, v in lb.labels_dict.items()]
    elif collection == 'ossai':
        valid_list = ['Femur Right', 'Femur Left', 'Pelvis Right', 'Pelvis Left', 'Sacrum', 'L6', 'L5', 'L4', 'L3', 'L2', 'L1']
    elif collection == 'skeleton':
        valid_list = ['Femur Right', 'Femur Left', 'Pelvis Right', 'Pelvis Left', 'Sacrum', 'L6', 'L5', 'L4', 'L3', 'L2', 'L1', 'T12', 'T11', 'T10', 'T9', 'T8', 'Humerus Right', 'Humerus Left', 'T7', 'T6', 'T5', 'T4', 'T3', 'T2', 'T1', 'C7', 'C6', 'C5', 'C4', 'C3', 'C2', 'C1', 'Rib Left 1', 'Rib Left 2', 'Rib Left 3', 'Rib Left 4', 'Rib Left 5', 'Rib Left 6', 'Rib Left 7', 'Rib Left 8', 'Rib Left 9', 'Rib Left 10', 'Rib Left 11', 'Rib Left 12', 'Rib Right 1', 'Rib Right 2', 'Rib Right 3', 'Rib Right 4', 'Rib Right 5', 'Rib Right 6', 'Rib Right 7', 'Rib Right 8', 'Rib Right 9', 'Rib Right 10', 'Rib Right 11', 'Rib Right 12', 'Scapula Left', 'Scapula Right', 'Clavicula Left', 'Clavicula Right']
    elif collection == 'cardio':
        valid_list = ['Heart Myocardium', 'Heart Atrium Left', 'Heart Ventricle Left', 'Heart Atrium Right', 'Heart Ventricle Right', 'Aorta', 'Inferior Vena Cava', 'Portal Vein and Splenic Vein', 'Pulmonary Artery', 'Iliac Artery Left', 'Iliac Artery Right', 'Iliac Vena Left', 'Iliac Vena Right']
    elif collection == 'organs':
        valid_list = ['Face', 'Brain', 'Trachea', 'Lung Upper Lobe Left', 'Lung Lower Lobe Left', 'Lung Upper Lobe Right', 'Lung Middle Lobe Right', 'Lung Lower Lobe Right', 'Adrenal Gland Right', 'Adrenal Gland Left', 'Spleen', 'Kidney Right', 'Kidney Left', 'Gallbladder', 'Liver', 'Pancreas']
    elif collection == 'gastro':
        valid_list = ['Esophagus', 'Stomach', 'Duodenum', 'Small Bowel', 'Colon', 'Urinary Bladder']
    elif collection == 'muscle':
        valid_list = ['Autochthon Left', 'Autochthon Right', 'Iliopsoas Left', 'Iliopsoas right', 'Gluteus Maximus Left', 'Gluteus Maximus Right', 'Gluteus Medius Left', 'Gluteus Medius Right', 'Gluteus Minimus Left', 'Gluteus Minimus Right']
    else:
        os.sys.exit('[ERROR] Unknown collection defined: {}'.format(collection))
    
    # We will use the sitk.MergeLabelMapFilter, but that only accepts LabelMap (which is a special data type).
    # We need to convert each image we read to LabelMaps using sitk.LabelImageToLabelMapFilter. Then after we 
    # complete all the merge operations we convert it back to LabelImage.
    labelmapfilt = sitk.LabelImageToLabelMapFilter()
    labelimagefilt = sitk.LabelMapToLabelImageFilter()
    changefilt = sitk.ChangeLabelImageFilter()
    mergefilt = sitk.MergeLabelMapFilter()
    successful_merge = False
    
    # Set type of merge
    ogo.message('Merge method is \'' + merge_method + '\'')
    if merge_method == 'aggregate':
        mergefilt.Aggregate
    elif merge_method == 'strict':
        mergefilt.Strict
    elif merge_method == 'keep':
        mergefilt.Keep
    elif merge_method == 'pack':
        mergefilt.Pack
    else:
        os.sys.exit('[ERROR] Unknown merge method.')
    
    # This option adds multi-label images together
    if add_multilabel:
        ogo.message('Adding multi-label images together.')
        if len(input_filenames)>2:
            ogo.message('[WARNING] Only first 2 of {} input images will be added together.'.format(len(input_filenames)))
            ogo.message('          The remaining input images are ignored.')

        ct1 = sitk.ReadImage(input_filenames[0], sitk.sitkUInt8)
        ct2 = sitk.ReadImage(input_filenames[1], sitk.sitkUInt8)
        ct1_labels = get_labels(ct1)
        ct2_labels = get_labels(ct2)
        ogo.message('  Image 1:')
        ogo.message('  {}'.format(input_filenames[0]))
        ogo.message('  labels are [' + ', '.join('{:d}'.format(i) for i in ct1_labels) + ']')
        ogo.message('  Image 2:')
        ogo.message('  {}'.format(input_filenames[1]))
        ogo.message('  labels are [' + ', '.join('{:d}'.format(i) for i in ct2_labels) + ']')
        ct1 = labelmapfilt.Execute(ct1)
        ct2 = labelmapfilt.Execute(ct2)
        ct_final = mergefilt.Execute(ct1, ct2)
        successful_merge = True

    else:
        # Cycle through all input images that are part of the collection. 
        ogo.message('Collection is \'' + collection + '\'')
        n_valid_images = 0
        first_image = True
        for fn in input_filenames:
            label_name,label_id = resolve_totalsegmentator_label_and_id(fn)
            if label_name in valid_list: # found a valid image
                ct = sitk.ReadImage(fn, sitk.sitkUInt8)
                input_labels = get_labels(ct)
                if len(input_labels)==1:
                    n_valid_images += 1
                    ogo.message('  {:>23s} (changing label {:3d} to {:3d})'.format(label_name,input_labels[0],label_id))
                    map = {input_labels[0]: label_id}
                    changefilt.SetChangeMap(map)
                    ct = changefilt.Execute(ct)
                    ct = labelmapfilt.Execute(ct)
                    
                    if first_image:
                        ct_final = ct
                        first_image = False
                    else:
                        ct_final = mergefilt.Execute(ct_final, ct)
                        successful_merge = True
                else:
                    ogo.message('[WARNING] Image with {} labels cannot be merged.'.format(len(input_labels)))
                    ogo.message('  {}'.format(fn))
            #else:
            #    ogo.message('[WARNING] Input file is not part of any collection. Most')
            #    ogo.message('          likely it is a multi-label image and should be')
            #    ogo.message('          combined with other images using --add_multilabel')
            #    ogo.message('          option.')
            #    ogo.message('  {}'.format(fn))
               
        if n_valid_images != len(valid_list):
            ogo.message('')
            ogo.message('[WARNING]: Found only {} of total possible {} images for collection \'{}\''.format(n_valid_images,len(valid_list),collection))
        
    if not successful_merge:
        os.sys.exit('[ERROR] Merging of listed input files is unsuccessful. Try --add_multilabel option?')
         
    ct_final = labelimagefilt.Execute(ct_final)
    ogo.message('')
    ogo.message('Final merged image contains the following labels:')
    final_labels = get_labels(ct_final)
    ogo.message('  [' + ', '.join('{:d}'.format(i) for i in final_labels) + ']')
            
    ogo.message('')
    ogo.message('Writing merged output image to file:')
    ogo.message('  {}'.format(output_filename))
    sitk.WriteImage(ct_final, output_filename)
        
    ogo.message('')
    ogo.message('Done.')

def main():
    # Setup description
    description='''
A utility to merge labels in an image and ensuring correct labels for each 
component. 

This is a particularly useful utility if using TotalSegmentator
https://arxiv.org/abs/2208.05868

All merged images are expected to have the same number of dimensions.

Merge methods are:
    aggregate: If the same label is found several times in the label maps, the
               label objects with the same label are merged.
    strict:    Keeps the labels unchanged and raises an exception if the same
               label is found in several images.
    keep:      Does its best to keep the label unchanged, but if a label is
               already used in a previous label map, a new label is assigned.
    pack:      Relabel all the label objects by order of processing.

'''
    epilog = '''
Example calls: 
  ogoMergeLabels femur_right.nii.gz femur_left.nii.gz femurs.nii.gz
  ogoMergeLabels *.nii.gz skeleton.nii.gz --collection skeleton

'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="MergeLabels",
        description=description,
        epilog=epilog
    )
    parser.add_argument('input_filenames', nargs='*', metavar='FILES IN', help='Input files')
    parser.add_argument('output_filename', metavar='FILE OUT', help='Output NIFTI file (default: %(default)s)')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite output without asking')
    parser.add_argument('--add_multilabel', action='store_true', help='Adds multi-label images')
    parser.add_argument('--merge_method', default='strict',choices=['aggregate', 'strict', 'keep', 'pack'],
                                                           help='Select merge rules (default: %(default)s)')
    parser.add_argument('--collection', default='all', choices=['all', 'ossai', 'skeleton', 'cardio', 'organs', 'muscle', 'gastro'],
                                                           help='Select collection (default: %(default)s)')

    # Parse and display
    args = parser.parse_args()
    #print(echo_arguments('MergeLabels', vars(args)))
        
    # Run program
    merge_labels(**vars(args))

if __name__ == '__main__':
    main()
