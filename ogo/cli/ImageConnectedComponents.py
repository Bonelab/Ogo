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
import yaml
import numpy as np
import SimpleITK as sitk
from datetime import date
from datetime import datetime
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
def ImageConnectedComponent(input_image, output_image, keep_parts, target_label, overwrite):

    # Check if output exists and should overwrite
    if output_image is not None:
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

    # Check that the target label is in the image
    if target_label == 0:
        os.sys.exit('[ERROR] Target label \"{}\" is not a valid label. Cannot select background.'.format(target_label))
    if target_label not in labels:
        os.sys.exit('[ERROR] Target label \"{}\" is not in image.'.format(target_label))
    
    # Check that the list of keep_parts is valid (postive integers)
    if keep_parts:
        ogo.message('Parts to keep:')
        ogo.message('  [' + ', '.join('{:d}'.format(i) for i in keep_parts) + ']')
        for part in keep_parts:
            if part < 1:
                os.sys.exit('[ERROR] Part numbers to keep must be a list of positive integers.')
    else:
        ogo.message('No parts to keep are defined.')
            
    filt = sitk.LabelShapeStatisticsImageFilter() # Used to get labels in image
    filt.Execute(ct)
    
    conn = sitk.ConnectedComponentImageFilter()
    conn.SetFullyConnected(True)
    stats = sitk.LabelIntensityStatisticsImageFilter()
    
    parts_list = [] # Collects the integer list of parts for the target label
    
    # Report the connected components
    report = ''
    report += '  {:>20s}\n'.format('_______________________________________________________________________Report')
    report += '  {:>20s} {:s}\n'.format('number of labels:',str(n_labels))
    report += '  {:>20s} {:>6s} {:>10s} {:>10s} {:>19s} {:>6s}\n'.format('LABEL','PART','VOL','VOX','CENTROID','OFFSET')
    report += '  {:>20s} {:>6s} {:>10s} {:>10s} {:>19s} {:>6s}\n'.format('#','#','mm3','#','(X,Y,Z)','mm')
    
    # Process the identification of labels and their parts
    for idx,label in enumerate(labels):
        try:
            desc = lb.labels_dict[label]['LABEL']
        except KeyError:
            desc = 'unknown label'
        ogo.message('  processing label {:>2d} ({})'.format(label,desc))
        
        centroid = filt.GetCentroid(label)
            
        # We are accelerating calculations by pulling an ROI for each label
        roi_bounding_box = filt.GetBoundingBox(label) # x_start, y_start, z_start, x_size, y_size, z_size
        start = (roi_bounding_box[0],roi_bounding_box[1],roi_bounding_box[2])
        size = (roi_bounding_box[3],roi_bounding_box[4],roi_bounding_box[5])
        roi_bounding_box_start = (start[0],start[1],start[2],0,0,0)
        ct_roi = sitk.RegionOfInterest(ct,size,start)
        
        ct_thres = ct_roi==label 
        ct_conn = conn.Execute(ct_roi,ct_thres) 
        ct_conn_sorted = sitk.RelabelComponent(ct_conn, sortByObjectSize=True) # could use minimumObjectSize
        stats.Execute(ct_conn_sorted,ct_roi)
        
        n_parts = stats.GetNumberOfLabels()

        # Generate report data
        desc = (desc[:11] + '..') if len(desc) > 13 else desc
        report += '  {:>15s}{:>5s} '.format('('+desc+')',str(label))        
        for part in stats.GetLabels(): # for each part of a given label
            centroid = stats.GetCentroid(part)
            bounding_box = stats.GetBoundingBox(part)
            bounding_box = tuple(map(sum, zip(bounding_box,(start[0],start[1],start[2],0,0,0)))) # We add the original start to the start of the part
            bb=[0]*3
            bb[0] = bounding_box[0] + int(math.ceil(bounding_box[3]/2))
            bb[1] = bounding_box[1] + int(math.ceil(bounding_box[4]/2))
            bb[2] = bounding_box[2] + int(math.ceil(bounding_box[5]/2))
            
            if part==1:
                ref_centroid = centroid
                report += '{:6d} {:10.1f} {:10d} ({:5d},{:5d},{:5d})\n'.format(part,stats.GetPhysicalSize(part),stats.GetNumberOfPixels(part),bb[0],bb[1],bb[2])
            else:
                distance_between_centroids = np.sqrt((ref_centroid[0] - centroid[0])**2 + (ref_centroid[1] - centroid[1])**2 + (ref_centroid[2] - centroid[2])**2)
                report += '  {:>15s}{:>5s} {:6d} {:10.1f} {:10d} ({:5d},{:5d},{:5d}) {:6.1f}\n'\
                          .format('','',part,stats.GetPhysicalSize(part),stats.GetNumberOfPixels(part),bb[0],bb[1],bb[2],distance_between_centroids)
            
            # Collect the full list of part numbers
            if label == target_label:
                parts_list.append(part)
                
    report += '  {:>20s}\n'.format('_______________________________________________________________________Report')
    
    ogo.message('Printing report')
    print('\n')
    print(report)
    
    # ---------------------------------------------------    
    # Removing parts by relabelling them to background

    if (not keep_parts) or (not output_image):
        ogo.message('No list of parts to keep is defined.')
        ogo.message('No output image filename is defined.')
        ogo.message('Done.')
        
    label = target_label
    
    for part in keep_parts:
        if part not in parts_list:
            os.sys.exit('[ERROR] Part {:d} is not one of the {:d} parts for label {:d}.'.format(part,len(parts_list),label))
    
    # Start with a copy of the original CT. We'll erase labels we don't want.
    ct_base = sitk.Image(ct)
    
    # Create a list of parts
    ct_thres = ct==label
    ct_conn = conn.Execute(ct,ct_thres)
    ct_conn_sorted = sitk.RelabelComponent(ct_conn, sortByObjectSize=True) # could use minimumObjectSize
    filt.Execute(ct_conn_sorted)

    # Erase all voxels of the label. In the next section we'll add the parts to keep back in.
    bin_part = ct==label
    mask = 1 - bin_part
    ct_base = sitk.Mask(ct_base, mask)
    ct_base = ct_base + 0*bin_part
    
    try:
        desc_label = lb.labels_dict[label]['LABEL']
    except KeyError:
        desc_label = 'unknown label'
    
    new_label = 0
    try:
        desc_new_label = lb.labels_dict[new_label]['LABEL']
    except KeyError:
        desc_new_label = 'unknown label'
    
    ogo.message('Label {:d} has {:d} parts: keeping {:d} parts and removing {:d} parts.'.format(label,len(parts_list),len(keep_parts),len(parts_list)-len(keep_parts)))

    for part in parts_list:
        
        if part in keep_parts:
            size = filt.GetPhysicalSize(part)
            ct_part = (ct_conn_sorted==part)
        
            ogo.message('  Keeping label {:d} ({:s}), part {} ({:.1f} mm3)'.format(\
                           label, desc_label, part, size))
            
            bin_part = ct_part>0
            mask = 1 - bin_part
            ct_base = sitk.Mask(ct_base, mask)
            ct_base = ct_base + label*bin_part
    
    ogo.message('')
    ogo.message('Writing output image to file:')
    ogo.message('  {}'.format(output_image))
    sitk.WriteImage(ct_base, output_image)            
    
    ogo.message('Done.')
    
    
def main():
    # Setup description
    description = '''
This tools takes a mask containing 1 or more labels as input and conducts
a connected component evaluation of the selected label. The components are
called parts that are ordered from largest to smallest. 

Define a list of parts to keep. The remaining parts will be set to 0 label.

Only the target label is affected. Other labels pass through to the output.

'''

    epilog = '''
Example calls: 
ogoImageConnectedComponent mask.nii.gz

'''
    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoImageConnectedComponent",
        description=description,
        epilog=epilog
    )
    parser.add_argument('input_image', help='Input image file (*.nii, *.nii.gz)')
    parser.add_argument('--output_image', metavar='FILE', default=None, help='Output image file (*.nii, *.nii.gz)')
    parser.add_argument('--target_label', type=int, default=1, metavar='LABEL', help='Label for connected component analysis (default: %(default)s)')
    parser.add_argument('--keep_parts', type=int, nargs='*', default=[], metavar='LABEL', help='Part numbers to keep (default: %(default)s)')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite output image without asking')
    
    # Parse and display
    args = parser.parse_args()
    print(echo_arguments('ImageConnectedComponent', vars(args)))

    # Run program
    ImageConnectedComponent(**vars(args))

if __name__ == '__main__':
    main()
