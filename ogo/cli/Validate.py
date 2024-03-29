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
from scipy.spatial import procrustes
from datetime import date
from datetime import datetime
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from ogo.util.echo_arguments import echo_arguments
import ogo.util.Helper as ogo
import ogo.dat.OgoMasterLabels as lb

def expected_bone_volumes():
    bone_volumes = {
         1:{ 'Min':50000.0, 'Max':1000000.0,  'Stdev':39632.3,  'Short femur':153110.0}, # Femur Right
         2:{ 'Min':50000.0, 'Max':1000000.0,  'Stdev':39494.4,  'Short femur':153110.0}, # Femur Left
         3:{ 'Min':80000.0, 'Max':600000.0,  'Stdev':70518.6},                           # Pelvis Right
         4:{ 'Min':80000.0, 'Max':600000.0,  'Stdev':70777.9},                           # Pelvis Left
         5:{'Min':100000.0, 'Max':390000.0,  'Stdev':43914.4},                           # Sacrum
         6:{ 'Min':30000.0, 'Max':110000.0,  'Stdev':13011.0},                           # L5
         7:{ 'Min':30000.0, 'Max':110000.0,  'Stdev':13444.9},                           # L4
         8:{ 'Min':30000.0, 'Max':110000.0,  'Stdev':13610.6},                           # L3
         9:{ 'Min':30000.0, 'Max':110000.0,  'Stdev':12253.7},                           # L2
        10:{ 'Min':30000.0, 'Max':110000.0,  'Stdev':11756.1},                           # L1
        11:{ 'Min':30000.0, 'Max':110000.0,  'Stdev':11756.1}                            # L6
    }
    return bone_volumes

def expected_Procrustes():
    # Used RETRO_00071.nii.gz as reference (top of femurs)
    # Labels 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
    positions = \
    [( -85.02785, -142.38381, 190.00000 ), \
     ( 134.24816, -147.52873, 190.00000 ), \
     ( -47.02655, -129.48297, 208.27110 ), \
     (  97.00727, -133.80435, 211.17284 ), \
     (  26.25469,  -84.00061, 237.57798 ), \
     (  23.43566, -121.04673, 287.29457 ), \
     (  25.05111, -127.58523, 319.54591 ), \
     (  26.85287, -123.60419, 355.15329 ), \
     (  25.53003, -116.45560, 390.07599 ), \
     (  24.17049, -109.25177, 425.37350 )]
    return positions
    
def expected_Procrustes_with_L6():
    # Used RETRO_01455.nii.gz as reference (top of femurs)
    # Labels 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
    positions = \
    [( -99.221,    -13.561,    147.75 ), \
     ( 125.225,     -1.464,    135.88 ), \
     ( -57.060,     -5.412,   176.046 ), \
     (  92.237,     -1.380,   165.237 ), \
     (  18.095,     40.348,   211.248 ), \
     (  24.245,     -1.285,   277.174 ), \
     (  25.745,      0.369,   311.009 ), \
     (  24.853,      6.145,   344.486 ), \
     (  24.297,     13.684,   376.100 ), \
     (  24.441,     19.705,   406.925 ), \
     (  21.298,      6.931,   246.802 )]
    return positions

def print_centroids(centroids):
    print('      centroids: {}'.format(len(centroids)))
    for c in centroids:
        #phys = [x - y for x, y in zip(c, origin)]
        #pos = [x / y for x, y in zip(phys, el_size_mm)]
        print('      [' + ', '.join('{:8.2f}'.format(i) for i in c) + ']')
        #print('      [' + ', '.join('{:8.2f}'.format(i) for i in c) + '], ['+ ', '.join('{:8d}'.format(int(i)) for i in pos) + ']')

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
    average_size = (size1+size2)/2.0
    thres = average_size * percent_threshold / 100.0
    
    if (size2 - size1) > thres:
        return True
    else:
        return False

def special_cases(n_parts,label,stats):
    
    if label==5: # Sacrum – sometimes the tail is separate from the rest
        if n_parts==2 and stats.GetPhysicalSize(2) > 400:
            return True
    
    elif label==6: # L5 – sometimes the pedicules are separate from the vertebral body (Pars defect)
        if n_parts==2 and stats.GetPhysicalSize(2) > 2000:
            return True

    else:
        return False

# +------------------------------------------------------------------------------+
def validate(input_image, report_file, yaml_file, expected_labels, overwrite, func):
        
    # Quality assurance dictionary
    QA_dict = {}
    
    # Collection information into dictionary for YAML file
    info_dict = {}
    info_dict['runtime']={}
    info_dict['runtime']['date']=str(date.today())
    info_dict['runtime']['time']=datetime.now().strftime("%H:%M:%S")
    info_dict['runtime']['script']=os.path.splitext(os.path.basename(sys.argv[0]))[0]
    info_dict['runtime']['version']=script_version
    info_dict['input_parameters']={}
    info_dict['input_parameters']['input_image']= input_image
    info_dict['input_parameters']['report_file']= report_file
    info_dict['input_parameters']['yaml_file']= yaml_file
    info_dict['input_parameters']['expected_labels']= expected_labels
    info_dict['input_parameters']['overwrite']= overwrite
    
    # Check if output exists and should overwrite
    if not report_file is None:
        if os.path.isfile(report_file) and not overwrite:
            result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(report_file))
            if result.lower() not in ['y', 'yes']:
                ogo.message('Not overwriting. Exiting...')
                os.sys.exit()

    # Check if output exists and should overwrite
    if not yaml_file is None:
        if os.path.isfile(yaml_file) and not overwrite:
            result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(yaml_file))
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
    
    el_size_mm = ct.GetSpacing()
    origin = ct.GetOrigin()

    # Create base filename (and some extra information)
    basename = os.path.basename(input_image)
    name, ext = os.path.splitext(input_image)
    if 'gz' in ext:
        name = os.path.splitext(name)[0]  # Manages files with double extension
        ext = '.nii' + ext
    dirname = os.path.dirname(input_image)
    
    # Start the report
    report = ''
    report += '  {:>20s}\n'.format('___________________________________________________________________Validation')
    report += '  {:>20s} {:s}\n'.format('ID:', os.path.basename(input_image))
    report += '  {:>20s} {:s}\n'.format('python script:', os.path.splitext(os.path.basename(sys.argv[0]))[0])
    report += '  {:>20s} {:.2f}\n'.format('version:', script_version)
    report += '  {:>20s} {:s}\n'.format('creation date:', str(date.today()))
    report += '\n'
    
    info_dict['fileinfo']={}
    info_dict['fileinfo']['basename']=basename
    info_dict['fileinfo']['name']=name
    info_dict['fileinfo']['extension']=ext
    info_dict['fileinfo']['dir']=dirname
    
    extra_msg = ''
          
    # Gather all labels in image and determine number of parts each label is broken into
    filt = sitk.LabelShapeStatisticsImageFilter() # Used to get labels in image
    filt.Execute(ct)
    
    n_labels = filt.GetNumberOfLabels()
    labels = filt.GetLabels()
    hasL6 = filt.HasLabel(11)
    label_repair_list = [] # list of [label, part, new_label]
    
    conn = sitk.ConnectedComponentImageFilter()
    conn.SetFullyConnected(True)
    stats = sitk.LabelIntensityStatisticsImageFilter()
    
    # Check if expected labels are in image
    test_labels = []
    testWereExpectedLabelsFound = True
    testDidOverallQualityAssessmentPass = True
    for label in expected_labels:
        try:
            desc = lb.labels_dict[label]['LABEL']
        except KeyError:
            desc = 'unknown label'
        
        QA_dict[label] = {}
        QA_dict[label]['name'] = desc
        QA_dict[label]['label'] = label

        if label in labels:
            test_labels.append(label)
        else:
            if label != 11:
                ogo.message('[WARNING] Expected label {} ({}) was not found in image.'.format(label,desc))
                testWereExpectedLabelsFound = False
                testDidOverallQualityAssessmentPass = False
            else: # we have a soft constraint on L6
                ogo.message('[NOTICE] Label {} ({}) was not found, but is not really expected.'.format(label,desc))

    if test_labels == []:
        os.sys.exit('[ERROR] None of the expected labels were found in image.')
        
    labels = test_labels

    ogo.message('Examine each label to determine number of parts and volume.')
    
    report += '  {:>20s}\n'.format('_______________________________________________________________________Report')
    report += '  {:>20s} {:s}\n'.format('number of labels:',str(n_labels))
    report += '  {:>20s} {:>6s} {:>10s} {:>10s} {:>19s} {:>6s}\n'.format('LABEL','PART','VOL','VOX','CENTROID','OFFSET')
    report += '  {:>20s} {:>6s} {:>10s} {:>10s} {:>19s} {:>6s}\n'.format('#','#','mm3','#','(X,Y,Z)','mm')
    
    current_centroids = []
    
    # Process the identification of labels and their parts
    for idx,label in enumerate(expected_labels):
        try:
            desc = lb.labels_dict[label]['LABEL']
        except KeyError:
            desc = 'unknown label'
        ogo.message('  processing label {:>2d} ({})'.format(label,desc))
        
        if label not in labels: # although we expected this label, there is none in the image (e.g. missing vertebrae, femurs, etc)
            QA_dict[label]['n_parts'] = 0
            QA_dict[label]['total_volume'] = 0
            QA_dict[label]['errors'] = ''
            QA_dict[label]['warnings'] = ''
            QA_dict[label]['status'] = True
            
            # Generate report data
            report += '  {:>15s}{:>5s} '.format('('+desc+')',str(label))        
            report += '{:>6s} {:>10s} {:>10s} ({:>5s},{:>5s},{:>5s})\n'.format('-','-','-','-','-','-')
            
        else: # label is in the image

            centroid = filt.GetCentroid(label)
            
            if label==1 or label==2: # For left and right femur, we take the top Z value instead of centroid for use with Procrustes
                bounding_box = filt.GetBoundingBox(label)
                height = (bounding_box[2]+bounding_box[5])*el_size_mm[2]+origin[2]
                lst = list(centroid) # to change a tuple we convert to list, edit, then change back
                lst[2] = height
                centroid = tuple(lst)
                
                #pos = [x / y for x, y in zip(centroid, el_size_mm)]
                #print('.  el_size_mm [' + ', '.join('{:8.2f}'.format(i) for i in el_size_mm) + ']')
                #print('.    centroid [' + ', '.join('{:8.2f}'.format(i) for i in centroid) + ']')
                #print('          pos [' + ', '.join('{:8d}'.format(int(i)) for i in pos) + ']')
                #print(' bounding_box [' + ', '.join('{:8d}'.format(int(i)) for i in bounding_box) + ']')
                
            current_centroids.append(centroid) # current_centroids is used for Procrustes
            
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
            QA_dict[label]['n_parts'] = n_parts
            QA_dict[label]['total_volume'] = 0
            QA_dict[label]['errors'] = ''
            QA_dict[label]['warnings'] = ''
            QA_dict[label]['status'] = True
            
            # Set some warning or error codes if necessary
            if (n_parts>1):
                if special_cases(n_parts,label,stats):
                    QA_dict[label]['warnings'] += 'wSPEC ' # Sacrum and L5 are sometimes two parts
                    msg = '[WARNING] Label {:d} ({}) exception with {} parts ({:.1f} mm3)'.format(label,desc,n_parts,stats.GetPhysicalSize(2))
                    ogo.message(msg)
                    extra_msg += '  ' + msg + '\n'
                else:
                    QA_dict[label]['errors'] += 'eFRAG ' # More than one part indicates the bone is fragmented
                    QA_dict[label]['status'] = False
                
            # Generate report data
            report += '  {:>15s}{:>5s} '.format('('+desc+')',str(label))        
            for part in stats.GetLabels(): # for each part of a given label
                centroid = stats.GetCentroid(part)
                bounding_box = stats.GetBoundingBox(part)
                bounding_box = tuple(map(sum, zip(bounding_box,(start[0],start[1],start[2],0,0,0)))) # We add the original start to the start of the part
                bb=[0]*3
                bb[0] = bounding_box[0] + int(math.ceil(bounding_box[3]/2))
                bb[1] = bounding_box[1] + int(math.ceil(bounding_box[4]/2))
                bb[2] = bounding_box[2] + int(math.ceil(bounding_box[5]/2))
                
                QA_dict[label]['total_volume'] += stats.GetPhysicalSize(part)
                
                if part==1:
                    ref_centroid = centroid
                    report += '{:6d} {:10.1f} {:10d} ({:5d},{:5d},{:5d})\n'.format(part,stats.GetPhysicalSize(part),stats.GetNumberOfPixels(part),bb[0],bb[1],bb[2])
                else:
                    distance_between_centroids = np.sqrt((ref_centroid[0] - centroid[0])**2 + (ref_centroid[1] - centroid[1])**2 + (ref_centroid[2] - centroid[2])**2)
                    report += '  {:>15s}{:>5s} {:6d} {:10.1f} {:10d} ({:5d},{:5d},{:5d}) {:6.1f}\n'\
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

    # Quality checks -------------------------------------------------------------------
    # 1. Examining femur symmetry
    if ((1 in labels) and (2 in labels)):
        ogo.message('Examining femur symmetry.')
        percent_femur = 10.0
        if is_smaller(filt,percent_femur,1,2): # Is the right femur smaller than left femur?
            QA_dict[1]['warnings'] += 'wSYMM '
            msg = '[WARNING] Right femur is more than {}% smaller than left femur.'.format(percent_femur)
            ogo.message(msg)
            extra_msg += '  ' + msg + '\n'
        if is_smaller(filt,percent_femur,2,1):
            QA_dict[2]['warnings'] += 'wSYMM '
            msg = '[WARNING] Left femur is more than {}% smaller than right femur.'.format(percent_femur)
            ogo.message(msg)
            extra_msg += '  ' + msg + '\n'
    
    # 2. Examining pelvis symmetry
    if ((3 in labels) and (4 in labels)):
        ogo.message('Examining pelvis symmetry.')
        percent_pelvis = 10.0
        if is_smaller(filt,percent_pelvis,3,4): # Is the right pelvis smaller than left pelvis?
            QA_dict[3]['warnings'] += 'wSYMM '
            msg = '[WARNING] Right pelvis is more than {}% smaller than left pelvis.'.format(percent_pelvis)
            ogo.message(msg)
            extra_msg += '  ' + msg + '\n'
        if is_smaller(filt,percent_pelvis,4,3):
            QA_dict[4]['warnings'] += 'wSYMM '
            msg = '[WARNING] Left pelvis is more than {}% smaller than right pelvis.'.format(percent_pelvis)
            ogo.message(msg)
            extra_msg += '  ' + msg + '\n'
    
    # 3. Check bone sizes are in normal range
    ogo.message('Examining bone volumes are within expected range.')
    bone_sizes_dict = expected_bone_volumes()
    for label in labels:
        if bone_sizes_dict.get(label) is not None:
            volume = filt.GetPhysicalSize(label)
            min_volume = bone_sizes_dict[label]['Min']
            max_volume = bone_sizes_dict[label]['Max']
            try:
                desc = lb.labels_dict[label]['LABEL']
            except KeyError:
                desc = 'unknown label'
            
            if volume < min_volume:
                QA_dict[label]['errors'] += 'eMINVOL '
                QA_dict[label]['status'] = False
                msg = '[ERROR] Label {:d} ({}) volume is below minimum ({} mm3).'.format(label,desc,min_volume)
                ogo.message(msg)
                extra_msg += '  ' + msg + '\n'
            
            if volume > max_volume:
                QA_dict[label]['errors'] += 'eMAXVOL '
                QA_dict[label]['status'] = False
                msg = '[ERROR] Label {:d} ({}) volume is above maximum ({} mm3).'.format(label,desc,max_volume)
                ogo.message(msg)
                extra_msg += '  ' + msg + '\n'
            
            
            if label == 1 or label == 2:
                short_femur_volume = bone_sizes_dict[label]['Short femur']
                if volume < short_femur_volume:
                    QA_dict[label]['warnings'] += 'wNOFEA '
                    msg = '[WARNING] FEA of label {:d} ({}) not recommended ({} mm3).'.format(label,desc,short_femur_volume)
                    ogo.message(msg)
                    extra_msg += '  ' + msg + '\n'
    
    # 4. Check relative positions by Procrustes analysis
    ogo.message('Examining bone positions by Procrustes analysis.')
    testDidProcrustesAnalysisPass = True
    sensitivity = 0.008 # adjust with smaller number to make more sensitive
    if hasL6:
        expected_centroids = expected_Procrustes_with_L6()
    else:
        expected_centroids = expected_Procrustes()

    #print('Expected:')
    #print_centroids(expected_centroids)
    #print('Current:')
    #print_centroids(current_centroids)
    
    disparity = 1.0
    
    if len(current_centroids)==len(expected_centroids):
        mtx1, mtx2, disparity = procrustes(expected_centroids, current_centroids)
        if disparity > sensitivity:
            msg = '[WARNING] Procrustes analysis detects position mismatch ({:.3f} > {:.3f}).'.format(disparity,sensitivity)
            ogo.message(msg)
            extra_msg += '  ' + msg + '\n'
            testDidProcrustesAnalysisPass = False
            testDidOverallQualityAssessmentPass = False

    line_hdr =   'line_{},{},'.format('hdr','name')
    line_units = 'line_{},{},'.format('units','[text]')
    line_data =  'line_{},{},'.format('data',basename)

    # Final report
    report += '\n'
    report += '  {:>20s}\n'.format('______________________________________________________Quality Assurance Check')
    report += '  {:>20s} {:>6s} {:>10s} {:>10s} {:>23s}\n'.format('LABEL','NPARTS','TOTAL_VOL','STATUS','ERRORS and WARNINGS')
    report += '  {:>20s} {:>6s} {:>10s} {:>10s} {:>23s}\n'.format('#','#','mm3','','')
    
    testAreAllLabelsIntact=True
    testIsRightFemurIntact=True
    testIsLeftFemurIntact=True
    testIsL1toL4Intact=True
    
    for label in expected_labels:
        
        desc = QA_dict[label]['name']
        n_parts = QA_dict[label]['n_parts']
        total_volume = QA_dict[label]['total_volume']
        warnings = QA_dict[label]['warnings']
        errors = QA_dict[label]['errors']
        status = QA_dict[label]['status']
        
        if n_parts == 0:
            report += '  {:>15s}{:>5s} {:>6s} {:>10s} {:>10s} {:>23s}\n'.format('('+desc+')',str(label),'-','-','-','-')        
        else:
            report += '  {:>15s}{:>5s} {:6d} {:10.1f} {:>10s} {:>23s}\n'.format('('+desc+')',str(label),n_parts,total_volume,str(status),warnings+errors)                    
        
        line_hdr += '{},{},{},{},{},{},{},'.format('desc','label','n_parts','total_vol','warnings','errors','status')
        line_units += '{},{},{},{},{},{},{},'.format('[text]','[#]','[#]','[mm3]','[msg]','[msg]','[test]')
        if n_parts == 0:
            line_data += '{},{},{},{},{},{},{},'.format(desc,label,'','','','','')
        else:
            line_data += '{},{},{},{:.1f},{},{},{},'.format(desc,label,n_parts,total_volume,warnings,errors,status)
        
        if status is False:
            testAreAllLabelsIntact=False
            testDidOverallQualityAssessmentPass = False
            if ((label >= 7) and (label <= 10)):
                testIsL1toL4Intact=False
            if label == 1:
                testIsRightFemurIntact=False
            if label == 2:
                testIsLeftFemurIntact=False
            
    report += '\n'
    pa_string = 'Pass procrustes analysis ({:.3f} < {:.3f})'.format(disparity,sensitivity)
    report += '  {:>42s}? {:>5s}\n'.format(pa_string,str(testDidProcrustesAnalysisPass))
    report += '\n'
    report += '  {:>42s}? {:>5s}\n'.format('Pass all labels found',str(testWereExpectedLabelsFound))
    report += '\n'
    report += '  {:>42s}? {:>5s}\n'.format('Pass right femur label intact',str(testIsRightFemurIntact))
    report += '  {:>42s}? {:>5s}\n'.format('Pass left femur label intact',str(testIsLeftFemurIntact))
    report += '  {:>42s}? {:>5s}\n'.format('Pass spine labels intact',str(testIsL1toL4Intact))
    report += '  {:>42s}? {:>5s}\n'.format('Pass all labels intact',str(testAreAllLabelsIntact))
    report += '\n'
    report += '  {:>42s}? {:>5s}\n'.format('Pass overall',str(testDidOverallQualityAssessmentPass))
    
    line_hdr += '{},{},{},{},{}\n'.format('disparity','procrustes','labels_found','intact_right_femur','intact_left_femur','intact_spine','intact_all_labels','final')
    line_units += '{},{},{},{},{}\n'.format('[1.0]','[test]','[test]','[test]','[test]','[test]','[test]','[test]')
    line_data += '{:.3f},{},{},{},{}\n'.format(disparity,str(testDidProcrustesAnalysisPass),\
                                                         str(testWereExpectedLabelsFound),\
                                                         str(testIsRightFemurIntact),\
                                                         str(testIsLeftFemurIntact),\
                                                         str(testIsL1toL4Intact),\
                                                         str(testAreAllLabelsIntact),\
                                                         str(testDidOverallQualityAssessmentPass))
    
    info_dict['check']={}
    info_dict['check']['labels']=QA_dict
    info_dict['check']['procrustes']={'n_expected': len(expected_centroids), 'n_current': len(current_centroids), 'disparity': float(disparity), 'sensitivity': sensitivity, 'status': testDidProcrustesAnalysisPass}
    info_dict['check']['all_found']={'n_labels': len(expected_labels), 'status': testWereExpectedLabelsFound}
    info_dict['check']['intact_right_femur']={'status': testIsRightFemurIntact}
    info_dict['check']['intact_left_femur']={'status': testIsLeftFemurIntact}
    info_dict['check']['intact_spine']={'status': testIsL1toL4Intact}
    info_dict['check']['intact_all_labels']={'status': testAreAllLabelsIntact}
    info_dict['check']['final']={'status': testDidOverallQualityAssessmentPass}
    
    report += '\n'
    report += '  Legend for errors and warnings:\n'
    report += '    eFRAG   – Error due to multiple connected components resulting in fragments\n'
    report += '    eMAXVOL – Error due to exceeding the expected maximum bone volume\n'
    report += '    eMINVOL – Error due to being lower than the expected minimum bone volume\n'
    report += '    wSPEC   – Warning that L5 or sacrum may have a disconnected component\n'
    report += '    wNOFEA  – Warning that the femurs are cut too short for FE analysis\n'
    
    # For grep in report
    grep_lines = ''
    grep_lines += '  {:>20s}\n'.format('_____________________________________________________Summary Results for GREP')
    grep_lines += line_hdr
    grep_lines += line_units
    grep_lines += line_data
    grep_lines += '\n'
    
    # Command line
    cmd_line = ''
    cmd_line += '  {:>20s}\n'.format('_________________________________________________________________Command line')
    cmd_line += '  {}\n'.format('Edit the command below by selecting replaced labels carefully.')
    cmd_line += '  {}\n'.format('Should some labels be set to zero (i.e. eliminated)?')
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
             
    if yaml_file:
        ogo.message('Saving Quality Assurance report to file:')
        ogo.message('      \"{}\"'.format(yaml_file))
        with open(yaml_file, 'w') as file:
            documents = yaml.dump(info_dict, file, sort_keys=False)
            
    if report_file:
        ogo.message('Saving report to file:')
        ogo.message('      \"{}\"'.format(report_file))
        txt_file = open(report_file, "w")
        txt_file.write(report)
        txt_file.write(cmd_line)
        txt_file.write(grep_lines)
        txt_file.close()

    ogo.message('Printing report')
    print('\n')
    print(report)
    if label_repair_list != []:
        print(cmd_line)
    
    ogo.message('Done.')

# +------------------------------------------------------------------------------+
def repair(input_image, output_image, relabel_parts, threshold_by_min_volume, \
           threshold_by_max_number_parts, min_volume, max_number_parts, \
           skip_labels, overwrite, func):

    an_action_was_defined = False # Changed to True if any action is defined.
        
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

    # Set up SITK filters
    filt = sitk.LabelShapeStatisticsImageFilter()
    filt.Execute(ct)
    
    conn = sitk.ConnectedComponentImageFilter()
    conn.SetFullyConnected(True)
    
    stats = sitk.LabelIntensityStatisticsImageFilter()
    
    # Start with a copy of the original CT. We'll erase labels we don't want.
    ct_base = sitk.Image(ct)
    #ct_base = sitk.Image(ct)<0 # with all zeros
    
    # ---------------------------------------------------    
    # Relabel parts
    if (relabel_parts):
        an_action_was_defined = True
        
        if np.remainder(len(relabel_parts),3) != 0:
            os.sys.exit('[ERROR] Incorrect definition for --relabel_parts. Tuples of 3 are expected.')

        n_relabel_parts = int(len(relabel_parts) / 3)
        relabel_parts = np.reshape(relabel_parts, (n_relabel_parts,3))
        
        ogo.message('Relabelling {:d} parts in image:'.format(n_relabel_parts))
        for swap in relabel_parts:
            label = int(swap[0])
            part = int(swap[1])
            new_label = int(swap[2])
            
            if label not in labels:
                ogo.message('[ERROR] The label {} is not in image.'.format(label))
                os.sys.exit()
                        
            ct_thres = ct==label
            ct_conn = conn.Execute(ct,ct_thres)
            ct_conn_sorted = sitk.RelabelComponent(ct_conn, sortByObjectSize=True) # could use minimumObjectSize
                        
            filt.Execute(ct_conn_sorted)
            size = filt.GetPhysicalSize(part)
            
            ct_part = (ct_conn_sorted==part)
            
            try:
                desc_label = lb.labels_dict[label]['LABEL']
            except KeyError:
                desc_label = 'unknown label'

            try:
                desc_new_label = lb.labels_dict[new_label]['LABEL']
            except KeyError:
                desc_new_label = 'unknown label'
            
            ogo.message('  label {:d} ({:s}), part {} ({:.1f} mm3) --> {} ({})'.format(\
                           label, desc_label, part, size, new_label, desc_new_label))
            
            bin_part = ct_part>0
            mask = 1 - bin_part
            ct_base = sitk.Mask(ct_base, mask)
            ct_base = ct_base + new_label*bin_part

    # ---------------------------------------------------    
    # Threshold by volume or maximum number of parts
    if threshold_by_min_volume or threshold_by_max_number_parts:
        an_action_was_defined = True
        ogo.message('')
        
        if threshold_by_min_volume:
            ogo.message('Thresholding by minimum volume {:.1f} mm3'.format(min_volume))
            if min_volume < 0.0:
                ogo.message('[ERROR] Minimum volume must be greater than zero.')
                os.sys.exit()
            
        if threshold_by_max_number_parts:
            ogo.message('Thresholding by maximum number of parts {:d}'.format(max_number_parts))
            if max_number_parts < 1:
                ogo.message('[ERROR] Maximum number of parts must be greater than zero.')
                os.sys.exit()
                
        ogo.message('  {:>20s} {:>6s} {:>10s} {:>10s}'.format('LABEL','PART','VOL','STATUS'))
        ogo.message('  {:>20s} {:>6s} {:>10s} {:>10s}'.format('#','#','mm3',' '))

        for label in labels:
            try:
                desc = lb.labels_dict[label]['LABEL']
            except KeyError:
                desc = 'unknown label'
            
            if label in skip_labels:
                ogo.message('  {:>15s}{:>5s} {:>6s} {:>10s} {:>10s}'.format('('+desc+')',str(label),'---','---','skip'))
                
            else:
                  
                ct_thres = ct_base==label
                ct_conn = conn.Execute(ct_base,ct_thres)
                ct_conn_sorted = sitk.RelabelComponent(ct_conn, sortByObjectSize=True) # could use minimumObjectSize
                
                stats.Execute(ct_conn_sorted,ct)
                n_parts = stats.GetNumberOfLabels()
                
                for part in stats.GetLabels():
                    keep=True
                    size = stats.GetPhysicalSize(part)
                    
                    status = ''
                    if threshold_by_min_volume:
                        if size<min_volume:
                            status = 'remove'
                            keep=False
                    if threshold_by_max_number_parts:
                        if part>max_number_parts:
                            status = 'remove'
                            keep=False
                            
                    if not keep:
                        ct_part = (ct_conn_sorted==part)
                        bin_part = ct_part>0
                        mask = 1 - bin_part
                        ct_base = sitk.Mask(ct_base, mask)
                        ct_base = ct_base + 0*bin_part # erase label
                    
                    if part==1:
                        ogo.message('  {:>15s}{:>5s} {:6d} {:10.1f} {:>10s}'.format('('+desc+')',str(label),part,size,status))
                    else:
                        ogo.message('  {:>15s}{:>5s} {:6d} {:10.1f} {:>10s}'.format(' ',' ',part,size,status))
                                            
    # Exit if no actions were defined
    if not an_action_was_defined:
        ogo.message('[ERROR]: No actions were defined. At least one option must be defined:')
        ogo.message('         --relabel_parts')
        ogo.message('         --threshold_by_min_volume')
        ogo.message('         --threshold_by_max_number_parts')
        os.sys.exit()
                
    final_labels = get_labels(ct_base)
    ogo.message('Final image contains the following labels:')
    ogo.message('  [' + ', '.join('{:d}'.format(i) for i in final_labels) + ']')
        
    ogo.message('')
    ogo.message('Writing merged output image to file:')
    ogo.message('  {}'.format(output_image))
    sitk.WriteImage(ct_base, output_image)            
    
    # Command line
    cmd_line = ''
    cmd_line += '  {:>20s}\n'.format('__________________________________________________________________Command line')
    cmd_line += '\n'
    cmd_line += '  {}\n'.format('ogoValidate validate \\')
    cmd_line += '    {} \n'.format(output_image)    
    
    print(cmd_line)
    
    ogo.message('Done.')
    
    
def main():
    # Setup description
    description = '''
This tools validates the segmented output of machine learning tools to ensure 
high quality and appropriate labelling of skeletal parts. As opposed to methods 
such as DICE scores and Hausdorff distance, it does quality assurance tests by
ensuring certain anatomical constraints are met. To learn more:

ogoValidate validate -h

The segmentation can be repaired using this tool if there is an existing label 
that is incorrectly assigned. The user can define the new label by defining 
which label and part is to be assigned a new label. To learn more:

ogoValidate repair -h

If some of the segmented output is not labelled (i.e. it didn't apply a label
to part of the bone, either manual labelling will be required, or the user is
directed to use Graph Cuts. See ogoGraphCuts to learn more.

Definitions:

Label – a value up to 255 that identifies the type of tissue.
Part – for each label it is represented by unconnected parts.

Validation is based on the anatomical constraints and checks that are tested:
  
- The input image contains all expected labels (the list may be user defined)
- A label has only one part (exceptions for Sacrum and L5, which may contain 2
  due to a Pars fracture)
- The left and right femur are the same volume (within a given tolerance)
- The left and right pelvis are the same volume (within a given tolerance)
- Each bone is within the expected min/max volume range that is pre-defined
  (A future implementation may refine the size limits relative to a reference
  bone, such as the sacrum. Large skeletons have large bones, small skeletons...)
- The position of each bone relative to other bones is appropriate as 
  assessed by a Procrustes analysis

Repairs can be performed the following ways:

- Label and part is defined and assigned a new label (relabel_parts)
– Remove parts beyond the maximum number of parts (threshold_by_max_number_parts)
- Remove any \'parts\' of a given volume (threshold_by_min_volume)
  
– A list of labels excluded from repairs using volume or parts can be defined
  (skip_labels)

Notes on how this works:

The validation cycles through the list of expected labels. For each label (e.g. 
Left Femur) it determines how many connected components there are (parts). If
a bone is correctly labelled there should only be one part. However, it is 
possible that only one part is defined, but that the label does not cover the
entire bone (e.g. a portion of Left Femur is labelled something else, like Right
Femur). Without visualizing the results of segmentation, this is difficult to 
detect. However, typically the volume of that bone label will be smaller than
normal, which is then used to identify that something is wrong. Furthermore, if
you are looking at a bone that has a symmetric pair (femur, pelvis), then a 
warning is applied if their volumes are different more than a given threshold.
Interestingly, the sacrum and L5 sometimes seem to be separated into two parts,
which is unusual but not impossible.

The YAML file output during validation captures the results of all tests performed
so that subsequent analysis of large datasets can (a) summarize the quality checks,
(b) identify which parts of a segmentation are valid, or (c) combine multiple 
model outputs from the same skeleton to capture the bones that have met the
quality checks.

'''

    epilog = '''
Example calls: 
ogoValidate validate image.nii.gz

# Relabels label 3, part 2 to label 5
ogoValidate repair image.nii.gz image_repair.nii.gz --relabel_parts 3 2 5

# Relabels label 3, part 2 to label 5 and label 3, part 3 to label 5
ogoValidate repair image.nii.gz image_repaired.nii.gz \\
    --relabel_parts \\
    3 2 5 \\
    3 3 5
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
    parser_validate.add_argument('--yaml_file', metavar='FILE', help='Write quality assurance (QA) outputs (*.yaml)')
    parser_validate.add_argument('--expected_labels', type=int, nargs='*', default=[1,2,3,4,5,6,7,8,9,10,11], metavar='LABEL', help='List of labels expected in image (default: %(default)s)')
    parser_validate.add_argument('--overwrite', action='store_true', help='Overwrite validation report without asking')
    parser_validate.set_defaults(func=validate)

    # Repair
    parser_repair = subparsers.add_parser('repair')
    parser_repair.add_argument('input_image', help='Input image file (*.nii, *.nii.gz)')
    parser_repair.add_argument('output_image', help='Output image file (*.nii, *.nii.gz)')
    parser_repair.add_argument('--relabel_parts', type=int, nargs='*', default=[], metavar='LABEL PART NEWLABEL', help='List of labels, parts, and new labels; space separated. Example shows relabelling L4 part 2 to  L3: (e.g. 7 2 8)')
    parser_repair.add_argument('--threshold_by_min_volume', action='store_true', help='Boolean setting to remove all parts smaller than the minimum volume (mm3)')
    parser_repair.add_argument('--threshold_by_max_number_parts', action='store_true', help='Boolean setting to remove all parts greater than maximum number of parts')
    parser_repair.add_argument('--min_volume', type=float, default=20000.0, metavar='VOL', help='Set minimum volume threshold (default: %(default)s mm3)')
    parser_repair.add_argument('--max_number_parts', type=int, default=1, metavar='N_PARTS', help='Set maximum number of parts (default: %(default)s)')
    parser_repair.add_argument('--skip_labels', type=int, nargs='*', default=[], metavar='LABEL', help='List of labels excluded for repair (default: %(default)s)')
    parser_repair.add_argument('--overwrite', action='store_true', help='Overwrite output image without asking')
    parser_repair.set_defaults(func=repair)

    # Parse and display
    args = parser.parse_args()
    print(echo_arguments('Validate', vars(args)))

    # Run program
    args.func(**vars(args))


if __name__ == '__main__':
    main()
