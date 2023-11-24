# /------------------------------------------------------------------------------+
# | 28-OCT-2023                                                                  |
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
import yaml
import json
import numpy as np
import SimpleITK as sitk
from datetime import date
from ogo.util.echo_arguments import echo_arguments
import ogo.util.Helper as ogo
import queue
import re
import collections
import pandas as pd
import ogo.dat.OgoMasterLabels as lb

def ValidateSelectBest(input_files,output_file,overwrite):

    best_dict = {}

    pass_final = []
    pass_intact_all_labels = []
    pass_intact_spine = []
    pass_intact_right_femur = []
    pass_intact_left_femur = []
    pass_procrustes = []
    
    # Check if output exists and should overwrite
    ogo.pass_check_if_output_exists(output_file,overwrite)

    # Set up to read image and mask inputs
    ogo.pass_check_if_file_exists(input_files[0])

    # Check endings
    ogo.pass_check_file_ending(input_files[0],['.yaml'])
    ogo.pass_check_file_ending(output_file,['.nii.gz'])
    
    # Read input
    idx=0
    for idx,f in enumerate(input_files):
        ogo.message('Reading {}'.format(f))
        with open(f, 'r') as file:
            result_dict = yaml.safe_load(file)
        
        check_dict = result_dict.get('check')
        if not check_dict:
            ogo.message('[ERROR] No \'check\' dictionary key in {}.'.format(f))
            os.sys.exit()
        
        fileinfo_dict = result_dict.get('fileinfo')
        source_basename = fileinfo_dict['basename']

        input_parameters_dict = result_dict.get('input_parameters')
        source_input_image = input_parameters_dict['input_image']
        
        if check_dict['final']['status']:
            pass_final.append(source_input_image)

        if check_dict['intact_all_labels']['status']:
            pass_intact_all_labels.append(source_input_image)
        
        if check_dict['intact_spine']['status']:
            pass_intact_spine.append(source_input_image)
        
        if check_dict['intact_left_femur']['status']:
            pass_intact_left_femur.append(source_input_image)

        if check_dict['intact_right_femur']['status']:
            pass_intact_right_femur.append(source_input_image)

        if check_dict['procrustes']['status']:
            pass_procrustes.append(source_input_image)
        
        labels_dict = check_dict.get('labels')
        for label_key, label_value in labels_dict.items():
            if best_dict.get(label_key) is None:
                best_dict[label_key] = {}
                best_dict[label_key]['label_name'] = label_value['name']
                best_dict[label_key]['label_label'] = label_value['label']
                best_dict[label_key]['label_status'] = False
                best_dict[label_key]['label_source'] = ''
                best_dict[label_key]['label_input_image'] = ''
            
            if best_dict[label_key]['label_status'] is False: # ensures we use first good label encountered
                if label_value['status'] and label_value['n_parts']>0:
                    best_dict[label_key]['label_status'] = True
                    best_dict[label_key]['label_source'] = source_basename
                    best_dict[label_key]['label_input_image'] = source_input_image
                
    ogo.message('Passed Final:')
    for entry in pass_final:
        ogo.message('  {}'.format(entry))
    ogo.message('Passed Intact All Labels:')
    for entry in pass_intact_all_labels:
        ogo.message('  {}'.format(entry))
    ogo.message('Passed Intact Spine:')
    for entry in pass_intact_spine:
        ogo.message('  {}'.format(entry))
    ogo.message('Passed Intact Right Femur:')
    for entry in pass_intact_right_femur:
        ogo.message('  {}'.format(entry))
    ogo.message('Passed Intact Left Femur:')
    for entry in pass_intact_left_femur:
        ogo.message('  {}'.format(entry))
    ogo.message('Passed Procrustes:')
    for entry in pass_procrustes:
        ogo.message('  {}'.format(entry))
        
    ogo.message('Frankenstein Possibility:')
    for k,v in best_dict.items():
        print('{:>20s} {:5d} {} {:s}'.format(v['label_name'],v['label_label'],v['label_status'],v['label_source']))
    
    exit()

    # Assemble output image
    if output_file:
        ogo.message('Assembling output image')
    
        ct_base = None
        
        for k,v in best_dict.items():
            ct_name = v['label_input_image']
            ct = sitk.ReadImage(ct_name, sitk.sitkUInt8)
            
            if ct_base is None:
                ct_base = ct<0
                
            if v['label_status']:
                
                label = v['label_label']
                ogo.message('  Gathering label {} from {}'.format(label,ct_name))
                bin_part = ct==label
                mask = 1 - bin_part
                ct_base = sitk.Mask(ct_base, mask)
                ct_base = ct_base + label*bin_part
                
    exit()
    # Calculate percentages, means and stdev
    for label_key, label_value in labels_dict.items():
        if summary_dict[label_key]['eFRAG_num']>0:
            summary_dict[label_key]['eFRAG_percent'] = 100.0 * summary_dict[label_key]['eFRAG_num'] / summary_dict[label_key]['label_num']
        if summary_dict[label_key]['eMINVOL_num']>0:
            summary_dict[label_key]['eMINVOL_percent'] = 100.0 * summary_dict[label_key]['eMINVOL_num'] / summary_dict[label_key]['label_num']
        if summary_dict[label_key]['eMAXVOL_num']>0:
            summary_dict[label_key]['eMAXVOL_percent'] = 100.0 * summary_dict[label_key]['eMAXVOL_num'] / summary_dict[label_key]['label_num']
        if summary_dict[label_key]['wNOFEA_num']>0:
            summary_dict[label_key]['wNOFEA_percent'] = 100.0 * summary_dict[label_key]['wNOFEA_num'] / summary_dict[label_key]['label_num']
        if summary_dict[label_key]['wSPEC_num']>0:
            summary_dict[label_key]['wSPEC_percent'] = 100.0 * summary_dict[label_key]['wSPEC_num'] / summary_dict[label_key]['label_num']
        if summary_dict[label_key]['wSYMM_num']>0:
            summary_dict[label_key]['wSYMM_percent'] = 100.0 * summary_dict[label_key]['wSYMM_num'] / summary_dict[label_key]['label_num']
        
        if summary_dict[label_key]['pass_num']>0:
            summary_dict[label_key]['pass_percent'] = 100.0 * summary_dict[label_key]['pass_num'] / summary_dict[label_key]['label_num']

        summary_dict[label_key]['volume_mean'] = np.mean(summary_dict[label_key]['volume_sum'])
        summary_dict[label_key]['volume_stdev'] = np.std(summary_dict[label_key]['volume_sum'])
        summary_dict[label_key]['volume_sum'].clear()
        
        test_result_dict['procrustes_percent'] = 100.0 * test_result_dict['procrustes_num'] / test_result_dict['total']        
        test_result_dict['all_found_percent'] = 100.0 * test_result_dict['all_found_num'] / test_result_dict['total']        
        test_result_dict['intact_femur_percent'] = 100.0 * test_result_dict['intact_femur_num'] / test_result_dict['total']
        test_result_dict['intact_spine_percent'] = 100.0 * test_result_dict['intact_spine_num'] / test_result_dict['total']        
        test_result_dict['intact_any_spine_percent'] = 100.0 * test_result_dict['intact_any_spine_num'] / test_result_dict['total']        
        test_result_dict['intact_all_labels_percent'] = 100.0 * test_result_dict['intact_all_labels_num'] / test_result_dict['total']        
        test_result_dict['final_percent'] = 100.0 * test_result_dict['final_num'] / test_result_dict['total']        
               
    ogo.message('Read {} files.'.format(idx))
    
    # Print successful and unsuccessful labels
    if lst_label_pass:
        ogo.message('Files containing successful labels '+' '.join('{:d}'.format(i) for i in seek_label_pass))
        ogo.message('  {:d} files'.format(len(lst_label_pass)))
        for fname in lst_label_pass:
            ogo.message ('  {}'.format(fname))
    if lst_label_fail:
        ogo.message('Files containing unsuccessful labels '+' '.join('{:d}'.format(i) for i in seek_label_fail))
        ogo.message('  {:d} files'.format(len(lst_label_fail)))
        for fname in lst_label_fail:
            ogo.message ('  {}'.format(fname))
    if seek_procrustes_pass:
        ogo.message('Files passing Procrustes result')
        ogo.message('  {:d} files'.format(len(lst_procrustes_pass)))
        for fname in lst_procrustes_pass:
            ogo.message ('  {} {:7.5f}'.format(fname[0],fname[1]))
    if seek_procrustes_fail:
        ogo.message('Files failing Procrustes result')
        ogo.message('  {:d} files'.format(len(lst_procrustes_fail)))
        for fname in lst_procrustes_fail:
            ogo.message ('  {} {:7.5f}'.format(fname[0],fname[1]))
    
    # Print to screen
    print('{:=>20s}={:=>12s}={:=>12s}={:=>12s}={:=>12s}={:=>12s}={:=>12s}={:=>12s}={:=>12s}'.format(\
          '','','','','','','','',''))
    print('{:>20s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s}'.format(\
          'Label','Pass','Volume','eMINVOL','eMAXVOL','eFRAG','wSYMM','wNOFEA','wSPEC'))
    print('{:>20s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s}'.format(\
          '(N)','%   (N)','mean (std)','%   (N)','%   (N)','%   (N)','%   (N)','%   (N)','%   (N)'))
    print('{:->20s} {:->12s} {:->12s} {:->12s} {:->12s} {:->12s} {:->12s} {:->12s} {:->12s}'.format(\
          '-','-','-','-','-','-','-','-','-'))
          
    for label_key, label_value in labels_dict.items():
        if label_key in list(lb.labels_dict.keys()):
            
            str_label = '{:s} {:4d}'.format(summary_dict[label_key]['label_name'].replace(" ", "_"),summary_dict[label_key]['label_num'])
            str_status = '{:5.1f} {:5d}'.format(summary_dict[label_key]['pass_percent'],summary_dict[label_key]['pass_num'])
            str_volume = '{:5.1f} {:5.1f}'.format(summary_dict[label_key]['volume_mean']/1000.0,summary_dict[label_key]['volume_stdev']/1000.0)
            str_eMINVOL = '{:5.1f} {:5d}'.format(summary_dict[label_key]['eMINVOL_percent'],summary_dict[label_key]['eMINVOL_num'])
            str_eMAXVOL = '{:5.1f} {:5d}'.format(summary_dict[label_key]['eMAXVOL_percent'],summary_dict[label_key]['eMAXVOL_num'])
            str_eFRAG = '{:5.1f} {:5d}'.format(summary_dict[label_key]['eFRAG_percent'],summary_dict[label_key]['eFRAG_num'])

            if label_key>=1 and label_key<=4:
                str_wSYMM = '{:5.1f} {:5d}'.format(summary_dict[label_key]['wSYMM_percent'],summary_dict[label_key]['wSYMM_num'])
            else:
                str_wSYMM = '{:>5s} {:>5s}'.format('-','-')
            
            if label_key>=1 and label_key<=2:                                        
                str_wNOFEA = '{:5.1f} {:5d}'.format(summary_dict[label_key]['wNOFEA_percent'],summary_dict[label_key]['wNOFEA_num'])
            else:
                str_wNOFEA = '{:>5s} {:>5s}'.format('-','-')
            
            if label_key>=5 and label_key<=6:                                                    
                str_wSPEC = '{:5.1f} {:5d}'.format(summary_dict[label_key]['wSPEC_percent'],summary_dict[label_key]['wSPEC_num'])
            else:
                str_wSPEC = '{:>5s} {:>5s}'.format('-','-')
            
            print('{:>20s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s}'.format(\
                  str_label,\
                  str_status,\
                  str_volume,\
                  str_eMINVOL,\
                  str_eMAXVOL,\
                  str_eFRAG,\
                  str_wSYMM,\
                  str_wNOFEA,\
                  str_wSPEC))
            

    print('{:->20s} {:->12s}'.format(\
          '-','-'))
          
    str_label = '{:s} {:>4d}'.format('Procrustes',test_result_dict['total'])
    str_status = '{:5.1f} {:5d}'.format(test_result_dict['procrustes_percent'],test_result_dict['procrustes_num'])
    print('{:>20s}  {:>10s}'.format(str_label,str_status))

    str_label = '{:s} {:>4d}'.format('All_Found',test_result_dict['total'])
    str_status = '{:5.1f} {:5d}'.format(test_result_dict['all_found_percent'],test_result_dict['all_found_num'])
    print('{:>20s}  {:>10s}'.format(str_label,str_status))
    
    str_label = '{:s} {:>4d}'.format('Intact_Femur',test_result_dict['total'])
    str_status = '{:5.1f} {:5d}'.format(test_result_dict['intact_femur_percent'],test_result_dict['intact_femur_num'])
    print('{:>20s}  {:>10s}'.format(str_label,str_status))

    str_label = '{:s} {:>4d}'.format('Spine_all',test_result_dict['total'])
    str_status = '{:5.1f} {:5d}'.format(test_result_dict['intact_spine_percent'],test_result_dict['intact_spine_num'])
    print('{:>20s}  {:>10s}'.format(str_label,str_status))
    
    str_label = '{:s} {:>4d}'.format('Spine_any',test_result_dict['total'])
    str_status = '{:5.1f} {:5d}'.format(test_result_dict['intact_any_spine_percent'],test_result_dict['intact_any_spine_num'])
    print('{:>20s}  {:>10s}'.format(str_label,str_status))
    
    str_label = '{:s} {:>4d}'.format('All',test_result_dict['total'])
    str_status = '{:5.1f} {:5d}'.format(test_result_dict['intact_all_labels_percent'],test_result_dict['intact_all_labels_num'])
    print('{:>20s}  {:>10s}'.format(str_label,str_status))
    
    str_label = '{:s} {:>4d}'.format('FINAL',test_result_dict['total'])
    str_status = '{:5.1f} {:5d}'.format(test_result_dict['final_percent'],test_result_dict['final_num'])
    print('{:>20s}  {:>10s}'.format(str_label,str_status))
                
    print('{:=>20s}={:=>12s}={:=>12s}={:=>12s}={:=>12s}={:=>12s}={:=>12s}={:=>12s}={:=>12s}'.format(\
          '','','','','','','','',''))

    print('Table caption:')
    print('  Procrustes analysis shows skeletons with unusual anatomic alignment (with or without L6).')
    print('  All found means all expected labels were in the image.')
    print('  Intact tests are positive if the bone(s) do not contain eFRAG error:')
    print('    Femur refers to left femur OR right femur.')
    print('    All spine refers to L1, L2, L3, AND L4.')
    print('    Any spine refers to L1, L2, L3, OR L4.')
    print('    All refers to all expected bone(s).')
    print('  Errors for each bone include:')
    print('    eMINVOL –- exceeds min volume limit')
    print('    eMAXVOL –- exceeds max volume limit')
    print('    eFRAG   –- fragmented')
    print('  Warnings for each bone include:')
    print('    wSYMM   –- not symmetric')
    print('    wNOFEA  –- femur is cut off too much')
    print('    wSPEC   –- Pars defect for L5 or sacrum tip is disassociated')
     
def main():
    # Setup description
    description = '''
    
Reads in a list of YAML files representing inference of the same inference but 
conducted with different ML models. It selects the combination of inferences 
so that the Frankenstein segmentation is the best possible combination.

Typically you would have more than one model to perform inference. For example, 
you may have three YAML files for the same skeleton. This will output the best 
inference result possible. 
             
'''

    epilog = '''
Example call: 
     
ogoValidateSelectBest RETRO_00196_model1.yaml RETRO_00196_model2.yaml --output RETRO_00196_best.yaml

'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoValidateSelectBest",
        description=description,
        epilog=epilog
    )

    parser.add_argument('input_files', nargs='+', metavar='YAML', help='Input YAML files (*.yaml)')

    parser.add_argument('--output_file', default=None, metavar='YAML',help='Best result (*.yaml, default: %(default)s)')
    parser.add_argument('--overwrite', action='store_true',help='Overwrite output file without asking')

    # Parse and display
    args = parser.parse_args()

    # Run program
    ValidateSelectBest(**vars(args))


if __name__ == '__main__':
    main()
