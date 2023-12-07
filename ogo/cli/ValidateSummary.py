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

def ValidateSummary(input_files,seek_label_pass,seek_label_fail,seek_procrustes_pass,seek_procrustes_fail):

    summary_dict = {}
    
    test_result_dict = {
        'procrustes_num':0,
        'procrustes_percent':0.0,
        'all_found_num' :0,
        'all_found_percent' :0.0,
        'intact_femur_num':0,
        'intact_femur_percent':0.0,
        'intact_spine_num':0,
        'intact_spine_percent':0.0,
        'intact_any_spine_num':0,
        'intact_any_spine_percent':0.0,
        'intact_all_labels_num':0,
        'intact_all_labels_percent':0.0,
        'final_num':0,
        'final_percent':0.0,
        'total':0
    }
        
    # Check if output exists and should overwrite
    #ogo.pass_check_if_output_exists(output_file,overwrite)

    # Set up to read image and mask inputs
    ogo.pass_check_if_file_exists(input_files[0])

    # Check endings
    ogo.pass_check_file_ending(input_files[0],['.yaml'])
    #ogo.pass_check_file_ending(output_file,['.csv'])
    
    lst_label_pass = []
    lst_label_fail = []
    lst_procrustes_pass = []
    lst_procrustes_fail = []
    message_trigger = np.ceil(np.linspace(0,len(input_files),20))
    
    # Read input
    idx=0
    for idx,f in enumerate(input_files):
        if idx in message_trigger:
            ogo.message('Reading {:4d} of {:d}: {}'.format(idx+1,len(input_files),f))
        with open(f, 'r') as file:
            result_dict = yaml.safe_load(file)
        
        fileinfo_dict = result_dict.get('fileinfo')
        
        check_dict = result_dict.get('check')
        if not check_dict:
            ogo.message('[ERROR] No \'check\' dictionary key in {}.'.format(f))
            os.sys.exit()
        
        test_result_dict['total'] += 1
        
        if check_dict['procrustes']['status']:
            test_result_dict['procrustes_num'] += 1
            
        if seek_procrustes_pass and check_dict['procrustes']['status']:
            if fileinfo_dict['basename'] not in lst_procrustes_pass:
                lst_procrustes_pass.append([fileinfo_dict['basename'],check_dict['procrustes']['disparity']])
        if seek_procrustes_fail and not check_dict['procrustes']['status']:
            if fileinfo_dict['basename'] not in lst_procrustes_fail:
                lst_procrustes_fail.append([fileinfo_dict['basename'],check_dict['procrustes']['disparity']])
        
        if check_dict['all_found']['status']:
            test_result_dict['all_found_num'] += 1

        if check_dict['intact_right_femur']['status'] or check_dict['intact_left_femur']['status']:
            test_result_dict['intact_femur_num'] += 1
        
        if check_dict['labels'][1]['status']:
            test_result_dict['intact_any_spine_num'] += 1
            
        if check_dict['intact_spine']['status']:
            test_result_dict['intact_spine_num'] += 1
            
        if check_dict['intact_all_labels']['status']:
            test_result_dict['intact_all_labels_num'] += 1

        if check_dict['final']['status']:
            test_result_dict['final_num'] += 1
                   
        labels_dict = check_dict.get('labels')
        for label_key, label_value in labels_dict.items():
            if summary_dict.get(label_key) is None:
                summary_dict[label_key] = {}
                summary_dict[label_key]['label_name'] = label_value['name']
                summary_dict[label_key]['label_num'] = 0
                summary_dict[label_key]['pass_num'] = 0
                summary_dict[label_key]['pass_percent'] = 0
                summary_dict[label_key]['volume_sum'] = []
                summary_dict[label_key]['volume_mean'] = 0.0
                summary_dict[label_key]['volume_stdev'] = 0.0
                summary_dict[label_key]['eFRAG_num'] = 0
                summary_dict[label_key]['eMINVOL_num'] = 0
                summary_dict[label_key]['eMAXVOL_num'] = 0
                summary_dict[label_key]['wNOFEA_num'] = 0
                summary_dict[label_key]['wSPEC_num'] = 0
                summary_dict[label_key]['wSYMM_num'] = 0
                summary_dict[label_key]['eFRAG_percent'] = 0.0
                summary_dict[label_key]['eMINVOL_percent'] = 0.0
                summary_dict[label_key]['eMAXVOL_percent'] = 0.0
                summary_dict[label_key]['wNOFEA_percent'] = 0.0
                summary_dict[label_key]['wSPEC_percent'] = 0.0
                summary_dict[label_key]['wSYMM_percent'] = 0.0
            
            if label_value['n_parts'] > 0: # checks only bones that existed (e.g., 0 parts means there was no bone)
                summary_dict[label_key]['label_num'] += 1
                summary_dict[label_key]['volume_sum'].append(label_value['total_volume'])
                
                if label_value['status']:
                    summary_dict[label_key]['pass_num'] += 1
                    if label_key in seek_label_pass:
                        if fileinfo_dict['basename'] not in lst_label_pass:
                            lst_label_pass.append(fileinfo_dict['basename'])
                else:
                    if label_key in seek_label_fail:
                        if fileinfo_dict['basename'] not in lst_label_fail:
                            lst_label_fail.append(fileinfo_dict['basename'])
                        
                
            if 'eFRAG' in label_value['errors']:
                summary_dict[label_key]['eFRAG_num'] += 1
            if 'eMINVOL' in label_value['errors']:
                summary_dict[label_key]['eMINVOL_num'] += 1
            if 'eMAXVOL' in label_value['errors']:
                summary_dict[label_key]['eMAXVOL_num'] += 1
            
            if 'wNOFEA' in label_value['warnings']:
                summary_dict[label_key]['wNOFEA_num'] += 1
            if 'wSPEC' in label_value['warnings']:
                summary_dict[label_key]['wSPEC_num'] += 1
            if 'wSYMM' in label_value['warnings']:
                summary_dict[label_key]['wSYMM_num'] += 1

        idx += 1
    
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
        
        if summary_dict[label_key]['volume_sum']:
            summary_dict[label_key]['volume_mean'] = np.mean(summary_dict[label_key]['volume_sum'])
            summary_dict[label_key]['volume_stdev'] = np.std(summary_dict[label_key]['volume_sum'])
        else:
            summary_dict[label_key]['volume_mean'] = 0.0
            summary_dict[label_key]['volume_stdev'] = 0.0
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
    
Reads in a list of YAML files that were created by ogoValidate on a series 
of inference results. 

The main outcome is a table printed to screen that summarizes the success
of inference.

Additional possibilities will identify specific inference results based on
a criteria:

  seek_label_pass: Lists inference results for files where at least one of the
             listed labels has passed.
  seek_label_fail: Lists inference results for files where at least one of the 
             listed labels has failed.
  seek_procrustes_pass: Lists inference results for files where Procrustes
             has passed (includes estimate of disparity)
  seek_procrustes_fail: Lists inference results for files where Procrustes
             has failed (includes estimate of disparity)
             
'''

    epilog = '''
Example call: 
     
ogoValidateSummary RETRO_000*.yaml
ogoValidateSummary RETRO_000*.yaml --seek_label_pass 7 8 9 10
ogoValidateSummary RETRO_000*.yaml --seek_procrustes_pass

'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoValidateSummary",
        description=description,
        epilog=epilog
    )

    parser.add_argument('input_files', nargs='+', metavar='YAML', help='Input YAML files (*.yaml)')

    #parser.add_argument('--output_file', default=None, metavar='CSV',help='Summary output table (*.csv, default: %(default)s)')
    parser.add_argument('--seek_label_pass', type=int, nargs='*', default=[], metavar='LABEL', help='List images containing pass labels (default: %(default)s)')
    parser.add_argument('--seek_label_fail', type=int, nargs='*', default=[], metavar='LABEL', help='List images containing fail labels (default: %(default)s)')
    parser.add_argument('--seek_procrustes_pass', action='store_true',help='List images containing pass labels (default: %(default)s)')
    parser.add_argument('--seek_procrustes_fail', action='store_true',help='List images containing fail labels (default: %(default)s)')
    #parser.add_argument('--overwrite', action='store_true',help='Overwrite output file without asking')

    # Parse and display
    args = parser.parse_args()

    # Run program
    ValidateSummary(**vars(args))


if __name__ == '__main__':
    main()
