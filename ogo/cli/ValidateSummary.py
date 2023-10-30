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

# ============================================================================================================================
#                Label         Pass       Volume      eMINVOL      eMAXVOL        eFRAG        eSYMM       wNOFEA        wSPEC
#             name (N)        % (N)   mean (std)        % (N)        % (N)        % (N)        % (N)        % (N)        % (N)
# -------------------- ------------ ------------ ------------ ------------ ------------ ------------ ------------ ------------
#   Femur Right (3993)  87.3( 3486) 142.9( 38.8)   1.4(   56)   0.7(   27)   4.5(  179)   8.1(  324)  60.9( 2432)             
#    Femur Left (3992)  91.2( 3640) 144.2( 39.1)   1.3(   51)   0.6(   24)   3.8(  153)   5.0(  199)  59.9( 2392)             
#  Pelvis Right (3994)  82.5( 3294) 322.0( 69.1)   0.9(   36)   0.2(    7)  16.5(  660)   0.2(    9)                          
#   Pelvis Left (3994)  82.9( 3312) 322.9( 69.6)   1.0(   38)   0.2(    9)  16.1(  642)   0.3(   11)                          
#        Sacrum (3994)  73.3( 2926) 215.6( 41.6)   1.8(   71)   0.2(    8)  25.4( 1014)                             2.9(  114)
#            L5 (3994)  64.6( 2579)  61.8( 14.3)   4.1(  162)   0.8(   33)  33.1( 1323)                             2.6(  104)
#            L4 (3993)  63.7( 2543)  61.2( 16.2)   5.2(  206)   1.2(   47)  34.7( 1386)                                       
#            L3 (3992)  66.5( 2655)  59.2( 16.2)   5.3(  213)   1.4(   56)  32.5( 1297)                                       
#            L2 (3991)  70.5( 2814)  54.5( 13.7)   5.3(  210)   1.1(   42)  27.6( 1101)                                       
#            L1 (3991)  62.9( 2509)  51.8( 12.4)   1.6(   65)   0.4(   15)  36.3( 1449)                                       
#           T12 (1455)  50.7(  738)  25.0( 27.1)   0.0(    0)   0.0(    0)  49.3(  717)                                       
# -------------------- ------------
#    Procrustes (3994)  56.6( 2259)
#  Labels Found (3994)  36.4( 1452)
#         FINAL (3994)   0.0(    0)
# ============================================================================================================================
# Table caption:
#   Pass rate is shown per bone for all the images.
#   Pass rate of all bone segmentations is summarized as FINAL.
#   Labels found means the image contained all expected labels.
#   Procrustes analysis shows skeletons with unusual anatomy (e.g., L6 present).
#   Errors for each bone include:
#     eMINVOL –- exceeds min volume limit
#     eMAXVOL –- exceeds max volume limit
#     eFRAG   –- fragmented
#   Warnings for each bone include:
#     wSYMM   –- not symmetric
#     wNOFEA  –- femur is cut off too much
#     wSPEC   –- Pars defect for L5 or sacrum tip is disassociated


# +------------------------------------------------------------------------------+
# Writer
def write(entry, deliminator, ofile):
    if ofile is None:
        output = os.sys.stdout
    else:
        if os.path.isfile(ofile):
            output = open(ofile, 'a')
        else:
            output = open(ofile, 'w')
            output.write(deliminator.join([str(x) for x in write.header]))
            output.write(os.linesep)

    output.write(deliminator.join([str(x) for x in entry]))
    output.write(os.linesep)

    if output is not os.sys.stdout:
        output.close()

# +------------------------------------------------------------------------------+
#def ValidateSummary(input_file, output_file, overwrite):
def ValidateSummary(input_files,seek_label_pass,seek_label_fail,seek_procrustes_pass,seek_procrustes_fail):

    verbose = False
    
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
     
ValidateSummary RETRO_000*.yaml
ValidateSummary RETRO_000*.yaml --seek_label_pass 7 8 9 10
ValidateSummary RETRO_000*.yaml --seek_procrustes_pass

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
