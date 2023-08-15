# /------------------------------------------------------------------------------+
# | 14-JUL-2023                                                                  |
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
import SimpleITK as sitk
from datetime import date
from ogo.util.echo_arguments import echo_arguments
import ogo.util.Helper as ogo
import queue
import re
import collections
import pandas as pd
import ogo.dat.OgoMasterLabels as lb

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
def SummarizeValidationData(input_file, output_file, overwrite):

    verbose = False
    
    # Define possible errors and warnings
    error_list = ['eFRAG','eSYMM','eMINVOL','eMAXVOL']
    warning_list = ['wNOFEA', 'wSPEC']
    
    # Check if output exists and should overwrite
    ogo.pass_check_if_output_exists(output_file,overwrite)

    # Set up to read image and mask inputs
    ogo.pass_check_if_file_exists(input_file)

    # Check endings
    ogo.pass_check_file_ending(input_file,['.csv'])
    ogo.pass_check_file_ending(output_file,['.csv'])
    
    # Read input
    df = pd.read_csv(input_file, skiprows=[1])

    column_list = df.columns.values.tolist()
    #print(column_list)
    
    # Get all the filenames in case helpful
    df_filenames = df.get('name')

    results_dict = {}
    
    # Loop through all labels (maximum 255)
    for i in range(255):
        if i==0:
            suffix=''
        else:
            suffix='.{:d}'.format(i)

        df_desc = df.get('desc'+suffix, default='unknown')
                
        if df_desc is not 'unknown':
            if verbose:
                print('Reading suffix {}:{:s}'.format(i,suffix))
            this_desc = df_desc.values[0]
            
            if df_desc.nunique() > 1:
                ogo.message('[ERROR] Column does not contain unique values: N={}\n{}'.format(df_desc.nunique(),df_desc))
                os.sys.exit(0)
                
            df_label = df.get('label'+suffix)
            this_label = df_label.values[0]
            this_count = df.loc[:,'total_vol'+suffix].notnull().sum()
            
            df_status = df.get('status'+suffix)
            if df_status.isnull().all():
                this_status_pass = 0
                this_percent_status_pass = 0
            else:
                this_status_pass = df_status.str.contains(r'PASS').sum()
                this_percent_status_pass = this_status_pass / this_count * 100.0
            
            # Parse the warnings
            df_warnings = df.get('warnings'+suffix)
            if df_warnings.isnull().all():
                this_wNOFEA_count = 0
                this_wSPEC_count = 0
            else:
                if this_label == 1 or this_label == 2: # only left or right femur
                    this_wNOFEA_count = df_warnings.str.contains(r'wNOFEA').sum()
                else:
                    this_wNOFEA_count = 0
                if this_label == 5 or this_label == 6: # only sacrum or L5
                    this_wSPEC_count = df_warnings.str.contains(r'wSPEC').sum()
                else:
                    this_wSPEC_count = 0
                
            # Parse the errors
            df_errors = df.get('errors'+suffix)
            if df_errors.isnull().all():
                this_eFRAG_count = 0
                this_eSYMM_count = 0
                this_eMINVOL_count = 0
                this_eMAXVOL_count = 0
            else:
                this_eFRAG_count = df_errors.str.contains(r'eFRAG').sum()
                this_eMINVOL_count = df_errors.str.contains(r'eMINVOL').sum()
                this_eMAXVOL_count = df_errors.str.contains(r'eMAXVOL').sum()
                if this_label == 1 or this_label == 2 or this_label == 3 or this_label == 4: # only femur or pelvis
                    this_eSYMM_count = df_errors.str.contains(r'eSYMM').sum()
                else:
                    this_eSYMM_count = 0
            
            volume_mean = df.loc[:,'total_vol'+suffix].mean()
            volume_std = df.loc[:,'total_vol'+suffix].std()
            volume_min = df.loc[:,'total_vol'+suffix].min()
            volume_max = df.loc[:,'total_vol'+suffix].max()
            
            label_dict={
                'desc': this_desc,
                'label': this_label,
                'count': this_count,
                'volume_mean': volume_mean,
                'volume_std': volume_std,
                'volume_min': volume_min,
                'volume_max': volume_max,
                'status': this_status_pass,
                'percent_status': this_percent_status_pass,
                'eFRAG': this_eFRAG_count,
                'eSYMM': this_eSYMM_count,
                'eMINVOL': this_eMINVOL_count,
                'eMAXVOL': this_eMAXVOL_count,
                'wNOFEA': this_wNOFEA_count,
                'wSPEC': this_wSPEC_count                
            }
            
            results_dict[this_label]=label_dict
            
            if verbose:
                print('  {:>20s} {}'.format('desc',this_desc))
                print('  {:>20s} {}'.format('label',this_label))
                print('  {:>20s} {}'.format('count',this_count))
                print('  {:>20s} {:<10.1f}'.format('volume_mean',volume_mean))
                print('  {:>20s} {:<10.1f}'.format('volume_std',volume_std))
                print('  {:>20s} {:<10.1f}'.format('volume_min',volume_min))
                print('  {:>20s} {:<10.1f}'.format('volume_max',volume_max))
                print('  {:>20s} {}'.format('status pass',this_status_pass))
                print('  {:>20s} {:<10.1f}'.format('percent status pass',this_percent_status_pass))
                print('  {:>20s} {}'.format('eFRAG count',this_eFRAG_count))
                print('  {:>20s} {}'.format('eSYMM count',this_eSYMM_count))
                print('  {:>20s} {}'.format('eMINVOL count',this_eMINVOL_count))
                print('  {:>20s} {}'.format('eMAXVOL count',this_eMAXVOL_count))
                print('  {:>20s} {}'.format('wNOFEA count',this_wNOFEA_count))
                print('  {:>20s} {}'.format('wSPEC count',this_wSPEC_count))

    df_procrustes = df.loc[:,'procrustes']
    procrustes_count = df_procrustes.notnull().sum()
    procrustes_pass = df_procrustes.str.contains(r'PASS').sum()
    procrustes_fail = df_procrustes.str.contains(r'FAIL').sum()
    
    if procrustes_pass+procrustes_fail != procrustes_count:
        ogo.message('[ERROR] Procrustes number of PASS + FAIL != TOTAL.')
        os.sys.exit()
    
    results_dict['procrustes']={'desc': 'Procrustes',\
                                'count': procrustes_count,\
                                'pass': procrustes_pass,\
                                'fail': procrustes_fail}
                                
    if verbose:
        print('  {:>20s}'.format('-- procrustes --'))
        print('  {:>20s} {}'.format('count',procrustes_count))
        print('  {:>20s} {}'.format('pass',procrustes_pass))
        print('  {:>20s} {}'.format('fail',procrustes_fail))
        print('  {:>20s} {:.1f}'.format('percent pass',(procrustes_pass/procrustes_count)*100.0))
        
    df_labels_found = df.loc[:,'labels_found']
    labels_found_count = df_labels_found.notnull().sum()
    labels_found_pass = df_labels_found.str.contains(r'PASS').sum()
    labels_found_fail = df_labels_found.str.contains(r'FAIL').sum()
    
    if labels_found_pass+labels_found_fail != labels_found_count:
        ogo.message('[ERROR] Labels found test number of PASS + FAIL != TOTAL.')
        os.sys.exit()
    
    results_dict['labels_found']={'desc': 'Labels Found',\
                                  'count': labels_found_count,\
                                  'pass': labels_found_pass,\
                                  'fail': labels_found_fail}
    if verbose:
        print('  {:>20s}'.format('-- labels_found --'))
        print('  {:>20s} {}'.format('count',labels_found_count))
        print('  {:>20s} {}'.format('pass',labels_found_pass))
        print('  {:>20s} {}'.format('fail',labels_found_fail))
        print('  {:>20s} {:.1f}'.format('percent pass',(labels_found_pass/labels_found_count)*100.0))
        
    df_final = df.loc[:,'final']
    final_count = df_final.notnull().sum()
    final_pass = df_final.str.contains(r'PASS').sum()
    final_fail = df_final.str.contains(r'FAIL').sum()
    
    if final_pass+final_fail != final_count:
        ogo.message('[ERROR] Final testnumber of PASS + FAIL != TOTAL.')
        os.sys.exit()
    
    results_dict['final']={'desc': 'FINAL',\
                           'count': final_count,\
                           'pass': final_pass,\
                           'fail': final_fail}
    if verbose:
        print('  {:>20s}'.format('-- final --'))
        print('  {:>20s} {}'.format('count',final_count))
        print('  {:>20s} {}'.format('pass',final_pass))
        print('  {:>20s} {}'.format('fail',final_fail))
        print('  {:>20s} {:.1f}'.format('percent pass',(final_pass/final_count)*100.0))
        
    # CSV file
    if output_file:
        write.header = ['Label','','Pass','','Volume','','eMINVOL','','eMAXVOL','','eFRAG','','eSYMM','','wNOFEA','','wSPEC','']
        write(write.header,',',output_file)
        write.header = ['name','N','%','N','mean','std','%','N','%','N','%','N','%','N','%','N','%','N']
        write(write.header,',',output_file)
    
    # Print to screen
    print('{:=>20s}={:=>12s}={:=>12s}={:=>12s}={:=>12s}={:=>12s}={:=>12s}={:=>12s}={:=>12s}'.format(\
          '','','','','','','','',''))
    print('{:>20s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s}'.format(\
          'Label','Pass','Volume','eMINVOL','eMAXVOL','eFRAG','eSYMM','wNOFEA','wSPEC'))
    print('{:>20s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s}'.format(\
          'name (N)','% (N)','mean (std)','% (N)','% (N)','% (N)','% (N)','% (N)','% (N)'))
    print('{:->20s} {:->12s} {:->12s} {:->12s} {:->12s} {:->12s} {:->12s} {:->12s} {:->12s}'.format(\
          '-','-','-','-','-','-','-','-','-'))
    for k,v in results_dict.items():
        if k in list(lb.labels_dict.keys()):
            if v['count'] == 0:
                str_label = '{:s} ({:->4s})'.format(v['desc'],'')
                str_status = ''
                str_volume = ''
                str_eMINVOL = ''
                str_eMAXVOL = ''
                str_eFRAG = ''
                str_eSYMM = ''
                str_wNOFEA = ''
                str_wSPEC = ''
                
                entry = [v['desc'],'',\
                         '','',\
                         '','',\
                         '','',\
                         '','',\
                         '','',\
                         '','',\
                         '','',\
                         '','']
                
            else:
                
                entry = [v['desc'],v['count'],\
                         v['percent_status'],v['status'],\
                         v['volume_mean']/1000.0,v['volume_std']/1000.0,\
                         v['eMINVOL']/v['count']*100.0,v['eMINVOL'],\
                         v['eMAXVOL']/v['count']*100.0,v['eMAXVOL'],\
                         v['eFRAG']/v['count']*100.0,v['eFRAG']]
                
                str_label = '{:s} ({:4d})'.format(v['desc'],v['count'])
                str_status = '{:5.1f}({:5d})'.format(v['percent_status'],v['status'])
                str_volume = '{:5.1f}({:5.1f})'.format(v['volume_mean']/1000.0,v['volume_std']/1000.0)
                str_eMINVOL = '{:5.1f}({:5d})'.format(v['eMINVOL']/v['count']*100.0,v['eMINVOL'])
                str_eMAXVOL = '{:5.1f}({:5d})'.format(v['eMAXVOL']/v['count']*100.0,v['eMAXVOL'])
                str_eFRAG = '{:5.1f}({:5d})'.format(v['eFRAG']/v['count']*100.0,v['eFRAG'])

                if v['eSYMM'] == 0:
                    str_eSYMM = ''
                    entry.append('')
                    entry.append('')
                else:
                    str_eSYMM = '{:5.1f}({:5d})'.format(v['eSYMM']/v['count']*100.0,v['eSYMM'])
                    entry.append(v['eSYMM']/v['count']*100.0)
                    entry.append(v['eSYMM'])

                if v['wNOFEA'] == 0:
                    str_wNOFEA = ''
                    entry.append('')
                    entry.append('')
                else:
                    str_wNOFEA = '{:5.1f}({:5d})'.format(v['wNOFEA']/v['count']*100.0,v['wNOFEA'])
                    entry.append(v['wNOFEA']/v['count']*100.0)
                    entry.append(v['wNOFEA'])
                    
                if v['wSPEC'] == 0:
                    str_wSPEC = ''
                    entry.append('')
                    entry.append('')
                else:
                    str_wSPEC = '{:5.1f}({:5d})'.format(v['wSPEC']/v['count']*100.0,v['wSPEC'])
                    entry.append(v['wSPEC']/v['count']*100.0)
                    entry.append(v['wSPEC'])
                
                
            print('{:>20s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s}'.format(\
                  str_label,\
                  str_status,\
                  str_volume,\
                  str_eMINVOL,\
                  str_eMAXVOL,\
                  str_eFRAG,\
                  str_eSYMM,\
                  str_wNOFEA,\
                  str_wSPEC))
            
            if output_file:
                write(entry,',',output_file)

    print('{:->20s} {:->12s}'.format(\
          '-','-'))
          
    for k,v in results_dict.items():
        if k not in list(lb.labels_dict.keys()):
            str_label = '{:s} ({:>4d})'.format(v['desc'],v['count'])
            str_status = '{:5.1f}({:5d})'.format(v['pass']/v['count']*100.0,v['pass'])
            print('{:>20s} {:>10s}'.format(\
                  str_label,\
                  str_status))
            
            if output_file:
                entry = [v['desc'],v['count'],\
                         v['pass']/v['count']*100.0,v['pass']]
                write(entry,',',output_file)
            
            
    print('{:=>20s}={:=>12s}={:=>12s}={:=>12s}={:=>12s}={:=>12s}={:=>12s}={:=>12s}={:=>12s}'.format(\
          '','','','','','','','',''))

    print('Table caption:')
    print('  Pass rate is shown per bone for all the images.')
    print('  Pass rate of all bone segmentations is summarized as FINAL.')
    print('  Labels found means the image contained all expected labels.')
    print('  Procrustes analysis shows skeletons with unusual anatomy (e.g., L6 present).')
    print('  Errors for each bone include:')
    print('    eMINVOL –- exceeds min volume limit')
    print('    eMAXVOL –- exceeds max volume limit')
    print('    eFRAG   –- fragmented')
    print('    eSYMM   –- not symmetric')
    print('  Warnings for each bone include:')
    print('    wNOFEA  –- femur is cut off too much')
    print('    wSPEC   –- Pars defect for L5 or sacrum tip is disassociated')

    exit()
 
def main():
    # Setup description
    description = '''
    
Reads in a CSV file produced by running ogoValidate for a large batch job.
The data typically looks like this:

  line_hdr,name,desc,label,n_parts,total_vol,warnings,errors,status,desc,
  line_units,[text],[text],[#],[#],[mm3],[msg],[msg],[test],[text],[#],[#
  line_data,RETRO_00001.nii.gz,Femur Right,1,1,185516.1,,,PASS,Femur Left
  line_data,RETRO_00002.nii.gz,Femur Right,1,1,145373.6,wNOFEA ,,PASS,Fem
  line_data,RETRO_00003.nii.gz,Femur Right,1,1,114049.3,wNOFEA ,,PASS,Fem
  line_data,RETRO_00004.nii.gz,Femur Right,1,1,142556.2,wNOFEA ,,PASS,Fem

It is created by running scripts something like this:

  less RETRO_04000_nnUNet_3d_fullres_QA.txt | grep line_hdr >> results.csv
  less RETRO_04000_nnUNet_3d_fullres_QA.txt | grep line_units >> results.csv
  less RETRO_*_nnUNet_3d_fullres_QA.txt | grep line_data >> results.csv

This program will read the results.csv file and produce output to the screen
as well as a CSV summary table. You could read the CSV summary table into
an Excel file to make a nice looking output.

'''

    epilog = '''
Example call: 
     
SummarizeValidationData results_Dataset508_3d_fullres.csv
SummarizeValidationData results_Dataset508_3d_fullres.csv --output_file summary.csv

'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoSummarizeValidationData",
        description=description,
        epilog=epilog
    )

    parser.add_argument('input_file', metavar='CSV', help='Input CSV file (*.csv)')
    parser.add_argument('--output_file', default=None, metavar='CSV', 
                                        help='Summary output table (*.csv, default: %(default)s)')
    parser.add_argument('--overwrite', action='store_true', 
                                        help='Overwrite output file without asking')

    # Parse and display
    args = parser.parse_args()
    print(echo_arguments('SummarizeValidationData', vars(args)))

    # Run program
    SummarizeValidationData(**vars(args))


if __name__ == '__main__':
    main()
