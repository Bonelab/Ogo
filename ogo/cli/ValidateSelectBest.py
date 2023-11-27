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

def check_dimensions_match(input_files):
    for idx,f in enumerate(input_files):
        with open(f, 'r') as file:
            result_dict = yaml.safe_load(file)
        fname = result_dict['input_parameters']['input_image']
        ct = sitk.ReadImage(fname, sitk.sitkUInt8)
        dim = ct.GetSize()
        ogo.message('  {}: {}'.format(f,dim))
        if idx==0:
            dim_base = ct.GetSize()
            fname_base = fname
        else:
            dim = ct.GetSize()
            if dim != dim_base:
                ogo.message('[ERROR] Input images do not have the same dimensions.')
                os.sys.exit()

def is_finished(mydict):
    finished = True
    for k,v in mydict.items():
        if v['found'] == False:
            finished = False
    #print('Finished = {}'.format(finished))
    return finished

def sort_input_files_by_procrustes(input_files):
    tmp_dict={}
    new_input_files=[]
    for f in input_files:
        with open(f, 'r') as file:
            result_dict = yaml.safe_load(file)
            check_dict = result_dict.get('check')
            tmp_dict[f]=check_dict['procrustes']['disparity']
            
    sorted_dict = dict(sorted(tmp_dict.items(), key=lambda item: item[1]))
    
    for k,v in sorted_dict.items():
        new_input_files.append(k)
    return new_input_files
    
def ValidateSelectBest(input_files,output_file,list_of,status_false,frankenstein,overwrite):

    best_dict = {
        1: {'desc': lb.labels_dict[1]['LABEL'], 'source': '', 'found': False},
        2: {'desc': lb.labels_dict[2]['LABEL'], 'source': '', 'found': False},
        3: {'desc': lb.labels_dict[3]['LABEL'], 'source': '', 'found': False},
        4: {'desc': lb.labels_dict[4]['LABEL'], 'source': '', 'found': False},
        5: {'desc': lb.labels_dict[5]['LABEL'], 'source': '', 'found': False},
        6: {'desc': lb.labels_dict[6]['LABEL'], 'source': '', 'found': False},
        7: {'desc': lb.labels_dict[7]['LABEL'], 'source': '', 'found': False},
        8: {'desc': lb.labels_dict[8]['LABEL'], 'source': '', 'found': False},
        9: {'desc': lb.labels_dict[9]['LABEL'], 'source': '', 'found': False},
       10: {'desc': lb.labels_dict[10]['LABEL'],'source': '', 'found': False},
       11: {'desc': lb.labels_dict[11]['LABEL'],'source': '', 'found': False}
       }
    
    results_list = []
    
    # Set up to read image and mask inputs
    for ifile in input_files:
        ogo.pass_check_if_file_exists(ifile)
        ogo.pass_check_file_ending(ifile,['.yaml'])
    
    # We find lists of files with status as requested
    if not frankenstein:
        ogo.message('Reading {} input files.'.format(len(input_files)))
        say_snip=True
        
        for idx,f in enumerate(input_files):
            if (idx<10 or idx==len(input_files)-1):
                ogo.message('Reading {}'.format(f))
            elif say_snip:
                ogo.message(' ...<snip>... ')
                say_snip=False
                
            with open(f, 'r') as file:
                result_dict = yaml.safe_load(file)
        
            check_dict = result_dict.get('check')
            labels_dict = check_dict.get('labels')
            input_parameters_dict = result_dict.get('input_parameters')
            
            # If a specific bone label, or 
            if 'label' in list_of:
                label = int(list_of.replace('label_',''))
                if labels_dict[label]['status'] == status_false:
                    results_list.append({'file': input_parameters_dict['input_image'], 'yaml': f})
            
            # If a general label like 'intact_spine', etc
            else:
                if check_dict[list_of]['status'] == status_false:
                    results_list.append({'file': input_parameters_dict['input_image'], 'yaml': f})
                    #if list_of == 'procrustes':
                    #    print('{:9s}{} --> {:8.6f}'.format('','procrustes',check_dict[list_of]['disparity']))
                    
        ogo.message('')
        ogo.message('Looking for a list of \'{}\' with status \'{}\'.'.format(list_of,status_false))
        
        n_found = len(results_list)
        n_total = len(input_files)
        ogo.message('  found {} of {}: {:.3f}%'.format(n_found,n_total,n_found/n_total*100.0))
        ogo.message('')
        if results_list:
            ogo.message('List output:')
            for elements in results_list:
                print('{:9s}{} --> {}'.format('',elements['file'],elements['yaml']))
                
    # Frankenstein method starts here
    if frankenstein:
        
        # We start with the lowest procrustes score as it is most likely the best
        input_files = sort_input_files_by_procrustes(input_files)
        
        # First pass to find 'intact_left_femur', intact_right_femur', 'intact_spine', 'final'
        for f in input_files:
            
            ogo.message('Reading {}'.format(f))
            with open(f, 'r') as file:
                result_dict = yaml.safe_load(file)
            
            check_dict = result_dict.get('check')
            labels_dict = check_dict.get('labels')
            input_parameters_dict = result_dict.get('input_parameters')
            
            ogo.message('  procrustes = {:8.6f}'.format(check_dict['procrustes']['disparity']))

            if check_dict['final']['status']:
                best_dict[1]['source']  = f
                best_dict[2]['source']  = f
                best_dict[3]['source']  = f
                best_dict[4]['source']  = f
                best_dict[5]['source']  = f
                best_dict[6]['source']  = f
                best_dict[7]['source']  = f
                best_dict[8]['source']  = f
                best_dict[9]['source']  = f
                best_dict[10]['source'] = f
                best_dict[11]['source'] = f
                
                best_dict[1]['found']  = True
                best_dict[2]['found']  = True
                best_dict[3]['found']  = True
                best_dict[4]['found']  = True
                best_dict[5]['found']  = True
                best_dict[6]['found']  = True
                best_dict[7]['found']  = True
                best_dict[8]['found']  = True
                best_dict[9]['found']  = True
                best_dict[10]['found'] = True
                best_dict[11]['found'] = True
                ogo.message('  --> taking \'{}\' labels from {}'.format('all',f))
            
            if is_finished(best_dict):
                break
            
            if check_dict['intact_spine']['status'] \
                    and not best_dict[7]['found'] \
                    and not best_dict[8]['found'] \
                    and not best_dict[9]['found'] \
                    and not best_dict[10]['found']:
                best_dict[7]['source']  = f
                best_dict[8]['source']  = f
                best_dict[9]['source']  = f
                best_dict[10]['source'] = f
                best_dict[7]['found']  = True
                best_dict[8]['found']  = True
                best_dict[9]['found']  = True
                best_dict[10]['found'] = True
                ogo.message('  --> taking \'{}\' labels from {}'.format('intact_spine',f))
            
            if check_dict['intact_right_femur']['status'] and not best_dict[1]['found']:
                best_dict[1]['source']  = f
                best_dict[1]['found']  = True
                ogo.message('  --> taking \'{}\' label from {}'.format('intact_right_femur',f))
                
            if check_dict['intact_left_femur']['status'] and not best_dict[2]['found']:
                best_dict[2]['source']  = f
                best_dict[2]['found']  = True
                ogo.message('  --> taking \'{}\' label from {}'.format('intact_left_femur',f))
        
            # Go bone by bone now
            for k,v in best_dict.items():
                if labels_dict[k]['status'] and not best_dict[k]['found']:
                    if k==11:
                        if labels_dict[k]['n_parts']==1:
                            best_dict[k]['source'] = f
                            best_dict[k]['found'] = True
                            ogo.message('  --> taking \'{}\' ({}) from {}'.format(best_dict[k]['desc'],k,f))
                    else:
                        best_dict[k]['source'] = f
                        best_dict[k]['found'] = True
                        ogo.message('  --> taking \'{}\' ({}) from {}'.format(best_dict[k]['desc'],k,f))
    
        ogo.message('')
        ogo.message('Frankenstein result:')
        for k,v in best_dict.items():
            print('{:>20s}: {:s}'.format(v['desc'],v['source']))
        ogo.message('')
        
        # Assemble output image
        if output_file:
            
            # Check if output exists and should overwrite
            ogo.pass_check_if_output_exists(output_file,overwrite)
            ogo.pass_check_file_ending(output_file,['.nii.gz'])

            ogo.message('Assembling output image')
            ogo.message('Checking that input images have the same dimensions')
            check_dimensions_match(input_files)
            ogo.message('    --> pass')
            ogo.message('')
            
            ct_base = None
            
            # Cycle through the input files
            for f in input_files:
                ct_file_opened = False
                for k,v in best_dict.items():
                    if v['source']==f: # if this input file will contribute
                        if not ct_file_opened:
                            with open(f, 'r') as file:
                                result_dict = yaml.safe_load(file)
                            ct_name = result_dict.get('input_parameters')['input_image']
                            ct = sitk.ReadImage(ct_name, sitk.sitkUInt8)
                            ct_file_opened = True
                        if ct_base is None:
                            ct_base = ct<0
                        if v['found']:
                            label = k
                            ogo.message('  gathering label {} ({}) from {}'.format(v['desc'],label,ct_name))
                            #print('yep: ' + f + ct_name)
                            bin_part = ct==label
                            mask = 1 - bin_part
                            ct_base = sitk.Mask(ct_base, mask)
                            ct_base = ct_base + label*bin_part
            
            ogo.message('')
            ogo.message('Writing output file:')
            ogo.message('  {}'.format(output_file))
            sitk.WriteImage(ct_base, output_file)            
            
        else:
            ogo.message('No output file defined, so no output written.')
        
    ogo.message('Done.')
     
def main():
    # Setup description
    description = '''
    
Reads in a list of YAML files representing inference of the same CT dataset but 
using different ML models. It can be used to perform a number of tasks:

1. Generate a list of model names based on a criteria (e.g. final status is pass, 
   intact_spine status is pass, intact_right_femur is pass, etc)

2. Create a new inference result by a Frankenstein of multiple inferences.
   This can be done by explicit instructions or automatically.
             
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
    parser.add_argument('--output_file', default=None, metavar='NIfTI',help='Best result (*.nii.gz, default: %(default)s)')
    parser.add_argument('--list_of', default='final', choices=['procrustes', 'all_found', 'intact_right_femur', 'intact_left_femur', 'intact_spine', 'intact_all_labels', 'final', 'label_1', 'label_2', 'label_3', 'label_4', 'label_5', 'label_6', 'label_7', 'label_8', 'label_9', 'label_10', 'label_11'],
                                help='Create a list (default: %(default)s)')
    parser.add_argument('--status_false', action='store_false',help='List conditions when status is false')
    parser.add_argument('--frankenstein', action='store_true',help='Assemble inference with best labels (default: %(default)s)')
    parser.add_argument('--overwrite', action='store_true',help='Overwrite output file without asking')

    # Parse and display
    args = parser.parse_args()

    # Run program
    ValidateSelectBest(**vars(args))


if __name__ == '__main__':
    main()
