# /------------------------------------------------------------------------------+
# | 8-FEB-2022                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+

# Imports
import argparse
import os
from batchgenerators.utilities.file_and_folder_operations import *
import numpy
#from ogo.util.echo_arguments import echo_arguments
#import ogo.util.Helper as ogo
from collections import OrderedDict

def print_dict(level,max_level,col1,col2,col3,key,value_in):
    txt = ''
    value_type, value_style = get_value_type_and_style(value_in)  
    txt += indent(level)+'{message1: <{width1}}{message2: >{width2}}\n'.format(\
        message2=str(key)+':', width1=col1, \
        message1=value_type,   width2=col2)
        
    level += 1
    if level > max_level:
        return txt
    
    for key,value in value_in.items():
        value_type, value_style = get_value_type_and_style(value)
        
        if (value_style == 'dict'):
            txt += print_dict(level,max_level,col1,col2,col3,key,value)
        elif (value_style == 'list'):
            txt += print_list(level,max_level,col1,col2,col3,key,value)
        elif (value_style == 'val'):
            txt += print_value(level,max_level,col1,col2,col3,key,value)
        else:
            print('key [{}]= {}\n'.format(value_type,key))
            print('value [{}] = {}\n'.format(value_type,value_in))
            os.sys.exit('[ERROR] Unknown data type.')

    return txt

def print_list(level,max_level,col1,col2,col3,key,value):
    txt = ''
    value_type, value_style = get_value_type_and_style(value)  
    if level > max_level:
        return txt
    txt += indent(level)+'{message1: <{width1}}{message2: >{width2}}\n'.format(\
        message2=str(key)+':', width1=col1, \
        message1=value_type,   width2=col2)
    for list_item in value:
        list_item_type, list_item_style = get_value_type_and_style(list_item)

        if (list_item_style == 'dict'): # Sometimes elements of a list are DICT
            txt += print_dict(level,max_level,col1,col2,col3,key,list_item)
        else:
            txt += indent(level)+'{message: >{width}}\n'.format(message=str(list_item), width=col1+col2+col3)
            
    return txt

def print_value(level,max_level,col1,col2,col3,key,value):
    txt = ''
    value_type, value_style = get_value_type_and_style(value)
    if level > max_level:
        return txt
    txt += indent(level)+'{message1: <{width1}}{message2: >{width2}} {message3: <{width3}}\n'.format(\
        message2=str(key)+':', width1=col1, \
        message1=value_type,   width2=col2, \
        message3=str(value), width3=col3)
    return txt

def get_value_type_and_style(value):
    
    # style dict
    if (type(value) is dict):
        return 'DICT','dict'
    elif (type(value) is OrderedDict):
        return 'ODICT','dict'
    
    # style list
    elif (type(value) is list):
        return 'LIST','list'
    elif (type(value) is numpy.ndarray):
        return 'NDARRAY','list'
    elif (type(value) is tuple):
        return 'TUPLE','list'

    # style val
    elif (type(value) is int):
        return 'INT','val'
    elif (type(value) is float):
        return 'FLOAT','val'
    elif (type(value) is str):
        return 'STR','val'
    elif (type(value) is type(None)):
        return 'NONE','val'
    elif (type(value) is bool):
        return 'BOOL','val'
    elif (type(value) is numpy.bool_):
        return 'NBOOL','val'
    elif (type(value) is numpy.str_):
        return 'NSTR','val'
    elif (type(value) is numpy.float32):
        return 'F32','val'
    elif (type(value) is numpy.float64):
        return 'F64','val'
    elif (type(value) is numpy.int64):
        return 'I64','val'
    
    else:
        print('[ERROR] Unexpected type: ',type(value))
        exit()

def indent(level):
    return ' '*level*4+'> '
            
def ReadPickle(input,max_level):

    if not input.lower().endswith('.pkl'):
        os.sys.exit('[ERROR] Input must be type PKL file: \"{}\"'.format(input))
    
    if max_level<0:
        os.sys.exit('[ERROR] A level less than 0 not allowed.')
        
    col1 = 12
    col2 = 30
    col3 = 10
    level = 0
    
    print('> Reading: {}'.format(input))
    print('-'*80)

    pkl = load_pickle(input)
    txt = ''        
    
    if type(pkl) is list:
        txt += print_list(level,max_level,col1,col2,col3,'list',pkl)
    elif (type(pkl) is dict) or (type(pkl) is OrderedDict):
        for key, value in pkl.items():
            value_type, value_style = get_value_type_and_style(value)
            
            #print('key [{}]= {}\n'.format(value_type,key))
            #print('value [{}] = {}\n'.format(value_type,value))
            
            if (value_style == 'dict'):
                txt += print_dict(level,max_level,col1,col2,col3,key,value)
            elif (value_style == 'list'):
                txt += print_list(level,max_level,col1,col2,col3,key,value)
            elif (value_style == 'val'):
                txt += print_value(level,max_level,col1,col2,col3,key,value)
            else:
                print('key [{}]= {}\n'.format(value_type,key))
                print('value [{}] = {}\n'.format(value_type,value))
                os.sys.exit('[ERROR] Unknown data type.')
    else:
        os.sys.exit('[ERROR] Unknown PKL input type ({}).'.format(type(pkl)))
        
    print(txt)
    print('-'*80)
    
def main():
    # Setup description
    description='''
Reads a pickle file.

If you need to minimize how much detail you see, select
a lower max_level (minimum is 0).

'''
    epilog='''
Example calls: 
ogoReadPickle plans.pkl
ogoReadPickle --max_level 1 plans.pkl

'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoReadPickle",
        description=description,
        epilog=epilog
    )
    parser.add_argument('input', help='Input pickle file (*.pkl)')
    parser.add_argument('--max_level', type=int, default=99, metavar='LEVEL', help='Set level limit (default: %(default)s)')

    print()

    # Parse and display
    args = parser.parse_args()
    #print(echo_arguments('ReadPickle', vars(args)))

    # Run program
    ReadPickle(**vars(args))

if __name__ == '__main__':
    main()