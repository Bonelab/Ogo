#####
# ogo_bmd_tabulate.py
#
# This script collects all the BMD results from the specified directory.
#
#####
#
# Andrew Michalski
# University of Calgary
# Biomedical Engineering Graduate Program
# September 6, 2019
# Modified to Py3: March 25, 2020
#####

script_version = 1.0

##
# Import the required modules
import os
import sys
import glob
import argparse
import time
import pandas as pd
from datetime import date
from collections import OrderedDict


from ogo.util.echo_arguments import echo_arguments


def bmdTabulate(args):

    ##
    # Collect the input arguments
    directory = args.analysis_directory
    bone = args.bone

    ##
    # Read the files in the Directory
    files = []

    if bone == 1:
        bone_filename ='*_RT_FEMUR_BMD_Results.txt'
        bone_filename2 ='_RT_FEMUR_BMD_Results.txt'
        fileName = "RT_FEMUR_BMD_Results.txt"
        bmd_type = 'RT FEMUR BMD (mg/cc)'
    elif bone == 2:
        bone_filename = '*_LT_FEMUR_BMD_Results.txt'
        bone_filename2 = '_LT_FEMUR_BMD_Results.txt'
        fileName = "LT_FEMUR_BMD_Results.txt"
        bmd_type = 'LT FEMUR BMD (mg/cc)'
    elif bone == 4:
        bone_filename = '*_L4_BMD_Results.txt'
        bone_filename2 = '_L4_BMD_Results.txt'
        fileName = "L4_SPINE_BMD_Results.txt"
        bmd_type = 'L4 SPINE BMD (mg/cc)'
    else:
        print("Bone value set is not defined. Ending script.")
        sys.exit()

    ##
    # Create dataframe for output
    df = pd.DataFrame(columns = ['ID', bmd_type])

    ##
    # Read and extract the data
    os.chdir(directory)
    files = sorted(glob.glob(bone_filename))
    k = 0
    for i in files:
        # open each txt file
        data = pd.read_csv(i, sep = "\t", header = None)
        bmd = data.loc[10][1]
        ID = i.replace(bone_filename2, "")
        df.loc[k] = [ID, bmd]
        k = k + 1

    ##
    # Write Output TXT file of dataframe
    df.to_csv(fileName, sep = '\t', index = False, header = True)

    print("Script Complete.")


def main():
    description = '''
    This script collects all the BMD results from the specified directory. 
    
    INPUT: Data_Directory; 
    OUTPUT: TXT file of results
    
    '''


    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoBMDTabulate",
        description=description
    )
    

    parser.add_argument("analysis_directory",
    help = "Path/To/Data/Directory")

    parser.add_argument("--bone", type = int,
    default = 1,
    help = "Set the bone value to analyze. 1 = RT_FEMUR, 2 = LT_FEMUR, 4 = L4_SPINE. (Default: %(default)s)")


    # Parse and display
    args = parser.parse_args()
    
    print(echo_arguments('BMD_Tabulate', vars(args)))

    # Run program
    bmdTabulate(args)

if __name__ == '__main__':
    main()