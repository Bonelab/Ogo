#####
# ogo_fe_pistoia_tabulate.py
#
# This script collects all the FE Pistoia results from the specified directory.
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

##
def fePistoiaTab(args):

    parser = argparse.ArgumentParser(
        description="""This script collects all the FE Pistoia results from the specified directory. INPUT: Data_Directory; OUTPUT: TXT file of results""")

    parser.add_argument("analysis_directory",
        help = "Path/To/Data/Directory")

    parser.add_argument("--model", type = int,
        default = 1,
        help = "Set the bone value to analyze. 1 = RT_FEMUR_SF;\n2 = RT_FEMUR_SLS;\n3 = LT_FEMUR_SF;\n4 = LT_FEMUR_SLS;\n6 = L4_SPINE_VC.\n(Default: %(default)s)")

    ##
    # Collect the input arguments
    args = parser.parse_args()
    directory = args.analysis_directory
    model = args.model

    ##
    # Read the files in the Directory
    files = []

    if model == 1:
        model_filename ='*_RT_FEMUR_SF_PISTOIA.txt'
        model_filename2 ='_RT_FEMUR_SF_PISTOIA.txt'
        fileName = "RT_FEMUR_SF_PISTOIA_Results.txt"
        model_type = 'RT FEMUR SF Failure Load (N)'
    elif model == 2:
        model_filename = '*_RT_FEMUR_SLS_PISTOIA.txt'
        model_filename2 = '_RT_FEMUR_SLS_PISTOIA.txt'
        fileName = "RT_FEMUR_SLS_PISTOIA_Results.txt"
        model_type = 'RT FEMUR SLS Failure Load (N)'
    elif model == 3:
        model_filename ='*_LT_FEMUR_SF_PISTOIA.txt'
        model_filename2 ='_LT_FEMUR_SF_PISTOIA.txt'
        fileName = "LT_FEMUR_SF_PISTOIA_Results.txt"
        model_type = 'LT FEMUR SF Failure Load (N)'
    elif model == 4:
        model_filename = '*_LT_FEMUR_SLS_PISTOIA.txt'
        model_filename2 = '_LT_FEMUR_SLS_PISTOIA.txt'
        fileName = "LT_FEMUR_SLS_PISTOIA_Results.txt"
        model_type = 'LT FEMUR SLS Failure Load (N)'
    elif model == 6:
        model_filename = '*_L4_FE_PISTOIA.txt'
        model_filename2 = '_L4_FE_PISTOIA.txt'
        fileName = "L4_SPINE_FE_PISTOIA_Results.txt"
        model_type = 'L4 SPINE Failure Load (N)'
    else:
        print("Model value set is not defined. Ending script.")
        sys.exit()

    ##
    # Create dataframe for output
    df = pd.DataFrame(columns = ['ID', model_type])

    ##
    # Read and extract the data
    os.chdir(directory)
    files = sorted(glob.glob(model_filename))
    k = 0
    for i in files:
        ID = i.replace(model_filename2, "")

        # Check to see if the file is empty and add to talbe if needed
        if os.stat(i).st_size == 0:
            x_fl = ""
            y_fl = ""
            z_fl = ""
            if model == 1:
                df.loc[k] = [ID, y_fl]
            if model == 2:
                df.loc[k] = [ID, z_fl]
            if model == 3:
                df.loc[k] = [ID, y_fl]
            if model == 4:
                df.loc[k] = [ID, z_fl]
            if model == 6:
                df.loc[k] = [ID, z_fl]
            k = k + 1
            continue

        # open each txt file and extract the failure loads
        lines = [line.rstrip('\n') for line in open(i)]
        lines = lines[11].replace(" ", "")
        lines = lines.replace("Failureload(RF*factor):", "")
        index_1 = lines.find("E")
        index_2 = lines.find("E", index_1 + 1)
        index_3 = lines.find("E", index_2 + 1)
        indices = [index_1, index_2, index_3]
        x_fl = lines[0:index_1+4]
        y_fl = lines[index_1+4:index_2+4]
        z_fl = lines[index_2+4:]


        if model == 1:
            df.loc[k] = [ID, y_fl]
        if model == 2:
            df.loc[k] = [ID, z_fl]
        if model == 3:
            df.loc[k] = [ID, y_fl]
        if model == 4:
            df.loc[k] = [ID, z_fl]
        if model == 6:
            df.loc[k] = [ID, z_fl]
        k = k + 1

    ##
    # Write Output TXT file of dataframe
    df.to_csv(fileName, sep = '\t', index = False, header = True)

    print("Script Complete.")


def main():
    description = '''
    This script collects all the FE Pistoia results from the specified directory. 
    
    INPUT: Data_Directory; 
    OUTPUT: TXT file of results
    
    '''


    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoFePistoiaTabulate",
        description=description
    )

    parser.add_argument("analysis_directory",
        help = "Path/To/Data/Directory")

    parser.add_argument("--model", type = int,
        default = 1,
        help = "Set the bone value to analyze. 1 = RT_FEMUR_SF;\n2 = RT_FEMUR_SLS;\n3 = LT_FEMUR_SF;\n4 = LT_FEMUR_SLS;\n6 = L4_SPINE_VC.\n(Default: %(default)s)")


    # Parse and display
    args = parser.parse_args()
    
    print(echo_arguments('fe_Pistoia_Tabulate', vars(args)))

    # Run program
    fePistoiaTab(args)

if __name__ == '__main__':
    main()