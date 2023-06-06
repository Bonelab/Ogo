# /------------------------------------------------------------------------------+
# | 19-APR-2023                                                                  |
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
import datetime

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
# Various selection criteria for sorting to find the best image series
def try_by_modality(df):
    keywords = ['CT']
    
    for keyword in keywords:
        tmp_reduced_df = df.loc[(df['Modality'].str.contains(keyword, na=False))]
        
        if tmp_reduced_df.shape[0] > 0:
            return tmp_reduced_df, keyword, tmp_reduced_df.shape[0]
            
    tmp_reduced_df = df
    return tmp_reduced_df, 'NOT FOUND', tmp_reduced_df.shape[0]

def try_by_image_type(df):
    keywords = ['PRIMARY'] # If we can't get PRIMARY...
    
    for keyword in keywords:
        tmp_reduced_df = df.loc[(df['ImageType'].str.contains(keyword, na=False))]
        
        if tmp_reduced_df.shape[0] > 0:
            return tmp_reduced_df, keyword, tmp_reduced_df.shape[0]
            
    tmp_reduced_df = df
    return tmp_reduced_df, 'NOT FOUND', tmp_reduced_df.shape[0]

def try_by_study_description(df):
    keywords = ['ABD','PELVIS','BODY','CHEST'] # If we can't get ABDOMEN, then get PELVIS, then BODY, etc
    
    for keyword in keywords:
        tmp_reduced_df = df.loc[(df['StudyDescription'].str.contains(keyword, na=False))]
        
        if tmp_reduced_df.shape[0] > 0:
            return tmp_reduced_df, keyword, tmp_reduced_df.shape[0]
            
    tmp_reduced_df = df
    return tmp_reduced_df, 'NOT FOUND', tmp_reduced_df.shape[0]
    
def try_by_series_description(df):
    keywords = ['ABD','PEL','CHEST','Chest','AP','WB','Lung','LUNG','BODY','Body'] # If we can't get ABDOMEN, then get PELVIS, then CHEST, etc.
    
    for keyword in keywords:
        tmp_reduced_df = df.loc[(df['SeriesDescription'].str.contains(keyword, na=False))]
        
        if tmp_reduced_df.shape[0] > 0:
            return tmp_reduced_df, keyword, tmp_reduced_df.shape[0]
            
    tmp_reduced_df = df
    return tmp_reduced_df, 'NOT FOUND', tmp_reduced_df.shape[0]
    
def try_by_slice_thickness(df): # If we can't get 1.0mm, then get 2.0mm, etc.
    #slice_thicknesses = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
    slice_thicknesses = [3.0]
    
    for slice_thickness in slice_thicknesses:
        tmp_reduced_df = df.loc[(df['SliceThickness'] <= slice_thickness)]
        
        if tmp_reduced_df.shape[0] > 0:
            return tmp_reduced_df, slice_thickness, tmp_reduced_df.shape[0]
            
    tmp_reduced_df = df
    return tmp_reduced_df, 'NOT FOUND', tmp_reduced_df.shape[0]

def try_by_series_number(df): # Only accept series 
    acceptable_series_numbers = [50]
    
    for series_number in acceptable_series_numbers:
        tmp_reduced_df = df.loc[(df['SeriesNumber'] < series_number)]
        
        if tmp_reduced_df.shape[0] > 0:
            return tmp_reduced_df, series_number, tmp_reduced_df.shape[0]
            
    tmp_reduced_df = df
    return tmp_reduced_df, 'NOT FOUND', tmp_reduced_df.shape[0]

def try_by_number_of_images(df): # If we can't get 1000 slices, then get 900, etc.
    #images_in_series = [1000, 900, 800, 700, 600, 500, 400, 300, 250, 200, 150, 100, 50]
    images_in_series = [300, 200, 100]
    
    for n_images in images_in_series:
        tmp_reduced_df = df.loc[(df['NumberOfReferences'] >= n_images)]
        
        if tmp_reduced_df.shape[0] > 0:
            return tmp_reduced_df, n_images, tmp_reduced_df.shape[0]
            
    tmp_reduced_df = df
    return tmp_reduced_df, 'NOT FOUND', tmp_reduced_df.shape[0]

def try_by_pixel_spacing(df): # If we can't get 0.5mm, then get 0.6mm, etc.
    #pixel_spacings = ['0.5', '0.6', '0.7', '0.8', '0.9', '1.0', '1.5']
    pixel_spacings = ['1.0', '1.5']
    
    for pixel_spacing in pixel_spacings:
        tmp_reduced_df = df.loc[(df['PixelSpacing'].str.contains(pixel_spacing, na=False))]
        
        if tmp_reduced_df.shape[0] > 0:
            return tmp_reduced_df, pixel_spacing, tmp_reduced_df.shape[0]
            
    tmp_reduced_df = df
    return tmp_reduced_df, 'NOT FOUND', tmp_reduced_df.shape[0]

def exclude_by_series_description(df,keywords):
    
    for keyword in keywords:
        tmp_reduced_df = df.loc[(df['SeriesDescription'].str.contains(keyword, na=False, case=False) == False)]
        df = tmp_reduced_df
    
    return tmp_reduced_df, '({} keys)'.format(len(keywords)), tmp_reduced_df.shape[0]  
    
def exclude_by_number_of_images(df,n_images):
    
    tmp_reduced_df = df.loc[(df['NumberOfReferences'] >= n_images)]
        
    return tmp_reduced_df, '(min={})'.format(n_images), tmp_reduced_df.shape[0]

# +------------------------------------------------------------------------------+
# Sort scans
def DicomSelector(csvfile,output,overwrite):
    
    minNumberOfSlices = 30
    exclude_keywords = ['NECK','CTNKE','CTHDE','HEAD','LUNG','CHEST'] # Add LEGS, MPR, LIVER?
    
    # Check for valid input
    if os.path.splitext(csvfile)[1].lower() not in '.csv':
        ogo.message('[ERROR] Valid CSV file has not been provided. File provided:')
        ogo.message('        {}'.format(csvfile))
        os.sys.exit()
    
    # Check if output report exists and should overwrite
    if output:
        if os.path.isfile(output) and not overwrite:
            result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(output))
            if result.lower() not in ['y', 'yes']:
                ogo.message('Not overwriting. Exiting...')
                os.sys.exit()
    else:
        output = os.path.splitext(csvfile)[0] + '_selected.csv'
        if os.path.isfile(output) and not overwrite:
            result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(output))
            if result.lower() not in ['y', 'yes']:
                ogo.message('Not overwriting. Exiting...')
                os.sys.exit()
    output_pull = os.path.splitext(output)[0] + '.sh'

    ogo.message('Results will be output to {}'.format(output))

    ogo.message('Begin sorting...')
    
    debugging = False
    
    # Get all the names of exams into a sorted list
    csvData = pd.read_csv(csvfile)
    name_list = csvData['Name'].tolist()
    name_list = list(set(name_list))
    name_list.sort()
    
    # Sort data frame
    ogo.message('Put the data into order by Name, Date, Slices, Slice Thickness...')
    csvData.sort_values(["Name", "AcquisitionDate", "NumberOfReferences", "SliceThickness"], 
                                    axis=0,
                                    ascending=[True, True, False, True], 
                                    inplace=True)
    
    # Remove all series with <minNumberOfSlices slices
    ogo.message('Dataframe size is {}.'.format(csvData.shape[0]))
    ogo.message('')

    sorted_df = pd.DataFrame([]) # Used to collect the selected series
    
    # Begin sorting one name at a time
    for idx, name in enumerate(name_list):
        
        name_csvData = csvData.loc[csvData['Name'].str.contains(name,case=False)]
        
        date_list = name_csvData['AcquisitionDate'].tolist()
        date_list = list(set(date_list))
        date_list = [x for x in date_list if str(x) != 'nan'] # cleans date list of nan values
        date_list.sort()
        
        for this_date in date_list:
            reduced_df = name_csvData.loc[(name_csvData['AcquisitionDate'] == this_date)]
            n_series = reduced_df.shape[0]
            ogo.message('{:s} on {}'.format(name,this_date))
            ogo.message('  {:20s} {:>10s} {:>10d}'.format('Start',' ',n_series))
            
            # Exclude based on Number of Images
            if n_series > 1:
                reduced_df, val, n_series = exclude_by_number_of_images(reduced_df,minNumberOfSlices)
                ogo.message('  {:20s} {:>10} {:>10d}'.format('NumberOfImages',val,reduced_df.shape[0]))
            if debugging:
                print(reduced_df.loc[:,["SeriesNumber","NumberOfReferences","SeriesDescription","SliceThickness","PixelSpacing"]])
                  
            # Exclude based on Series Description
            if n_series > 1:
                reduced_df, val, n_series = exclude_by_series_description(reduced_df,exclude_keywords)
                ogo.message('  {:20s} {:>10} {:>10d}'.format('SeriesDescription',val,reduced_df.shape[0]))
            if debugging:
                print(reduced_df.loc[:,["SeriesNumber","NumberOfReferences","SeriesDescription","SliceThickness","PixelSpacing"]])
            
            # Include by Modality
            if n_series > 1:
                reduced_df, val, n_series = try_by_modality(reduced_df)
                ogo.message('  {:20s} {:>10s} {:>10d}'.format('Modality',val,reduced_df.shape[0]))
            if debugging:
                print(reduced_df.loc[:,["SeriesNumber","NumberOfReferences","SeriesDescription","SliceThickness","PixelSpacing"]])
             
            # Include by Slice Thickness
            if n_series > 1:
                reduced_df, val, n_series = try_by_slice_thickness(reduced_df)
                ogo.message('  {:20s} {:>10} {:>10d}'.format('SliceThickness',val,reduced_df.shape[0]))
            if debugging:
                print(reduced_df.loc[:,["SeriesNumber","NumberOfReferences","SeriesDescription","SliceThickness","PixelSpacing"]])
            
            # Pixel Spacing
            if n_series > 1:
                reduced_df, val, n_series = try_by_pixel_spacing(reduced_df)
                ogo.message('  {:20s} {:>10} {:>10d}'.format('PixelSpacing',val,reduced_df.shape[0]))
            if debugging:
                print(reduced_df.loc[:,["SeriesNumber","NumberOfReferences","SeriesDescription","SliceThickness","PixelSpacing"]])
            
            ogo.message('  {:20s} {:>10s} {:>10d}'.format('Finish',' ',n_series))
                  
              
            if False:
                # Study Description
                if n_series > 1:
                    reduced_df, val, n_series = try_by_study_description(reduced_df)
                    ogo.message('  {:20s} {:>10} {:>10d}'.format('StudyDescription',val,reduced_df.shape[0]))
                if debugging:
                   print(reduced_df.loc[:,["SeriesNumber","NumberOfReferences","SeriesDescription","SliceThickness","PixelSpacing"]])
            
                # Series Description
                if n_series > 1:
                    reduced_df, val, n_series = try_by_series_description(reduced_df)
                    ogo.message('  {:20s} {:>10} {:>10d}'.format('SeriesDescription',val,reduced_df.shape[0]))
                if debugging:
                   print(reduced_df.loc[:,["SeriesNumber","NumberOfReferences","SeriesDescription","SliceThickness","PixelSpacing"]])
            
                # Series Number
                if n_series > 1:
                    reduced_df, val, n_series = try_by_series_number(reduced_df)
                    ogo.message('  {:20s} {:>10} {:>10d}'.format('SeriesNumber',val,reduced_df.shape[0]))
                if debugging:
                    print(reduced_df.loc[:,["SeriesNumber","NumberOfReferences","SeriesDescription","SliceThickness","PixelSpacing"]])
                
                # Image Type
                if n_series > 1:
                    reduced_df, val, n_series = try_by_image_type(reduced_df)
                    ogo.message('  {:20s} {:>10} {:>10d}'.format('ImageType',val,reduced_df.shape[0]))
                if debugging:
                    print(reduced_df.loc[:,["SeriesNumber","NumberOfReferences","SeriesDescription","SliceThickness","PixelSpacing"]])
                
                # Number of Images
                if n_series > 1:
                    reduced_df, val, n_series = try_by_number_of_images(reduced_df)
                    ogo.message('  {:20s} {:>10} {:>10d}'.format('NumberOfImages',val,reduced_df.shape[0]))
                if debugging:
                   print(reduced_df.loc[:,["SeriesNumber","NumberOfReferences","SeriesDescription","SliceThickness","PixelSpacing"]])
                
                # Lowest Series Number
                if n_series > 1:
                    reduced_df = reduced_df.head(1)
                    ogo.message('  {:20s} {:>10} {:>10d}'.format('First series',val,reduced_df.shape[0]))
                    n_series = reduced_df.shape[0]
                    ogo.message('  {:20s} {:>10} {:>10d} **'.format('MORE CRITERIA NEEDED',' ',n_series))
                    #print(reduced_df.loc[:,["SeriesNumber","NumberOfReferences","SeriesDescription","SliceThickness","PixelSpacing"]])
                if debugging:
                    print(reduced_df.loc[:,["SeriesNumber","NumberOfReferences","SeriesDescription","SliceThickness","PixelSpacing"]])
                
                ogo.message('  {:16} (Series Number = {:3d}) {:>6d}'.format('FINAL SELECTION',int(reduced_df.iloc[0]['Series Number']),reduced_df.shape[0]))
                print(reduced_df.loc[:,["SeriesNumber","NumberOfReferences","SeriesDescription","SliceThickness","PixelSpacing"]])
            
            ogo.message('')

            # Collect the selected series
            sorted_df = pd.concat([sorted_df,reduced_df])
    
    # Final selection
    ogo.message('Dataframe size is reduced to {}.'.format(sorted_df.shape[0]))
    ogo.message('')
    
    # Write csv file
    ogo.message('Writing csv output file:')
    ogo.message('  {}'.format(output))
    sorted_df.to_csv(output)
    ogo.message('')

    # Write script to pull the dicom files
    ogo.message('Writing file for pulling dicom files:')
    ogo.message('  {}'.format(output_pull))
    ogo.message('')

    
    fpull = open(output_pull, 'w')
    sorted_df = sorted_df.reset_index()  # make sure indexes pair with number of rows
    current_dir = os.getcwd()

    fpull.write('#!/bin/bash\n')
    fpull.write('exit\n')
    fpull.write('# \n')
    fpull.write('# Generated by ogoDicomScanner: {}\n'.format(datetime.datetime.now()))
    fpull.write('# \n')
    fpull.write('# Script to pull a dicom series, create a NIfTI, and visualize with offscreen rendering.\n')
    fpull.write('# Prior to running, remove line 2 after checking the following:\n')
    fpull.write('#   1. Check BASE_DIR, OUTPUT_DICOM_DIR, OUTPUT_NIFTI_DIR for accuracy.\n')
    fpull.write('#   2. Check that conda is installed on the path described below.\n')
    fpull.write('#   3. Check that ogo is installed in your conda environment.\n')
    fpull.write('#   4. Remove line #2 once all checks are complete.\n')
    fpull.write(os.linesep)
    fpull.write('BASE_DIR=\'{}\'\n'.format(current_dir))
    fpull.write('OUTPUT_DICOM_DIR=\'{}\'\n'.format(current_dir+'/dicom'))
    fpull.write('OUTPUT_NIFTI_DIR=\'{}\'\n'.format(current_dir+'/nifti'))
    fpull.write(os.linesep)
    fpull.write('source /Users/skboyd/opt/miniconda3/etc/profile.d/conda.sh\n')
    fpull.write('conda activate ogo\n')
    fpull.write(os.linesep)

    for index, row in sorted_df.iterrows():
        fname = '{}_{}_{}_{}'.format( \
            row['Name'], \
            str(row['StudyDate']).replace('-',''), \
            str(row['SeriesDescription']).replace(',','_').replace(' ','_').replace('__','_'), \
            int(row['SeriesNumber']) \
            )
        dicompull_cmd = 'dicompull -k SeriesInstanceUID={} {}/\'{}\' -o {}\'/{}\'\n'.format( \
            row['SeriesInstanceUID'], \
            '${BASE_DIR}', \
            os.path.split(row['ReferencedFileID'])[0], \
            '${OUTPUT_DICOM_DIR}', \
            fname \
            )
        dicomtonifti_cmd = 'dicomtonifti -brz --fsl {}/\'{}\' -o {}\'/{}.nii.gz\'\n'.format( \
            '${OUTPUT_DICOM_DIR}', \
            fname, \
            '${OUTPUT_NIFTI_DIR}', \
            fname \
            )
        ogoVisualize_cmd = 'ogoVisualize vis2d --offscreen {}/\'{}.nii.gz\' --outfile {}/\'{}.tif\'\n'.format( \
            '${OUTPUT_NIFTI_DIR}', \
            fname, \
            '${OUTPUT_NIFTI_DIR}', \
            fname \
            )
        fpull.write(dicompull_cmd)
        fpull.write(dicomtonifti_cmd)
        fpull.write(ogoVisualize_cmd)
        fpull.write(os.linesep)
        
    fpull.write('exit\n')
    fpull.close()
    
    ogo.message('Done!')
    
def main():
    # Setup description
    description = '''
Takes a CSV file of meta data where each row represents a dicom series within
dicom studies. The selection criteria is used to narrow the search for the 
best series. 

'''

    epilog = '''
Example call: 
     
ogoDicomSelector list.csv --output narrowed_list.csv

'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoDicomSelector",
        description=description,
        epilog=epilog
    )

    parser.add_argument('csvfile', metavar='CSV',
                                        help='Input CSV file (*.csv)')
    parser.add_argument('--output', default=None, metavar='CSV', 
                                        help='Reduced list of dicom series (*.csv, default: %(default)s)')
    parser.add_argument('--overwrite', action='store_true', 
                                        help='Overwrite output file without asking')

    # Parse and display
    args = parser.parse_args()
    print(echo_arguments('DicomSelector', vars(args)))

    # Run program
    DicomSelector(**vars(args))


if __name__ == '__main__':
    main()
