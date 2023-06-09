# /------------------------------------------------------------------------------+
# | 23-MAR-2023                                                                  |
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

# +------------------------------------------------------------------------------+
# Create dictionary of required tags
def get_dicom_tags():  
    meta_dict = {
        # General meta data
        '0008|0008': 'Image Type',
        '0008|0060': 'Modality',
        '0008|0022': 'Acquisition Date',
        '0008|1030': 'Study Description',
        '0008|103e': 'Series Description', # added
        '0020|000e': 'Series Instance UID', # added
        '0020|0011': 'Series Number', # added
        '0008|0050': 'Accession Number', # added
        '0054|0081': 'Number of Slices', # added
        '0010|21b0': 'Additional Patient History', # added
        '0008|0070': 'Manufacturer',
        '0008|0100': 'Code Value',
        '0008|0102': 'Coding Scheme Designator',
        '0008|0104': 'Code Meaning',
    
        # Scanner info
        '0008|1090': 'Manufacturer\'s Model Name',
        '0009|0010': 'Scanner private 1',
        '0009|1001': 'Scanner private 2',
        '0009|1002': 'Scanner private 3',
        '0009|1004': 'Scanner private 4',
        '0009|1027': 'Scanner private 5',
    
        # Patient Data
        '0010|0010': 'Patient\'s Name',
        '0010|0020': 'Patient\'s ID',
        '0010|0030': 'Patient\'s Birth Date',
        '0010|0040': 'Patient\'s Sex',
        '0010|1010': 'Patient\'s Age',
        '0008|1050': 'Performing Physician\'s Name', # added
        '0008|1060': 'Name of Physician(s) Reading Study', # added
    
        # Scan options
        '0018|0050': 'Slice Thickness',
        '0018|0060': 'kVp',
        '0018|0088': 'Spacing Between Slices',
        '0018|0090': 'Data Collection Diameter',
        '0018|1000': 'Device Serial Number',
        '0018|1020': 'Software Versions(s)',
        '0018|1030': 'Protocal Name',
        '0018|1100': 'Reconstruction Diameter',
        '0018|1110': 'Distance Source to Detector',
        '0018|1111': 'Distance Source to Patient',
        '0018|1120': 'Gantry/Detector Tilt',
        '0018|1130': 'Table Height',
        '0018|1140': 'Rotation Direction',
        '0018|1150': 'Exposure Time',
        '0018|1151': 'X-ray Tube Current',
        '0018|1152': 'Exposure',
        '0018|1160': 'Filter Type',
        '0018|1170': 'Generator Power',
        '0018|1190': 'Focal Spot(s)',
        '0018|1210': 'Convolution Kernel',
        '0018|5100': 'Patient Position',
        '0018|9305': 'Revolution Time',
        '0018|9036': 'Single Collimation Width',
        '0018|9037': 'Total Collimation Width',
        '0018|9309': 'Table Speed',
        '0018|9310': 'Table Feed per Rotation',
        '0018|9311': 'Spiral Pitch Factor',
    
        # Positioning
        '0020|0032': 'Image Position (Patient)',
        '0020|0037': 'Image Orientation (Patient)',
        '0020|0052': 'Frame of Reference UID',
        '0020|1040': 'Position Reference Indicator',
        '0020|1041': 'Slice Location',
    
        # Pixel Data
        '0028|0002': 'Samples per Pixel',
        '0028|0004': 'Photometric Interpretation',
        '0028|0010': 'Rows',
        '0028|0011': 'Columns',
        '0028|0030': 'Pixel Spacing',
        '0028|0100': 'Bits Allocated',
        '0028|0101': 'Bits Stored',
        '0028|0102': 'High Bit',
        '0028|0103': 'Pixel Representation',
        '0028|0120': 'Pixel Padding Value',
        '0028|1050': 'Window Center',
        '0028|1051': 'Window Width',
        '0028|1052': 'Rescale Intercept',
        '0028|1053': 'Rescale Slope',
        '0028|1054': 'Rescale Type',
    }
    return collections.OrderedDict(meta_dict)

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
def try_by_image_type(df):
    keywords = ['PRIMARY'] # If we can't get PRIMARY...
    
    for keyword in keywords:
        tmp_reduced_df = df.loc[(df['Image Type'].str.contains(keyword))]
        
        if tmp_reduced_df.shape[0] > 0:
            return tmp_reduced_df, keyword, tmp_reduced_df.shape[0]
            
    tmp_reduced_df = df
    return tmp_reduced_df, 'NOT FOUND', tmp_reduced_df.shape[0]

def try_by_study_description(df):
    keywords = ['ABD','PELVIS','BODY','CHEST'] # If we can't get ABDOMEN, then get PELVIS, then BODY, etc
    
    for keyword in keywords:
        tmp_reduced_df = df.loc[(df['Study Description'].str.contains(keyword))]
        
        if tmp_reduced_df.shape[0] > 0:
            return tmp_reduced_df, keyword, tmp_reduced_df.shape[0]
            
    tmp_reduced_df = df
    return tmp_reduced_df, 'NOT FOUND', tmp_reduced_df.shape[0]
    
def try_by_series_description(df):
    keywords = ['ABD','PEL','CHEST','Chest','AP','WB','Lung','LUNG','BODY','Body'] # If we can't get ABDOMEN, then get PELVIS, then CHEST, etc.
    
    for keyword in keywords:
        tmp_reduced_df = df.loc[(df['Series Description'].str.contains(keyword, na=False))]
        
        if tmp_reduced_df.shape[0] > 0:
            return tmp_reduced_df, keyword, tmp_reduced_df.shape[0]
            
    tmp_reduced_df = df
    return tmp_reduced_df, 'NOT FOUND', tmp_reduced_df.shape[0]
    
def try_by_slice_thickness(df): # If we can't get 1.0mm, then get 2.0mm, etc.
    #slice_thicknesses = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
    slice_thicknesses = [1.0, 1.5, 2.0, 2.5, 3.0]
    
    for slice_thickness in slice_thicknesses:
        tmp_reduced_df = df.loc[(df['Slice Thickness'] <= slice_thickness)]
        
        if tmp_reduced_df.shape[0] > 0:
            return tmp_reduced_df, slice_thickness, tmp_reduced_df.shape[0]
            
    tmp_reduced_df = df
    return tmp_reduced_df, 'NOT FOUND', tmp_reduced_df.shape[0]

def try_by_series_number(df): # Only accept series 
    acceptable_series_numbers = [50]
    
    for series_number in acceptable_series_numbers:
        tmp_reduced_df = df.loc[(df['Series Number'] < series_number)]
        
        if tmp_reduced_df.shape[0] > 0:
            return tmp_reduced_df, series_number, tmp_reduced_df.shape[0]
            
    tmp_reduced_df = df
    return tmp_reduced_df, 'NOT FOUND', tmp_reduced_df.shape[0]

def try_by_number_of_images(df): # If we can't get 1000 slices, then get 900, etc.
    #images_in_series = [1000, 900, 800, 700, 600, 500, 400, 300, 250, 200, 150, 100, 50]
    images_in_series = [200, 100]
    
    for n_images in images_in_series:
        tmp_reduced_df = df.loc[(df['Number of Images'] >= n_images)]
        
        if tmp_reduced_df.shape[0] > 0:
            return tmp_reduced_df, n_images, tmp_reduced_df.shape[0]
            
    tmp_reduced_df = df
    return tmp_reduced_df, 'NOT FOUND', tmp_reduced_df.shape[0]

def try_by_pixel_spacing(df): # If we can't get 1.0mm, then get 2.0mm, etc.
    pixel_spacings = ['0.5', '0.6', '0.7', '0.8', '0.9', '1.0', '1.5']
    
    for pixel_spacing in pixel_spacings:
        tmp_reduced_df = df.loc[(df['Pixel Spacing'].str.contains(pixel_spacing))]
        
        if tmp_reduced_df.shape[0] > 0:
            return tmp_reduced_df, pixel_spacing, tmp_reduced_df.shape[0]
            
    tmp_reduced_df = df
    return tmp_reduced_df, 'NOT FOUND', tmp_reduced_df.shape[0]

# +------------------------------------------------------------------------------+
# Sort scans
def sort_scans(csvfile,input_dir,overwrite):
    
    # Check for valid input
    if not csvfile or (os.path.splitext(csvfile)[1].lower() not in '.csv'):
        ogo.message('[ERROR] Valid CSV file has not been provided. File provided:')
        ogo.message('        {}'.format(csvfile))
        os.sys.exit()
    
    # Define output file
    csvfile_sorted = os.path.splitext(csvfile)[0] + '_sorted.csv'
    c3dfile_sorted = os.path.splitext(csvfile)[0] + '_sorted.sh'
    if os.path.isfile(csvfile_sorted) and not overwrite:
        result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(csvfile_sorted))
        if result.lower() not in ['y', 'yes']:
            ogo.message('Not overwriting. Exiting...')
            os.sys.exit()
    if os.path.isfile(c3dfile_sorted) and not overwrite:
        result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(c3dfile_sorted))
        if result.lower() not in ['y', 'yes']:
            ogo.message('Not overwriting. Exiting...')
            os.sys.exit()

    ogo.message('Begin sorting...')
    
    debugging = False
    
    # Get all the names of exams into a sorted list
    csvData = pd.read_csv(csvfile)
    name_list = csvData['Name'].tolist()
    name_list = list(set(name_list))
    name_list.sort()
    
    # Sort data frame
    ogo.message('Put the data into order by Name, Date, Slices, Slice Thickness...')
    csvData.sort_values(["Name", "Acquisition Date", "Number of Images", "Slice Thickness"], 
                                    axis=0,
                                    ascending=[True, True, False, True], 
                                    inplace=True)
    # Remove and series with <5 slices
    ogo.message('Dataframe size is {}.'.format(csvData.shape[0]))
    csvData = csvData.loc[(csvData['Number of Images'] > 5)]
    ogo.message('Dataframe size is reduced to {} for series with more than 5 slices.'.format(csvData.shape[0]))

    sorted_df = pd.DataFrame([]) # Used to collect the selected series
    
    # Begin sorting one name at a time
    for idx, name in enumerate(name_list):
        
        name_csvData = csvData.loc[csvData['Name'].str.contains(name,case=False)]
        
        date_list = name_csvData['Acquisition Date'].tolist()
        date_list = list(set(date_list))
        date_list.sort()
        date_list = [int(x) for x in date_list if str(x) != 'nan'] # cleans date list of nan values
        
        for this_date in date_list:
            reduced_df = name_csvData.loc[(name_csvData['Acquisition Date'] == this_date)]
            n_series = reduced_df.shape[0]
            ogo.message('{} on {} has {} series'.format(name,this_date,n_series))
            
            # Accept only low Series Numbers that avoid dose reports, etc.
            if n_series > 1:
                reduced_df, val, n_series = try_by_series_number(reduced_df)
                ogo.message('  [{:>20s} == {:>10}] {:>8d}'.format('Series Number',val,reduced_df.shape[0]))
                
            # Select only series that are primary (not reformats)
            if n_series > 1:
                reduced_df, val, n_series = try_by_image_type(reduced_df)
                ogo.message('  [{:>20s} == {:>10s}] {:>8d}'.format('Image Type',val,reduced_df.shape[0]))
            
            if debugging:
                print(reduced_df.loc[:,["Series Number","Number of Images","Series Description","Slice Thickness","Pixel Spacing"]])
           
            # Use series description to identify
            if n_series > 1:
                reduced_df, val, n_series = try_by_series_description(reduced_df)
                ogo.message('  [{:>20s} == {:>10s}] {:>8d}'.format('Series Description',val,reduced_df.shape[0]))

            if debugging:
                print(reduced_df.loc[:,["Series Number","Number of Images","Series Description","Slice Thickness","Pixel Spacing"]])
            
            # Use slice thickness to identify highest resolution
            if n_series > 1:
                reduced_df, val, n_series = try_by_slice_thickness(reduced_df)
                ogo.message('  [{:>20s} == {:>10}] {:>8d}'.format('Slice Thickness',val,reduced_df.shape[0]))
            
            if debugging:
                print(reduced_df.loc[:,["Series Number","Number of Images","Series Description","Slice Thickness","Pixel Spacing"]])
            
            # Use number of images to identify largest scan
            if n_series > 1:
                reduced_df, val, n_series = try_by_number_of_images(reduced_df)
                ogo.message('  [{:>20s} == {:>10}] {:>8d}'.format('Number of Images',val,reduced_df.shape[0]))
            
            if debugging:
                print(reduced_df.loc[:,["Series Number","Number of Images","Series Description","Slice Thickness","Pixel Spacing"]])
            
            # User pixel spacing to select
            if n_series > 1:
                reduced_df, val, n_series = try_by_pixel_spacing(reduced_df)
                ogo.message('  [{:>20s} == {:>10}] {:>8d}'.format('Pixel Spacing',val,reduced_df.shape[0]))
            
            if debugging:
                print(reduced_df.loc[:,["Series Number","Number of Images","Series Description","Slice Thickness","Pixel Spacing"]])
            
            # At this point we can't see much difference based on meta data, so we pick lowest series number
            if n_series > 1:
                reduced_df = reduced_df.head(1)
                ogo.message('  [{:>20s} == {:>10}] {:>8d}'.format('First series',n_series,reduced_df.shape[0]))
                n_series = reduced_df.shape[0]
                #ogo.message('  {:25} {:>19d} **'.format('MORE CRITERIA NEEDED',n_series))
                #print(reduced_df.loc[:,["Series Number","Number of Images","Study Description","Slice Thickness","Pixel Spacing"]])
                
            if debugging:
                print(reduced_df.loc[:,["Series Number","Number of Images","Series Description","Slice Thickness","Pixel Spacing"]])
            
            #print(reduced_df.iloc[0]['Series Number'])
            ogo.message('  {:16} (Series Number = {:3d}) {:>6d}'.format('FINAL SELECTION',int(reduced_df.iloc[0]['Series Number']),reduced_df.shape[0]))
            print(reduced_df.loc[:,["Series Description","Number of Images","Slice Thickness","Pixel Spacing"]])
            ogo.message('')
            
            # Collect the selected series
            sorted_df = pd.concat([sorted_df,reduced_df])
    
    # Write csv file
    ogo.message('Writing csv output file:')
    ogo.message('  {}'.format(csvfile_sorted))
    sorted_df.to_csv(csvfile_sorted)

    # Write file to execute c3d (convert3d)
    ogo.message('Writing file for c3d execution:')
    ogo.message('  {}'.format(c3dfile_sorted))
    
    output = open(c3dfile_sorted, 'w')
    sorted_df = sorted_df.reset_index()  # make sure indexes pair with number of rows

    output.write('#!/bin/bash\n')
    output.write('BASE_DIR=\'{}\'\n'.format(input_dir))
    output.write('OUTPUT_DIR=\'{}\'\n'.format(input_dir))
    output.write(os.linesep)
    output.write('# NEXT STEPS:\n')
    output.write('#   Run this script that calls c3d to \'grep\' the Series UID. The problem is that the simplitk\n')
    output.write('#   class doesn\'t capture the entire UID needed to convert the series. When you run the script\n')
    output.write('#   the results will provide you the full Series UID, which can then be substituted into the\n')
    output.write('#   commented out lines below. Remember to use the TAB as a search key when manipulating the\n')
    output.write('#   outputs into a list. Also, note that some images capture more than one GREP output (usually\n')
    output.write('#   the series in the 600 range). You\'ll have to remove those before copy and paste. Boy, this\n')
    output.write('#   seems unnecessarily complex!\n')
    output.write(os.linesep)
    for index, row in sorted_df.iterrows():
        txt = '# c3d -dicom-series-read {}\'{}\' {} -o {}\'/{}_{:d}.nii\' # Series number {}\n'.format( \
            '${BASE_DIR}', \
            row['Directory'].replace(input_dir,''), \
            row['Series'], \
            '${OUTPUT_DIR}', \
            row['Name'], \
            int(row['Acquisition Date']), \
            int(row['Series Number']) \
            )
        output.write(txt)
    output.write(os.linesep)
    for index, row in sorted_df.iterrows():
        txt = 'c3d -dicom-series-list {}\'{}\' | grep {} >> get_correct_series_uid.txt\n'.format( \
            '${BASE_DIR}', \
            row['Directory'].replace(input_dir,''), \
            row['Series'], \
            )
        output.write(txt)
        txt_splitter = 'echo \'#--split--\' >> get_correct_series_uid.txt\n'
        output.write(txt_splitter)
    output.write('exit\n')
    output.close()
    
    ogo.message('Done!')
            

# +------------------------------------------------------------------------------+
def DicomScanner(input_dir, output_csv, deliminator, selection, overwrite):

    # Assumes a CSV output has already been created
    if selection:
        ogo.message('Performing selection of best series from available exams')
        sort_scans(output_csv,input_dir,overwrite)
        exit()
        
    # Check if output report exists and should overwrite
    if output_csv:
        if os.path.isfile(output_csv) and not overwrite:
            result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(output_csv))
            if result.lower() not in ['y', 'yes']:
                ogo.message('Not overwriting. Exiting...')
                os.sys.exit()
    
    # Check input directory
    if not os.path.isdir(input_dir):
        ogo.message('[ERROR] Input directory must be a directory.')
        os.sys.exit()
        
    meta_dict =  get_dicom_tags()

    # Write header
    write.header = ['Directory', 'Series', 'Name', 'Number of Images']
    for key, value in meta_dict.items():
        write.header.append(value)
    write(write.header,deliminator,output_csv)
    
    dir_queue = queue.Queue()
    dir_queue.put(os.path.abspath(input_dir))
    
    ogo.message('Starting with directory\n\n     {}\n'.format(input_dir))
    
    while not dir_queue.empty():
        # Pop a directory
        this_dir = dir_queue.get()

        #ogo.message('Working on directory:\n         {}'.format(this_dir))

        dcm_processed = False
        for name in os.listdir(this_dir): # finds all dirs in this_dir

            # Check if directory
            full_path = os.path.abspath(os.path.join(this_dir, name))
            if os.path.isdir(full_path):
                # Push onto queue
                ogo.message('Adding directory to queue\n\n     ../{}\n'.format(full_path.replace(input_dir,'')))
                dir_queue.put(full_path) # A list of *all* sub-directories in input_directory (input1)
                continue
                            
            # Work with a DICOM image (all directories have been identified; working on files)
            if not dcm_processed and os.path.isfile(full_path):
                
                ogo.message('Processing\n\n     ../{}\n'.format(full_path.replace(input_dir,'')))
            
                # Catch cases where file is not DICOM (e.g. desktop.ini)
                try:
                    sitk.ReadImage(full_path)
                except:
                    continue
                    
                # Create the DICOM reader and find all series within DICOM directory
                reader = sitk.ImageSeriesReader()
                reader.LoadPrivateTagsOn()
                reader.MetaDataDictionaryArrayUpdateOn()
                
                #useSeriesDetails = True
                #series_ids = reader.GetGDCMSeriesIDs(this_dir,useSeriesDetails)
                series_ids = reader.GetGDCMSeriesIDs(this_dir)
                n_series_ids = len(series_ids)
                
                # Loop through each DICOM series
                for idx,this_series in enumerate(series_ids):
                    
                    ogo.message('')
                    ogo.message('Working on series {} of {}'.format(idx+1,n_series_ids))
                    
                    dicom_names = reader.GetGDCMSeriesFileNames(this_dir,this_series)
                    n_dicom_names = len(dicom_names)
                    
                    # Skip DICOM series that are clearly not useful
                    if n_dicom_names == 0:
                        continue
                    
                    # Parse name, DICOM directory
                    this_dicom_dir = this_dir.split(os.path.sep)[-1]
                    
                    # Flip switch at last series
                    if (idx+1)==n_series_ids:
                        dcm_processed = True
                    ogo.message('Processing {} DICOM images'.format(len(dicom_names)))
                    ogo.message('      dir: {} '.format(this_dicom_dir))
                    ogo.message('   series: {} '.format(this_series))
                
                    # Read meta data
                    reader.SetFileNames(dicom_names)
                    try:
                        reader.Execute()
                    except:
                        ogo.message('[ERROR] Problem reading: {}'.format(full_path))
                        continue
                        
                    #this_name = re.match(r'.*\/(RETRO_\d+)\/.*', this_dir)[1]
                    #this_name = re.match(r'.*\/(CTDXAICI_\d+_V\d+)\/.*', this_dir)[1]
                    this_name = re.match(r'.*\/(CTDXAICI_\d+)\/.*', this_dir)[1] # Help for regular expressions: https://regexr.com
                    entry = [this_dir, this_series, this_name, n_dicom_names]
                    
                    # Read meta data in DICOM
                    if False:
                        try:
                            all_data_keys = reader.GetMetaDataKeys(5)
                            for k in all_data_keys:
                                print(k,reader.GetMetaData(5,k))
                        except:
                            ogo.message('[WARNING] NO META KEYS for {}, series {}'.format(this_name,this_series))
                    
                    # Parse meta data
                    for key, value in meta_dict.items():
                        this_entry = ''
                        if len(dicom_names) > 5 and reader.HasMetaDataKey(5,key):
                            this_entry = str(reader.GetMetaData(5,key)).replace(',',' ') # reader.GetMetaDataKeys() will get all keys
                        entry.append(this_entry)
                    
                    write(entry,deliminator,output_csv)

    ogo.message('Done.')
    
    
def main():
    # Setup description
    description = '''
    
DO NOT USE. Use OgoDicomScanner instead. S Boyd, June 7, 2023

Traverses DICOM directories from a given root path and finds all exams and
associated series. The output is a CSV file of DICOM meta data that can be 
read with tools such as Excel.

Once a CSV file is created, use the --selection option to narrow down the 
list of DICOM series that are useful for quantitative analysis.

The last step is to convert selected series (based on Series Instance UID)
to NIfTI format. This is done using third-party convert3D, developed as part
of the ITK-Snap project by ITK, and requires downloading and installing.

'''

    epilog = '''
Example call: 
     
ogoDicomScanner ./root_dicom_data_dir --output_csv ./processed_headers.csv

'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoDicomScanner",
        description=description,
        epilog=epilog
    )

    parser.add_argument('input_dir', 
                                        help='Input directory to traverse')
    parser.add_argument('--output_csv', default=None, metavar='CSV', 
                                        help='Rows of DICOM meta data per series (*.csv, default: %(default)s)')
    parser.add_argument('--deliminator', type=str, default=',', metavar='DELIM',
                                        help='Output CSV deliminator (default: %(default)s)')
    parser.add_argument('--selection', action='store_true', 
                                        help='Select DICOM series by criteria (default: %(default)s)')
    parser.add_argument('--overwrite', action='store_true', 
                                        help='Overwrite output file without asking')

    # Parse and display
    args = parser.parse_args()
    print(echo_arguments('DicomScanner', vars(args)))

    # Run program
    DicomScanner(**vars(args))


if __name__ == '__main__':
    main()
