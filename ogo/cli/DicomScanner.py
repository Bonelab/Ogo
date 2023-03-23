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

# Create dictionary of required tags
def get_dicom_tags():  
    meta_dict = {
        # General meta data
        '0008|0008': 'Image Type',
        '0008|0060': 'Modality',
        '0008|0022': 'Acquisition Date',
        '0008|1030': 'Study Description',
        '0008|103E': 'Series Description', # added
        '0020|0011': 'Series Number', # added
        '0008|0050': 'Accession Number', # added
        '0054|0081': 'Number of Slices', # added
        '0010|21b0': 'Additional Patient History', # added
        '0008|0070': 'Manufacturer',
        '0008|0100': 'Code Value',
        '0008|103E': 'Series Description',
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
        '0008|1050': 'Patient\'s Name 1', # added
        '0008|1060': 'Patient\'s Name 2', # added
    
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
def DicomScanner(input_dir, output_csv, deliminator, selection, overwrite):

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
    
    # Use a FIFO queue to process all files
    # We assume no one has made loops using symlinks...
    
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
                ogo.message('Adding directory to queue\n\n     {}\n'.format(full_path))
                dir_queue.put(full_path) # A list of *all* sub-directories in input_directory (input1)
                continue
                
            #n_dir = dir_queue.qsize()
            #ogo.message('Number of directories on queue: {}'.format(n_dir))
    
            # Work with a DICOM image (all directories have been identified; working on files)
            if not dcm_processed and os.path.isfile(full_path):
                
                ogo.message('Processing\n\n     {}\n'.format(full_path))
            
                # Catch cases where file is not DICOM (e.g. desktop.ini)
                try:
                    sitk.ReadImage(full_path)
                except:
                    continue
                    
                # Create the DICOM reader and find all series within DICOM directory
                reader = sitk.ImageSeriesReader()
                series_ids = reader.GetGDCMSeriesIDs(this_dir)
                n_series_ids = len(series_ids)
                
                # Loop through each DICOM series
                for idx,this_series in enumerate(series_ids):
                    
                    ogo.message('')
                    ogo.message('Working on series {} of {}'.format(idx,n_series_ids))
                    
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
                    reader.MetaDataDictionaryArrayUpdateOn()
                    reader.LoadPrivateTagsOn()
                    try:
                        reader.Execute()
                    except:
                        print('[ERROR] Problem reading: {}'.format(full_path))
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
                            print(' ** NO META KEYS for {}, series {}'.format(this_name,this_series))
                    
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
Traverses DICOM directories from a given root path and finds all exams and
associated series. The output is a CSV file that can be read with tools such
as Excel.

A selection tool helps narrow down the useful DICOM within a given exam based
on fixed criteria from the meta data.
'''

    epilog = '''
Example call: 
     
ogoDicomScanner ./root_dicom_data_dir --output_csv ./processed_headers.csv
ogoDicomScanner /Users/skboyd/Library/CloudStorage/OneDrive-UniversityofCalgary/ML/data/CTDXAICI/dicom/2023-03-06/data/controls
ogoDicomScanner /Users/skboyd/Library/CloudStorage/OneDrive-UniversityofCalgary/ML/data/CTDXAICI/dicom/2023-03-06/data/controls/CTDXAICI_0030

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
