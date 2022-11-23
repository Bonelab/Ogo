# /------------------------------------------------------------------------------+
# | 22-NOV-2022                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+

script_version = 1.0

# Imports
import argparse
import os
import sys
import dicom2nifti
from datetime import date
from collections import OrderedDict
from ogo.util.echo_arguments import echo_arguments
import ogo.util.Helper as ogo

# +------------------------------------------------------------------------------+
def dcm2nii(dicom_directory, output_folder, report_only, overwrite):
    
    report_file = output_folder + '/report.txt'
    
    # Check if report exists and should overwrite
    if os.path.isfile(report_file) and not overwrite:
        result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(report_file))
        if result.lower() not in ['y', 'yes']:
            ogo.message('Not overwriting. Exiting...')
            os.sys.exit()
    
    # Check output folder
    if not os.path.isdir(dicom_directory):
        os.sys.exit('[ERROR] DICOM folder does not exist: \"{}\"'.format(dicom_directory))
    
    # Read the DICOM directory
    ogo.message('Reading DICOM directory:')
    ogo.message('  {}'.format(dicom_directory))
    dicom_dir_data = dicom2nifti.common.read_dicom_directory(dicom_directory, stop_before_pixels=True)
    
    # Scan all the images
    ogo.message('Scanning the DICOM images for meta data.')
    series_dict = OrderedDict()
    count = 0
    for FileDataset in dicom_dir_data:
        count += 1
        SeriesDescription = FileDataset.data_element('SeriesDescription') # returns a pydicom.dataelem.DataElement
        StudyID = FileDataset.data_element('StudyID')
        SeriesNumber = FileDataset.data_element('SeriesNumber')
        InstanceNumber = FileDataset.data_element('InstanceNumber')
        Rows = FileDataset.data_element('Rows')
        Columns = FileDataset.data_element('Columns')
        PixelSpacing = FileDataset.data_element('PixelSpacing')
        #SliceThickness = FileDataset.data_element('SliceThickness')
        SpacingBetweenSlices = FileDataset.data_element('SpacingBetweenSlices')
        StudyInstanceUID = FileDataset.data_element('StudyInstanceUID')
        SeriesInstanceUID = FileDataset.data_element('SeriesInstanceUID')
        
        sn = SeriesNumber.value
        if sn not in series_dict.keys():
            series_dict[sn] = {
                'SeriesDescription':    SeriesDescription.value,
                'StudyID':              StudyID.value,
                'SeriesNumber':         SeriesNumber.value,
                'InstanceNumber':       InstanceNumber.value,
                'Rows':                 Rows.value,
                'Columns':              Columns.value,
                'Slices':               1,
                'PixelSpacing':         PixelSpacing.value,
                'SpacingBetweenSlices': SpacingBetweenSlices.value,
                'StudyInstanceUID':     StudyInstanceUID.value,
                'SeriesInstanceUID':    SeriesInstanceUID.value
                }
        else:
            series_dict[sn]['Slices'] += 1
        
    ogo.message('Examined {} images belonging to {} image series.'.format(len(dicom_dir_data),len(series_dict.keys())) )
    
    # Build the report
    ogo.message('Build the report.')
    series_dict = OrderedDict(sorted(series_dict.items(), key=lambda x: x[1]['SeriesNumber']))
    el = list(series_dict.keys())[0]
    v = series_dict.get(el)
    cmds = ''
    report = ''
    report += '{:30s} {}\n'.format('Directory',dicom_directory)
    report += '{:30s} {}\n'.format('Instance Number',v['InstanceNumber'])
    report += '{:30s} {}\n'.format('Study Instance UID',v['StudyInstanceUID'])
    report += '\n'
    report += '{:15s} {:30s} {:20s} {:10s} {:20s} {:20s}\n'.format('Series Number',\
                                                            'Description',\
                                                            'Dimensions',\
                                                            'Images',\
                                                            'Pixel Spacing',\
                                                            'Filename')
    
    for k,v in series_dict.items():
        dims = '{} x {} x {}'.format(v['Rows'],v['Columns'],v['Slices'])
        el_size = '{:.2f} x {:.2f} x {:.2f}'.format(v['PixelSpacing'][0],v['PixelSpacing'][1],v['SpacingBetweenSlices'])
        # Filename
        series_filename = str(v['SeriesNumber']) + \
                          '_' + \
                          v['SeriesDescription'].lower().replace(' ','_').replace('.','').replace('/','') + \
                          '.nii.gz'
        # Highlight
        if ('axial' in v['SeriesDescription'].lower() or 'torso' in v['SeriesDescription'].lower()) and \
            v['Slices'] > 100 and \
            v['SpacingBetweenSlices'] < 3.0:
            highlight = '*'
            cmds += 'ls ' + output_folder + '/' + series_filename + '\n'
        else:
            highlight = ' '

        report += '{:<10d}{:5s} {:30s} {:20s} {:<10d} {:20s} {}\n'.format(v['SeriesNumber'],\
                                                              highlight,\
                                                              v['SeriesDescription'],\
                                                              dims,\
                                                              v['Slices'],\
                                                              el_size,\
                                                              series_filename)
    
    report += '\n'
    report += 'Helpful commands:\n'
    report += cmds
    
    # Save report
    if report_file:
        ogo.message('Saving report to file:')
        ogo.message('  {}'.format(report_file))
        txt_file = open(report_file, "w")
        txt_file.write(report)
        txt_file.close()

    # Generate the NIFTI files
    if not report_only:
        ogo.message('Generating NIFTI files:')
        ogo.message('  {}'.format(output_folder))
        dicom2nifti.convert_directory(dicom_directory, output_folder, compression=True, reorient=True)
    else:
        ogo.message('No NIFTI generated.')
    
    ogo.message('Done.')
    
def main():
    # Setup description
    description = '''
Utility that converts all series in a DICOM exam into NIFTI files.

The directory where .dcm files are located is provided as input and
this utility will scan all images and collate them into individual
NIFTI files for each image series.

A report is generated that provides basic characteristics of all the
image series. This report can be useful to decide which image series
to use.

Usually there are more image series than needed. The report is intended
to help decide which files to utilize.

'''

    epilog = '''
This utility is based on the package dicom2nifti:
https://dicom2nifti.readthedocs.io/en/latest/index.html

Example call:
ogodcm2nii /input/Patient0008 /nifti/Patient0008 \\
           --overwrite

'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="dcm2nii",
        description=description,
        epilog=epilog
    )
    parser.add_argument('dicom_directory', default=".", metavar='DICOM_DIR', help='Directory of DCM files')
    parser.add_argument('output_folder', default=".", metavar='NIFTI_DIR', help='Directory to write NIFTI files')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite output report without asking')
    parser.add_argument('--report_only', action='store_true', help='Generate report without writing NIFTI (default: %(default)s)')

    print()

    args = parser.parse_args()
    print(echo_arguments('dcm2nii', vars(args)))

    # Run program
    dcm2nii(**vars(args))


if __name__ == '__main__':
    main()
