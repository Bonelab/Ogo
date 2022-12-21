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
import dicom2nifti.settings as settings

# +------------------------------------------------------------------------------+
def dcm2nii(dicom_directory, output_folder, report_only, skip_report, overwrite):
    
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
    
    if not skip_report:
        
        # Read the DICOM directory
        ogo.message('Reading DICOM directory:')
        ogo.message('  {}'.format(dicom_directory))
        dicom_dir_data = dicom2nifti.common.read_dicom_directory(dicom_directory, stop_before_pixels=True)
        
        # Scan all the images
        ogo.message('Scanning the DICOM images for meta data.')
        series_dict = OrderedDict()
        for FileDataset in dicom_dir_data:

            try:
                SeriesDescription = FileDataset.data_element('SeriesDescription').value # returns a pydicom.dataelem.DataElement
            except KeyError:
                ogo.message('WARNING: Could not find SeriesDescription')
                SeriesDescription = 'n/a'
                
            try:
                StudyID = FileDataset.data_element('StudyID').value
            except KeyError:
                ogo.message('WARNING: Could not find StudyID')
                StudyID = 0
            
            try:
                PatientBirthDate = FileDataset.data_element('PatientBirthDate').value
            except KeyError:
                ogo.message('WARNING: Could not find PatientBirthDate')
                PatientBirthDate = 0
               
            try:
                PatientID = FileDataset.data_element('PatientID').value
            except KeyError:
                ogo.message('WARNING: Could not find PatientID')
                PatientID = 'n/a'
            
            try:
                PatientName = FileDataset.data_element('PatientName').value
            except KeyError:
                ogo.message('WARNING: Could not find PatientName')
                PatientName = 'n/a'
            
            try:
                SeriesNumber = FileDataset.data_element('SeriesNumber').value
            except KeyError:
                ogo.message('WARNING: Could not find SeriesNumber')
                SeriesNumber = 0
            
            try:
                InstanceNumber = FileDataset.data_element('InstanceNumber').value
            except KeyError:
                ogo.message('WARNING: Could not find InstanceNumber')
                InstanceNumber = 0

            try:
                AccessionNumber = FileDataset.data_element('AccessionNumber').value
            except KeyError:
                ogo.message('WARNING: Could not find AccessionNumber')
                AccessionNumber = 0

            try:
                Rows = FileDataset.data_element('Rows').value
            except KeyError:
                ogo.message('WARNING: Could not find Rows')
                Rows = 0
                
            try:
                Columns = FileDataset.data_element('Columns').value
            except KeyError:
                ogo.message('WARNING: Could not find Columns')
                Columns = 0
            
            try:
                SliceThickness = FileDataset.data_element('SliceThickness').value
            except KeyError:
                ogo.message('WARNING: Could not find SliceThickness')
                SliceThickness = 0
            
            try:
                SpacingBetweenSlices = FileDataset.data_element('SpacingBetweenSlices').value
            except KeyError:
                ogo.message('WARNING: Could not find SpacingBetweenSlices')
                SpacingBetweenSlices = 0
            
            try:
                StudyInstanceUID = FileDataset.data_element('StudyInstanceUID').value
            except KeyError:
                ogo.message('WARNING: Could not find StudyInstanceUID')
                StudyInstanceUID = 'n/a'
            
            try:
                SeriesInstanceUID = FileDataset.data_element('SeriesInstanceUID').value
            except KeyError:
                ogo.message('WARNING: Could not find SeriesInstanceUID')
                SeriesInstanceUID = 'n/a'

            try:
                PixelSpacing = FileDataset.data_element('PixelSpacing').value
            except KeyError:
                ogo.message('WARNING: Could not find PixelSpacing')
                PixelSpacing = 0

            sn = SeriesNumber
            if sn not in series_dict.keys():
                series_dict[sn] = {
                    'SeriesDescription':    SeriesDescription,
                    'StudyID':              StudyID,
                    'PatientBirthDate':     PatientBirthDate,
                    'PatientID':            PatientID,
                    'PatientName':          PatientName,
                    'SeriesNumber':         SeriesNumber,
                    'InstanceNumber':       InstanceNumber,
                    'AccessionNumber':      AccessionNumber,
                    'Rows':                 Rows,
                    'Columns':              Columns,
                    'Slices':               1,
                    'PixelSpacing':         PixelSpacing,
                    'SpacingBetweenSlices': SpacingBetweenSlices,
                    'SliceThickness':       SliceThickness,
                    'StudyInstanceUID':     StudyInstanceUID,
                    'SeriesInstanceUID':    SeriesInstanceUID
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
        cmds_remove = ''
        report = ''
        report += '{:30s} {}\n'.format('Directory',dicom_directory)
        report += '{:30s} {}\n'.format('AccessionNumber',v['AccessionNumber'])
        report += '{:30s} {}\n'.format('PatientBirthDate',v['PatientBirthDate'])
        report += '{:30s} {}\n'.format('PatientID',v['PatientID'])
        report += '{:30s} {}\n'.format('PatientName',v['PatientName'])
        report += '{:30s} {}\n'.format('InstanceNumber',v['InstanceNumber'])
        report += '{:30s} {}\n'.format('StudyInstanceUID',v['StudyInstanceUID'])
        report += '\n'
        report += '{:12s} {:40s} {:20s} {:10s} {:15s} {:20s} {:15s} {:20s}\n'.format('SeriesNumber',\
                                                                       'Description',\
                                                                       'Dimensions',\
                                                                       'Slices',\
                                                                       'PixelSpacing',\
                                                                       'SpacingBetweenSlices',\
                                                                       'SliceThickness',\
                                                                       'Filename')
        
        for k,v in series_dict.items():
            
            dims = '{} x {} x {}'.format(v['Rows'],v['Columns'],v['Slices'])
            pixelspacing = '{:.2f} x {:.2f}'.format(v['PixelSpacing'][0],v['PixelSpacing'][1])
            
            if v['SpacingBetweenSlices'] == 0 or v['SpacingBetweenSlices'] == None:
                spacingbetweenslices = '-'
            else:
                spacingbetweenslices = '{:.2f}'.format(v['SpacingBetweenSlices'])
            
            if v['SliceThickness'] == 0 or v['SliceThickness'] == None:
                slicethickness = '-'
            else:
                slicethickness = '{:.2f}'.format(v['SliceThickness'])
            
            # Filename
            series_filename = str(v['SeriesNumber']) + \
                              '_' + \
                              v['SeriesDescription'].lower().replace(' ','_').replace('.','').replace('/','').replace(',','') + \
                              '.nii.gz'
            
            # Highlight
            if ('axial' in v['SeriesDescription'].lower() or 'abdomen' in v['SeriesDescription'].lower() or 'pel' in v['SeriesDescription'].lower() or 'torso' in v['SeriesDescription'].lower()) and \
                v['Slices'] > 100 and \
                v['SpacingBetweenSlices'] < 3.0:
                highlight = '*'
                cmds += 'ls ' + output_folder + '/' + series_filename + '\n'
            else:
                highlight = ' '
        
            report += '{:<7d}{:5s} {:40s} {:20s} {:<10d} {:<15s} {:<20s} {:<15s} {}\n'.format(v['SeriesNumber'],\
                                                                  highlight,\
                                                                  v['SeriesDescription'],\
                                                                  dims,\
                                                                  v['Slices'],\
                                                                  pixelspacing,\
                                                                  spacingbetweenslices,\
                                                                  slicethickness,\
                                                                  series_filename)
        
        report += '\n'
        report += 'Candidate image series: criteria: (\'axial\' or \'pel\' or \'abdomen\' or \'torso\') and (>100 slices) and (slice spacing <3.0mm)\n'
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
        settings.disable_validate_slice_increment()
        settings.disable_validate_slicecount()

        #settings.disable_validate_slice_increment()
        #settings.enable_resampling()
        #settings.set_resample_spline_interpolation_order(1)
        #settings.set_resample_padding(-1000)

        dicom2nifti.convert_directory(dicom_directory, output_folder, compression=True, reorient=True)
    else:
        ogo.message('No NIFTI generated.')
    
    # Sometimes not all image series can be extracted.
    # Try installing the code and running:
    #   $ pip install dicom2nifti
    #   $ dicom2nifti --allow-inconsistent-slice-increment input_directory output_directory
    #
    # Here are valid settings to use in python:
    #    print(dir(dicom2nifti.settings))
    # disable_pydicom_read_force
    # disable_resampling
    # disable_validate_instance_number
    # disable_validate_multiframe_implicit
    # disable_validate_orientation
    # disable_validate_orthogonal
    # disable_validate_slice_increment
    # disable_validate_slicecount
    # enable_pydicom_read_force
    # enable_resampling
    # enable_validate_instance_number
    # enable_validate_multiframe_implicit
    # enable_validate_orientation
    # enable_validate_orthogonal
    # enable_validate_slice_increment
    # enable_validate_slicecount
    # pydicom_read_force
    # resample
    # resample_padding
    # resample_spline_interpolation_order
    # set_resample_padding
    # set_resample_spline_interpolation_order
    # validate_instance_number
    # validate_multiframe_implicit
    # validate_orientation
    # validate_orthogonal
    # validate_slice_increment
    # validate_slicecount

    
    
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

To view DICOM standard data elements:
https://dicom.innolitics.com/ciods/cr-image/patient/00100010

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
    parser.add_argument('--skip_report', action='store_true', help='Skip writing the report (default: %(default)s)')

    print()

    args = parser.parse_args()
    print(echo_arguments('dcm2nii', vars(args)))

    # Run program
    dcm2nii(**vars(args))


if __name__ == '__main__':
    main()
