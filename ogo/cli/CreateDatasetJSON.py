# /------------------------------------------------------------------------------+
# | 06-JUL-2022                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+

# Imports
import argparse
import os
import numpy as np
from ogo.util.echo_arguments import echo_arguments
from datetime import date
import json

def printFiles(imagesTr_files,labelsTr_files,imagesTs_files):
    guard = '!-------------------------------------------------------------------------------'
    print(guard)
    print('imagesTr:')
    for fname in imagesTr_files:
        print('{:s}'.format(fname))
    print(guard)
    print('labelsTr:')
    for fname in labelsTr_files:
        print('{:s}'.format(fname))
    print(guard)
    print('imagesTs:')
    for fname in imagesTs_files:
        print('{:s}'.format(fname))
    print()
        
def getOnlyNIFTI(file_list):
    for fname in file_list:
        if 'nii' not in fname:
            file_list.remove(fname)
    return file_list
        
def getSortedFileList(dir_name):
    list_of_files = sorted( filter( lambda x: os.path.isfile(os.path.join(dir_name, x)),os.listdir(dir_name) ) )
    return list_of_files

def appendSubdirectory(file_list,subdir):
    idx = 0
    for fname in file_list:
        file_list[idx] = os.path.join(subdir,fname)
        idx += 1
    return file_list
 
#!-------------------------------------------------------------------------------
def CreateDatasetJSON(task_directory, name, description, reference, license, release, modality, tensorImageSize, overwrite=False):

    # Set some variables
    dataset_name = name
    dataset_description = description
    dataset_reference = reference
    dataset_license = license
    dataset_release = release
    dataset_modality = modality
    dataset_tensorImageSize = tensorImageSize

    json_filename = os.path.join(task_directory, "dataset.json")

    # Check if dataset.json file exists and should overwrite
    if os.path.isfile(json_filename) and not overwrite:
        result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(json_filename))
        if result.lower() not in ['y', 'yes']:
            print('Not overwriting. Exiting...')
            os.sys.exit()
    
    # Check for mandatory imagesTr and labelsTr directories
    imagesTr_path = os.path.join(task_directory,'imagesTr/')
    labelsTr_path = os.path.join(task_directory,'labelsTr/')
    if not os.path.exists(imagesTr_path):
        print('ERROR: Expected path to imagesTr. This is a required path.')
        print('       {:s}'.format(imagesTr_path))
        os.sys.exit()
    if not os.path.exists(labelsTr_path):
        print('ERROR: Expected path to labelsTr. This is a required path.')
        print('       {:s}'.format(labelsTr_path))
        os.sys.exit()
    
    # Check for optional imagesTs directory
    imagesTs_path = os.path.join(task_directory,'imagesTs/')
    if not os.path.exists(imagesTs_path):
        print('WARNING: Path to imagesTs does not exist. This is an optional path.')
        print('       {:s}'.format(imagesTs_path))
        print()
    
    # Generate a list of files in each directory
    imagesTr_files = getOnlyNIFTI(getSortedFileList(imagesTr_path))
    labelsTr_files = getOnlyNIFTI(getSortedFileList(labelsTr_path))
    if os.path.exists(imagesTs_path):
        imagesTs_files = getOnlyNIFTI(getSortedFileList(imagesTs_path))
    else:
        imagesTs_files = []
            
    # Check there are files in each mandatory directory
    if len(imagesTr_files)<1:
        print('ERROR: No files found in imagesTr.')
        os.sys.exit()
    if len(labelsTr_files)<1:
        print('ERROR: No files found in labelsTr.')
        os.sys.exit()
        
    # Check that number of imagesTr files is equal number of labelsTr files
    if len(imagesTr_files) != len(labelsTr_files):
        print('ERROR: Number of files in imagesTr must equal labelsTr: {:d} != {:d}.'\
        .format(len(imagesTr_files),len(labelsTr_files)))
        os.sys.exit()
    
    # Check that the files are named the same and in order
    idx = 0
    for fname in imagesTr_files:
        name,ext = os.path.splitext(fname)
        if 'gz' in ext:
            name = os.path.splitext(name)[0] # Manages files with double extension
        if '_0000' in name:
            name = name.replace('_0000','')
        if name not in labelsTr_files[idx]:
            print('ERROR: Expected matching names in imagesTr and labelsTr.')
            print('       File {:s} does not match {:s}'.format(name,labelsTr_files[idx]))
            os.sys.exit()
        idx += 1
    
    # Complete the path for each file by adding subdirectory
    appendSubdirectory(imagesTr_files,'./imagesTr')
    appendSubdirectory(labelsTr_files,'./labelsTr')
    appendSubdirectory(imagesTs_files,'./imagesTs')
    
    # Show the final list file files
    printFiles(imagesTr_files,labelsTr_files,imagesTs_files)

    # Create labels
    labels = {
        0: "Clear Label",
        1: "Femur Right",
        2: "Femur Left",
        3: "Pelvis Right",
        4: "Pelvis Left",
        5: "Sacrum",
        6: "L5",
        7: "L4",
        8: "L3",
        9: "L2",
       10: "L1"
   }

    # Create an array of dict for image/label pairs
    # Help! I cannot figure out this more compact version:
    #     json_dict['training'] = [{'image': "%s" % i, 'label': "%s" % i} for i in imagesTr_files]
    # So instead I do it the long way...
    training = []
    for idx,fname in enumerate(imagesTr_files):
        tmp_dict = {}
        tmp_dict["image"] = imagesTr_files[idx]
        tmp_dict["label"] = labelsTr_files[idx]
        training.append(tmp_dict)
    
    # Create the json file
    json_dict = {}
    json_dict['name'] = dataset_name
    json_dict['description'] = dataset_description
    json_dict['tensorImageSize'] = dataset_tensorImageSize
    json_dict['reference'] = dataset_reference
    json_dict['licence'] = dataset_license
    json_dict['release'] = dataset_release
    json_dict['modality'] = {str(0): dataset_modality}
    json_dict['labels'] = {str(i): labels[i] for i in labels.keys()}
    
    json_dict['numTraining'] = len(imagesTr_files)
    json_dict['numTest'] = len(imagesTs_files)
    json_dict['training'] = training
    json_dict['test'] = ["%s" % i for i in imagesTs_files]

    print('Saving JSON file:\n{:s}'.format(json_filename))
    print(json_dict)
    
    print('Saving JSON file:')
    print(' {:s}'.format(json_filename))
  
    # Dump directly from dictionary
    with open(json_filename, 'w') as outfile:
        json.dump(json_dict, outfile)
    
    #json_string = json.dumps(json_dict)
    #print(json_string)
    
def main():
    # Setup description
    description='''
Data wrangling to prepare for running nnUNet. 

The output of this program is the dataset.json file required for nnUNet.
This defines the directories and files used for training and testing. It
searches the expected directories for the image file names.

Each image file is stored in a Task using the numbers 500 and greater.

nnUNet_raw_data_base/nnUNet_raw_data/
├── Task501_Spine
├── Task502_LFemur
├── Task503_RFemur
├── ...

Task501_Spine/
├── dataset.json
├── imagesTr
├── imagesTs
└── labelsTr

Formatting of file names is crucial:
   imagesTr --> RETRO_005_0000.nii
   labelsTr --> RETRO_005.nii.gz
   imagesTs --> RETRO_070_0000.nii.gz
   
The dataset.json file refers to all three categories of files without the 
trailing "_0000". So, for above files they would be referred to as:
   imagesTr --> RETRO_005.nii
   labelsTr --> RETRO_005.nii.gz
   imagesTs --> RETRO_070.nii.gz

'''

    epilog='''
Prior to running this program, it is expected the appropriate directories, 
file naming, and pre-processing has been completed. Use a script to aid in 
moving raw data into this format.

More details can be found in the publication:
Isensee, F., Jaeger, P.F., Kohl, S.A.A. et al. "nnU-Net: a self-configuring 
method for deep learning-based biomedical image segmentation." Nat Methods 
(2020). https://doi.org/10.1038/s41592-020-01008-z

And see their GIT repository:
https://github.com/MIC-DKFZ/nnUNet
 
'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="CreateDatasetJSON",
        description=description,
        epilog=epilog
    )
    parser.add_argument('task_directory', default=".", metavar='DIR', help='Directory path for Task5XX_...')
    parser.add_argument('--name', metavar='TEXT', default="RETROData", help='Name of dataset (default: %(default)s)')
    parser.add_argument('--description', metavar='TEXT', default="KUB CT abdomen scans including L4 and both femurs", help='Description of dataset (default: %(default)s)')
    parser.add_argument('--reference', metavar='TEXT', default="McCaig Institute, University of Calgary", help='Reference of dataset (default: %(default)s)')
    parser.add_argument('--license', metavar='TEXT', default="CC-BY-SA 4.0", help='License for dataset (default: %(default)s)')
    parser.add_argument('--release', metavar='TEXT', default=date.today().strftime("%B %d, %Y"), help='Date of dataset (default: %(default)s)')
    parser.add_argument('--modality', metavar='TEXT', default="CT", help='Modality of dataset (default: %(default)s)')
    parser.add_argument('--tensorImageSize', metavar='TEXT', default="3D", help='TensorImageSize of dataset (default: %(default)s)')
    parser.add_argument('-o', '--overwrite', action='store_true', help='Overwrite output without asking')

    print()

    args = parser.parse_args()
    print(echo_arguments('CreateDatasetJSON', vars(args)))

    # Run program
    CreateDatasetJSON(**vars(args))

if __name__ == '__main__':
    main()
