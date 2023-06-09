The files in this folder have been deprecated for various reasons. Any of them can still be used, but you must run them
with a command such as:

$ python ReadPickle.py -h

A brief explanation of why they were deprecated is provided here:

1. DicomScanner.py (June 8, 2023) – This was replaced by DicomSelector as a tool to search through a CSV file for a dicom dataset. There is no need for duplication.

2. ReadPickle.py (June 8, 2023) – This function was previously used to look at pkl files, but instead it is more convenient to simply 

3. dcm2nii.py (June 8, 2023) – This function uses dicom2nifti, which is problematic in terms of the proper translation of the SFORM and QFORM transformations from the original DICOM to the NIfTI. It is recommended instead to use this external package: https://github.com/dgobbi/vtk-dicom/tags
