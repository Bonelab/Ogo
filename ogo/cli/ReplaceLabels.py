# /------------------------------------------------------------------------------+
# | 19-JUL-2022                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+

# Imports
import argparse
import os
import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy
from ogo.util.echo_arguments import echo_arguments
import ogo.cli.Helper as ogo

def histogram(image):
    array = vtk_to_numpy(image.GetPointData().GetScalars()).ravel()
    guard = '!-------------------------------------------------------------------------------'

    if (array.min() < -128):
      range_min = -32768
    elif (array.min() < 0):
      range_min = -128
    else:
      range_min = 0
    
    if (array.max() > 255):
      range_max = 32767
    elif (array.max() > 127):
      range_max = 255
    else:
      range_max = 127

    nRange = [range_min, range_max]
    nBins = 128
    
    # https://numpy.org/doc/stable/reference/generated/numpy.histogram.html
    hist,bin_edges = np.histogram(array,nBins,nRange,None,None,False)
    nValues = sum(hist)

    print(guard)
    print('!>  {:4s} ({:.3s}) : {:s}'.format('Lab','Qty','#Voxels'))
    for bin in range(nBins):
      index = nRange[0] + int(bin * (nRange[1]-nRange[0])/(nBins-1))
      count = hist[bin]/nValues # We normalize so total count = 1
      nStars = int(count*100)
      if (count>0 and nStars==0): # Ensures at least one * if the histogram bin is not zero
        nStars = 1
      if (nStars > 60):
        nStars = 60 # just prevents it from wrapping in the terminal
      if (count>0):
        print('!> {:4d} ({:.3f}): {:d}'.format(index,count,hist[bin]))
    print(guard)

def ReplaceLabels(input_filename, output_filename, inputLabels, outputLabels, overwrite=False):

    # Check if output exists and should overwrite
    if os.path.isfile(output_filename) and not overwrite:
        result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(output_filename))
        if result.lower() not in ['y', 'yes']:
            ogo.message('Not overwriting. Exiting...')
            os.sys.exit()

    # Check that length of input labels equals length of output labels
    if len(inputLabels) is not len(outputLabels):
        ogo.message('ERROR: Number of input labels defined must equal number of output labels.')
        ogo.message('       [input #{:d} != output #{:d}]'.format(len(inputLabels),len(outputLabels)))
        os.sys.exit()
    
    # Check that at least least one set of input/output labels is defined
    if len(inputLabels) < 1:
        ogo.message('ERROR: At least one label must be defined.')
        os.sys.exit()
        
    # Check that labels exist and are in range
    print('Labels to be replaced:')
    for idx,ids in enumerate(inputLabels):
        print('!> {:3d} --> {:3d}'.format(inputLabels[idx],outputLabels[idx]))
        if (inputLabels[idx]<0 or inputLabels[idx]>255 or outputLabels[idx]<0 or outputLabels[idx]>255):
            ogo.message('ERROR: Labels out of range. Must be 0-255.')
            ogo.message('       [Suspects are input label {:d} or output label {:d}]'.format(inputLabels[idx],outputLabels[idx]))
            os.sys.exit()
    
    # Making a map for the labels
    map_for_labels = dict(zip(inputLabels, outputLabels))

    # Read input
    if not os.path.isfile(input_filename):
        os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(input_filename))

    if input_filename.lower().endswith('.nii'):
        reader = vtk.vtkNIFTIImageReader()
    elif input_filename.lower().endswith('.nii.gz'):
        reader = vtk.vtkNIFTIImageReader()
    else:
        os.sys.exit('[ERROR] Cannot find reader for file \"{}\"'.format(input_filename))

    print()
    ogo.message('Reading input image ' + input_filename)
    reader.SetFileName(input_filename)
    reader.Update()
    
    image = reader.GetOutput() # pointer to image
    
    scalarType = image.GetScalarType()
    ogo.message('Input image scalar type: {:s}'.format(image.GetScalarTypeAsString()))
    
    ogo.message('Input image labels:')
    histogram(image)
    
    array = vtk_to_numpy(image.GetPointData().GetScalars())
    
    # Replace each label by cycling through image data; one cycle per label
    for idx,lab in enumerate(map_for_labels):
        ogo.message('!> Replacing {:d} with {:d}'.format(inputLabels[idx],outputLabels[idx]))
        count = np.count_nonzero(array == inputLabels[idx])
        ogo.message('!> --> replaced {:d} labels'.format(count))
        array[array == inputLabels[idx]] = outputLabels[idx]
    
    ogo.message('Output image labels:')
    histogram(image)

    # Create writer
    if output_filename.lower().endswith('.nii'):
        writer = vtk.vtkNIFTIImageWriter()
    elif output_filename.lower().endswith('.nii.gz'):
        writer = vtk.vtkNIFTIImageWriter()
    else:
        os.sys.exit('[ERROR] Cannot find writer for file \"{}\"'.format(output_filename))
          
    ogo.message('Saving output image ' + output_filename)
    writer.SetInputData(image)
    writer.SetFileName(output_filename)
    writer.SetTimeDimension(reader.GetTimeDimension())
    writer.SetTimeSpacing(reader.GetTimeSpacing())
    writer.SetRescaleSlope(reader.GetRescaleSlope())
    writer.SetRescaleIntercept(reader.GetRescaleIntercept())
    writer.SetQFac(reader.GetQFac())
    writer.SetQFormMatrix(reader.GetQFormMatrix())
    writer.SetNIFTIHeader(reader.GetNIFTIHeader())
    writer.Update()

def main():
    # Setup description
    description='''
Utility to read segmented image and replace specified labels with new labels.
It can be used to remove labels from an image (replace with background of 0) or 
to change a label. It works through the list provided sequentially.

Valid input and output file formats include: 
.nii, .nii.gz

Currently only accepts NIFTI file formats as input and output.

Valid input and output labels should be within 0 - 255 (inclusive)
'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoReplaceLabels",
        description=description
    )
    parser.add_argument('input_filename', help='Input image file (*.nii, *.nii.gz)')
    parser.add_argument('output_filename', help='Output image file (*.nii, *.nii.gz)')
    parser.add_argument('--inputLabels', type=int, nargs='*', default=[0], metavar='ID', help='Target label input IDs; space separated (e.g. 1 2 3)')
    parser.add_argument('--outputLabels', type=int, nargs='*', default=[0], metavar='ID', help='Target label output IDs; space separated (e.g. 4 5 6)')
    parser.add_argument('-o', '--overwrite', action='store_true', help='Overwrite output without asking')

    print()

    # Parse and display
    args = parser.parse_args()
    print(echo_arguments('ReplaceLabels', vars(args)))

    # Run program
    ReplaceLabels(**vars(args))

if __name__ == '__main__':
    main()
