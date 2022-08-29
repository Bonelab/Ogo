# MDSC689.03
# Advanced Medical Image Processing
#
# Assignment #2 solution
# Written by Nathan Neeteson
# Created:          2022-Jan-25
# Last Revision:    2022-Jan-25
# --------------------------------------------
#
# Read in a DICOM or NIFTI, apply a median or gaussian filter,
# then render a slice from the image.
#
# Usage instructions:
#
#   python W03.py -h
#
# Example command line to run the script
#
#    python W03.py <path-to-image> -f G -t 100 200
#
# --------------------------------------------
from __future__ import annotations

import vtk
import numpy as np

from mdsc689_03_util import (
    create_reader, get_filename_GUI,
    create_adjust_image_callback,
    create_mappers_mousewheel_callbacks
)

from mdsc689_03_segmentation import double_threshold_vtk, dilate_vtk

from mdsc689_03_parser import create_parser

def get_segmentation(image: vtk.vtkImageData, args: Namespace) -> vtk.vtkImageData:

    # just create a default dilation kernel.. you could make this more
    # configurable if you wanted
    # this shape only looks at adjacent neighbours, not corners
    k = np.zeros((3,3,3))
    k[:,2,2] = 1
    k[2,:,2] = 1
    k[2,2,:] = 1

    if args.filter == 'G':

        # pass through a gaussian filter object
        vtk_filter = vtk.vtkImageGaussianSmooth()
        vtk_filter.SetStandardDeviation(1)
        vtk_filter.SetRadiusFactors(1,1,1)
        vtk_filter.SetDimensionality(3)

        vtk_filter.SetInputData(image)

        vtk_filter.Update()

        segmentation = double_threshold_vtk(vtk_filter.GetOutput(), *args.thresh)

    elif args.filter == 'M':

        vtk_filter = vtk.vtkImageMedian3D()
        vtk_filter.SetKernelSize(3,3,3)

        vtk_filter.SetInputData(image)

        vtk_filter.Update()

        segmentation = double_threshold_vtk(vtk_filter.GetOutput(), *args.thresh)

    else:

        segmentation = double_threshold_vtk(image, *args.thresh)

    for _ in range(args.dilations):

        segmentation = dilate_vtk(segmentation, k)

    return segmentation

def main():

    TEXT_FONT = 20
    TEXT_R, TEXT_G, TEXT_B = 1, 1, 1
    TEXT_POS_X, TEXT_POS_Y = 5, 5

    SEG_WINDOW = 1500
    SEG_LEVEL = 0

    parser = create_parser()
    args = parser.parse_args()

    # if the user did not give a filename on command-line, use a GUI dialog to ask for it
    args.image_filename = args.image_filename if args.image_filename else get_filename_GUI()

    # print the filename to the terminal, so the user could copy and paste it and use
    # it as a command line argument for expedience if they want to re-run the program
    print('Filename selected:')
    print(args.image_filename)


    # load in the image
    reader = create_reader(args.image_filename)

    # create the segmentation
    segmentation = get_segmentation(reader.GetOutput(),args)

    # create all of the onjects for the pipeline

    # image objects
    img_resizer = vtk.vtkImageResize()
    img_mapper = vtk.vtkImageMapper()
    img_actor = vtk.vtkActor2D()

    # segmentation objects
    lookup_table = vtk.vtkLookupTable()
    color_mapper = vtk.vtkImageMapToColors()
    seg_resizer = vtk.vtkImageResize()
    seg_mapper = vtk.vtkImageMapper()
    seg_actor = vtk.vtkActor2D()

    # shared objects
    window = vtk.vtkRenderWindow()
    renderer = vtk.vtkRenderer()
    interactor = vtk.vtkRenderWindowInteractor()
    txt = vtk.vtkTextActor()

    # hook the shared stuff together
    window.AddRenderer(renderer)
    interactor.SetRenderWindow(window)

    # connect the image pipeline
    img_resizer.SetInputConnection(reader.GetOutputPort())
    img_mapper.SetInputConnection(img_resizer.GetOutputPort())
    img_actor.SetMapper(img_mapper)
    renderer.AddActor(img_actor)

    # connect the segmentation pipeline
    color_mapper.SetInputData(segmentation)
    color_mapper.SetLookupTable(lookup_table)
    seg_resizer.SetInputConnection(color_mapper.GetOutputPort())
    seg_mapper.SetInputConnection(seg_resizer.GetOutputPort())
    seg_actor.SetMapper(seg_mapper)
    renderer.AddActor(seg_actor)

    # check what slice index to start at
    starting_slice_index = args.z_slice if args.z_slice else (img_mapper.GetWholeZMax()-img_mapper.GetWholeZMin())//2

    # set image pipeline properties
    img_mapper.SetColorWindow(args.window)
    img_mapper.SetColorLevel(args.level)
    img_mapper.SetZSlice(starting_slice_index)

    # set segmentation pipeline properties
    lookup_table.SetNumberOfTableValues(2)
    lookup_table.SetTableValue(0,0.0,0.0,0.0,0.0)
    lookup_table.SetTableValue(1,0.0,1.0,0.0,1.0) # green
    lookup_table.Build()
    color_mapper.SetOutputFormatToRGBA()
    color_mapper.PassAlphaToOutputOn()
    seg_mapper.SetColorWindow(100)
    seg_mapper.SetColorLevel(0)
    seg_mapper.SetZSlice(starting_slice_index)

    # set text properties
    txt.SetInput(f'Slice: {starting_slice_index}')
    txt.GetTextProperty().SetFontFamilyToArial()
    txt.GetTextProperty().SetFontSize(TEXT_FONT)
    txt.GetTextProperty().SetColor(TEXT_R,TEXT_G,TEXT_B)
    txt.SetDisplayPosition(TEXT_POS_X,TEXT_POS_Y)
    renderer.AddActor(txt)




    # add observers / callbacks / styles
    scroll_forward, scroll_backward = create_mappers_mousewheel_callbacks([img_mapper,seg_mapper],txt)
    interactor.AddObserver('MouseWheelForwardEvent',scroll_forward)
    interactor.AddObserver('MouseWheelBackwardEvent',scroll_backward)

    adjust_image = create_adjust_image_callback(
        reader.GetOutput().GetDimensions(),
        [img_resizer, seg_resizer],
        [img_actor, seg_actor],
        window
    )
    window.AddObserver('ModifiedEvent', adjust_image)

    # start the interactor running
    interactor.Initialize()
    window.Render()
    interactor.Start()

if __name__ == '__main__':
    main()
