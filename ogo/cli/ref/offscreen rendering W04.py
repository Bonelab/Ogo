# MDSC689.03
# Advanced Medical Image Processing
#
# Assignment #4 solution
# Written by Nathan Neeteson
# Created:          2022-Feb-08
# Last Revision:    2022-Feb-08
# --------------------------------------------
#
# Read in a DICOM or NIFTI, then create a 3D surface and volume
# visualization along with a plane widget
#
# Usage instructions:
#
#   python W04.py -h
#
# Example command line to run the script
#
#    python W04.py <path-to-image>
#
# --------------------------------------------
from __future__ import annotations

import vtk

from mdsc689_03_util import (
    create_reader, get_filename_GUI
)

from mdsc689_03_parser import create_parser


def main():

    # constants / defaults

    BACKGROUND_COLOR_NAME = "AliceBlue"
    OUTLINE_COLOR_NAME = "Black"
    SURFACE_COLOR_NAME = "White"
    VOLUME_COLOR_NAME = "Red"

    GAUSS_STD = 1.5
    GAUSS_RADIUS_FACTOR = 1.5
    GAUSS_DIMENSIONALITY = 3

    MEDIAN_KERNEL_SIZE = 3

    PADDER_WIDTH = 3
    PADDER_VALUE = -500

    PLANE_WIDGET_LEVEL = 300
    PLANE_WIDGET_WINDOW = 1000
    PLANE_WIDGET_SLICE_INDEX = 0
    PLANE_WIDGET_OPACITY = 0.9

    SURFACE_THRESHOLD = 200

    VOLUME_THRESHOLD_LOWER = -50
    VOLUME_THRESHOLD_UPPER = 100
    VOLUME_OPACITY = 0.03

    # parsing

    parser = create_parser()
    args = parser.parse_args()

    # if the user did not give a filename on command-line, use a GUI dialog to ask for it
    args.image_filename = args.image_filename if args.image_filename else get_filename_GUI()

    # print the filename to the terminal, so the user could copy and paste it and use
    # it as a command line argument for expedience if they want to re-run the program
    print('Filename selected:')
    print(args.image_filename)

    # create a colors object

    colors = vtk.vtkNamedColors()

    # create the shared start of the pipeline:
    # reader > (filter) > padder

    reader = create_reader(args.image_filename)
    padder = vtk.vtkImageConstantPad()

    if args.filter == 'G':
        filter = vtk.vtkImageGaussianSmooth()
        filter.SetStandardDeviation(GAUSS_STD)
        filter.SetRadiusFactors(
            GAUSS_RADIUS_FACTOR,
            GAUSS_RADIUS_FACTOR,
            GAUSS_RADIUS_FACTOR
        )
        filter.SetDimensionality(GAUSS_DIMENSIONALITY)
        filter.SetInputConnection(reader.GetOutputPort())

        padder.SetInputConnection(filter.GetOutputPort())

    elif args.filter == 'M':
        filter = vtk.vtkImageMedian3D()
        filter.SetKernelSize(
            MEDIAN_KERNEL_SIZE,
            MEDIAN_KERNEL_SIZE,
            MEDIAN_KERNEL_SIZE
        )
        filter.SetInputConnection(reader.GetOutputPort())

        padder.SetInputConnection(filter.GetOutputPort())

    else:
        padder.SetInputConnection(reader.GetOutputPort())

    # to make the padder do anything, we have to set the output extent of the
    # padder to be one larger than the image on all sides

    padder.Update()
    image_extent = list(padder.GetOutput().GetExtent())

    for i in range(6):
        if i%2:
            image_extent[i] += PADDER_WIDTH
        else:
            image_extent[i] -= PADDER_WIDTH

    padder.SetOutputWholeExtent(image_extent)
    padder.SetConstant(PADDER_VALUE)
    padder.Update()

    # create the shared rendering components of the pipeline:
    # renderer > window > interactor

    window = vtk.vtkRenderWindow()
    renderer = vtk.vtkRenderer()
    interactor = vtk.vtkRenderWindowInteractor()

    window.AddRenderer(renderer)
    interactor.SetRenderWindow(window)

    interactor.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())

    renderer.SetBackground(
        colors.GetColor3d(BACKGROUND_COLOR_NAME)[0],
        colors.GetColor3d(BACKGROUND_COLOR_NAME)[1],
        colors.GetColor3d(BACKGROUND_COLOR_NAME)[2]
    )

    # add a plane widget for each dimension

    plane_widgets = []

    for i in range(3):
        plane_widgets.append(vtk.vtkImagePlaneWidget())
        plane_widgets[-1].SetInputConnection(padder.GetOutputPort())
        plane_widgets[-1].SetPlaneOrientation(i)
        plane_widgets[-1].SetWindowLevel(PLANE_WIDGET_WINDOW,PLANE_WIDGET_LEVEL)
        plane_widgets[-1].SetSliceIndex(PLANE_WIDGET_SLICE_INDEX)
        plane_widgets[-1].GetColorMap().SetOutputFormatToRGB()
        plane_widgets[-1].GetColorMap().PassAlphaToOutputOff()
        plane_widgets[-1].GetTexturePlaneProperty().SetOpacity(PLANE_WIDGET_OPACITY)
        plane_widgets[-1].SetInteractor(interactor)

    # add a surface rendering pipeline
    # marching cubes > poly data mapper > actor

    marching_cubes = vtk.vtkImageMarchingCubes()  # input is vtkDataObject, output is vtkPolyData
    surface_decimator = vtk.vtkDecimatePro()
    surface_mapper = vtk.vtkPolyDataMapper()
    surface_actor = vtk.vtkActor()

    marching_cubes.ComputeGradientsOn()
    marching_cubes.ComputeNormalsOn()
    marching_cubes.ComputeScalarsOff()
    marching_cubes.SetNumberOfContours(1)
    marching_cubes.SetValue(0,SURFACE_THRESHOLD)
    
    surface_decimator.PreserveTopologyOn()
    surface_decimator.SetTargetReduction(0.5)

    surface_actor.GetProperty().SetColor(
        colors.GetColor3d(SURFACE_COLOR_NAME)[0],
        colors.GetColor3d(SURFACE_COLOR_NAME)[1],
        colors.GetColor3d(SURFACE_COLOR_NAME)[2]
    )

    marching_cubes.SetInputConnection(padder.GetOutputPort())
    #surface_mapper.SetInputConnection(marching_cubes.GetOutputPort())
    surface_decimator.SetInputConnection(marching_cubes.GetOutputPort())
    surface_mapper.SetInputConnection(surface_decimator.GetOutputPort())
    surface_actor.SetMapper(surface_mapper)

    renderer.AddActor(surface_actor)

    # add a volume rendering pipeline
    # opacity function, color function > volume property > ray cast mapper > volume

    opacity_function = vtk.vtkPiecewiseFunction()
    color_function = vtk.vtkColorTransferFunction()
    volume_property = vtk.vtkVolumeProperty()
    volume_mapper = vtk.vtkOpenGLGPUVolumeRayCastMapper()
    volume = vtk.vtkVolume()

    opacity_function.AddPoint(VOLUME_THRESHOLD_LOWER, VOLUME_OPACITY)
    opacity_function.AddPoint(VOLUME_THRESHOLD_UPPER, VOLUME_OPACITY)
    opacity_function.ClampingOff()

    color_function.AddRGBPoint(
        VOLUME_THRESHOLD_LOWER,
        colors.GetColor3d(VOLUME_COLOR_NAME)[0],
        colors.GetColor3d(VOLUME_COLOR_NAME)[1],
        colors.GetColor3d(VOLUME_COLOR_NAME)[2]
    )

    color_function.AddRGBPoint(
        VOLUME_THRESHOLD_UPPER,
        colors.GetColor3d(VOLUME_COLOR_NAME)[0],
        colors.GetColor3d(VOLUME_COLOR_NAME)[1],
        colors.GetColor3d(VOLUME_COLOR_NAME)[2]
    )

    volume_property.ShadeOn()
    volume_property.SetInterpolationTypeToLinear()
    volume_mapper.SetBlendModeToComposite()

    volume_mapper.SetInputConnection(padder.GetOutputPort())
    volume_property.SetColor(color_function)
    volume_property.SetScalarOpacity(opacity_function)
    volume.SetProperty(volume_property)
    volume.SetMapper(volume_mapper)

    renderer.AddVolume(volume)

    # add an outline

    outline = vtk.vtkOutlineFilter()
    outline_mapper = vtk.vtkPolyDataMapper()
    outline_actor = vtk.vtkActor()

    outline.SetInputConnection(padder.GetOutputPort())
    outline_mapper.SetInputConnection(outline.GetOutputPort())
    outline_actor.SetMapper(outline_mapper)
    outline_actor.GetProperty().SetColor(
        colors.GetColor3d(OUTLINE_COLOR_NAME)[0],
        colors.GetColor3d(OUTLINE_COLOR_NAME)[1],
        colors.GetColor3d(OUTLINE_COLOR_NAME)[2]
    )

    renderer.AddActor(outline_actor)

    # start the interactor running

    renderer.ResetCamera()
    interactor.Initialize()
    window.Render()
    for plane_widget in plane_widgets:
        plane_widget.On()
    interactor.Start()


if __name__ == '__main__':
    main()
