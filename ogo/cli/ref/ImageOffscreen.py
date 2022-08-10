# /------------------------------------------------------------------------------+
# | 2-AUG-2022                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+

# Imports
import argparse
import os
import vtk
from ogo.util.echo_arguments import echo_arguments
import ogo.cli.Helper as ogo
import ogo.dat.OgoMasterLabels as lb

def ImageOffscreen(input_filename, output_filename, overwrite=False):

    ogo.message('Hello')
    
    graphics_factory = vtk.vtkGraphicsFactory()
    graphics_factory.SetOffScreenOnlyMode(1)
    graphics_factory.SetUseMesaClasses(1)
    
    #imaging_factory = vtk.vtkImagingFactory()
    #imaging_factory.SetUseMesaClasses(1)
    #
    sphereSource = vtk.vtkSphereSource()
    
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(sphereSource.GetOutputPort())
    
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    
    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.SetOffScreenRendering(1)
    renderWindow.AddRenderer(renderer)
    
    renderer.AddActor(actor)
    renderer.SetBackground(1, 1, 1) # Background color white
    
    renderWindow.Render()

    windowToImageFilter = vtk.vtkWindowToImageFilter()
    windowToImageFilter.SetInput(renderWindow)
    windowToImageFilter.Update()
    
    writer = vtk.vtkPNGWriter()
    writer.SetFileName(output_filename)
    writer.SetInputConnection(windowToImageFilter.GetOutputPort())
    writer.Write()
    
   
def main():
    # Setup description
    description='''
Generates an offscreen rendering of a NIFTI file.
'''
    epilog='''
USAGE: 
ogoImageOffscreen input.nii.gz input.png --overwrite
python ImageOffscreen.py /Users/skboyd/Desktop/ML/kub_mask.nii.gz /Users/skboyd/Desktop/image.png --overwrite
'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoImageOffscreen",
        description=description,
        epilog=epilog
    )
    parser.add_argument('input_filename', help='Input image file (*.nii, *.nii.gz)')
    parser.add_argument('output_filename', help='Output image file (*.png)')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite output without asking')

    print()

    # Parse and display
    args = parser.parse_args()
    print(echo_arguments('ImageOffscreen', vars(args)))

    # Run program
    ImageOffscreen(**vars(args))

if __name__ == '__main__':
    main()
