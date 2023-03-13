# /------------------------------------------------------------------------------+
# | 23-OCT-2022                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+

# Imports
import argparse
import os
import vtk
import numpy as np
import SimpleITK as sitk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from ogo.util.echo_arguments import echo_arguments
import ogo.util.Helper as ogo
import ogo.dat.OgoMasterLabels as lb

def getLabels(image):
    
    array = vtk_to_numpy(image.GetPointData().GetScalars()).ravel()
    labels = []
    for k, v in lb.labels_dict.items():
        count = np.count_nonzero(array == k)
        if count>0 and k>0:
            labels.append(k)
    
    return labels
    
def resolve_key_from_value(label_name):
    #labels = [v['LABEL'] for k, v in lb.labels_dict.items() if ((v['TOTSEG']+'.nii.gz') in filename)]
    labels = [v['LABEL'] for k, v in lb.labels_dict.items()]
    if labels:
        keys = [k for k, v in lb.labels_dict.items() if v['LABEL'] == label_name]
        if keys:
            return int(keys[0])
    return None

def boundingbox(image):
    
    dim = image.GetDimensions()
    img = vtk_to_numpy(image.GetPointData().GetScalars()).reshape(dim,order='F')

    
    r = np.any(img, axis=(1, 2))
    c = np.any(img, axis=(0, 2))
    z = np.any(img, axis=(0, 1))
    
    rmin, rmax = np.where(r)[0][[0, -1]]
    cmin, cmax = np.where(c)[0][[0, -1]]
    zmin, zmax = np.where(z)[0][[0, -1]]
    
    # Add some padding if possible
    if zmin>0:
        zmin = zmin-1
    if cmin>0:
        cmin = cmin-1
    if rmin>0:
        rmin = rmin-1
    if zmax<img.shape[2]-1:
        zmax = zmax+1
    if cmax<img.shape[1]-1:
        cmax = cmax+1
    if rmax<img.shape[0]-1:
        rmax = rmax+1
    
    return rmin, rmax, cmin, cmax, zmin, zmax

def Visualize(input_filename, outfile, offscreen, select, collection, gaussian, radius, isosurface, elevation, azimuth, flip, overwrite=False):
    
    # Check if output exists and should overwrite
    if os.path.isfile(outfile) and not overwrite:
        result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(outfile))
        if result.lower() not in ['y', 'yes']:
            ogo.message('Not overwriting. Exiting...')
            os.sys.exit()
    
    if outfile is not 'None' and not (outfile.lower().endswith('.tif')):
        os.sys.exit('[ERROR] Output must be a TIF file ending with .tif: \"{}\"'.format(outfile))
    
    if offscreen and outfile is 'None':
        # Create default output filename for image
        name, ext = os.path.splitext(input_filename)
        if 'gz' in ext:
            name = os.path.splitext(name)[0]  # Manages files with double extension
        outfile = '{}_e{:03d}_a{:03d}.tif'.format(name,elevation,azimuth)
        ogo.message('[WARNING]: No filename specified for output image.')
        ogo.message('           Creating a filename based on input file.')
        ogo.message('           {}'.format(outfile))
    
    # Read input
    if not os.path.isfile(input_filename):
        os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(input_filename))

    if not (input_filename.lower().endswith('.nii') or input_filename.lower().endswith('.nii.gz')):
        os.sys.exit('[ERROR] Input must be type NIFTI file: \"{}\"'.format(input_filename))
    
    ogo.message('Reading input image ' + input_filename)
    reader = vtk.vtkNIFTIImageReader()
    reader.SetFileName(input_filename)
    reader.Update()
    
    image = reader.GetOutput()
    ibounds = image.GetBounds()
    
    ogo.message('Scanning for available labels:')
    labels = getLabels(image)

    for idx,label in enumerate(labels):
        rgb = lb.labels_dict[label]['RGB']
        desc = lb.labels_dict[label]['LABEL']
        ogo.message('{:5d}: ({:3d},{:3d},{:3d}) – {}'.format(label,rgb[0],rgb[1],rgb[2],desc))
    
    # Specify a selection of labels (helpful for TotalSegmentator results)
    if collection:
        if collection == 'all':
            valid_list = [v['LABEL'] for k, v in lb.labels_dict.items()]
        elif collection == 'ossai':
            valid_list = ['Femur Right', 'Femur Left', 'Pelvis Right', 'Pelvis Left', 'Sacrum', 'L6', 'L5', 'L4', 'L3', 'L2', 'L1']
        elif collection == 'skeleton':
            valid_list = ['Femur Right', 'Femur Left', 'Pelvis Right', 'Pelvis Left', 'Sacrum', 'L6', 'L5', 'L4', 'L3', 'L2', 'L1', 'T12', 'T11', 'T10', 'T9', 'T8', 'Humerus Right', 'Humerus Left', 'T7', 'T6', 'T5', 'T4', 'T3', 'T2', 'T1', 'C7', 'C6', 'C5', 'C4', 'C3', 'C2', 'C1', 'Rib Left 1', 'Rib Left 2', 'Rib Left 3', 'Rib Left 4', 'Rib Left 5', 'Rib Left 6', 'Rib Left 7', 'Rib Left 8', 'Rib Left 9', 'Rib Left 10', 'Rib Left 11', 'Rib Left 12', 'Rib Right 1', 'Rib Right 2', 'Rib Right 3', 'Rib Right 4', 'Rib Right 5', 'Rib Right 6', 'Rib Right 7', 'Rib Right 8', 'Rib Right 9', 'Rib Right 10', 'Rib Right 11', 'Rib Right 12', 'Scapula Left', 'Scapula Right', 'Clavicula Left', 'Clavicula Right']
        elif collection == 'cardio':
            valid_list = ['Heart Myocardium', 'Heart Atrium Left', 'Heart Ventricle Left', 'Heart Atrium Right', 'Heart Ventricle Right', 'Aorta', 'Inferior Vena Cava', 'Portal Vein and Splenic Vein', 'Pulmonary Artery', 'Iliac Artery Left', 'Iliac Artery Right', 'Iliac Vena Left', 'Iliac Vena Right']
        elif collection == 'organs':
            valid_list = ['Face', 'Brain', 'Trachea', 'Lung Upper Lobe Left', 'Lung Lower Lobe Left', 'Lung Upper Lobe Right', 'Lung Middle Lobe Right', 'Lung Lower Lobe Right', 'Adrenal Gland Right', 'Adrenal Gland Left', 'Spleen', 'Kidney Right', 'Kidney Left', 'Gallbladder', 'Liver', 'Pancreas']
        elif collection == 'gastro':
            valid_list = ['Esophagus', 'Stomach', 'Duodenum', 'Small Bowel', 'Colon', 'Urinary Bladder']
        elif collection == 'muscle':
            valid_list = ['Autochthon Left', 'Autochthon Right', 'Iliopsoas Left', 'Iliopsoas right', 'Gluteus Maximus Left', 'Gluteus Maximus Right', 'Gluteus Medius Left', 'Gluteus Medius Right', 'Gluteus Minimus Left', 'Gluteus Minimus Right']
        else:
            os.sys.exit('[ERROR] Unknown collection defined: {}'.format(collection))
        
        for item in valid_list:
            select.append(resolve_key_from_value(item))
    
    # Reduce list to user selected labels
    if len(select)>0:
        tmp = []
        for label in select:
            if label in labels:
                tmp.append(label)
            else:
                ogo.message('[WARNING]: Selected label {} [{}] not in input image.'.format(lb.labels_dict[label]['LABEL'],label))
        labels = tmp
    
    # Create lists for VTK classes
    thres = []
    extract = []
    gauss = []
    mcube = []
    actor = []
    mapper = []
    
    # Create renderer
    renderWindow = vtk.vtkRenderWindow()
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(0.0, 0.0, 0.0) # 0.4, 0.5, 0.7
    renderWindow.AddRenderer(renderer)
    renderWindow.SetSize(2000,2000)
    if offscreen:
        renderWindow.SetOffScreenRendering(1)
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow( renderWindow )
    interactor = vtk.vtkInteractorStyleTrackballCamera()
    renderWindowInteractor.SetInteractorStyle(interactor)
    
    for idx,label in enumerate(labels):
        
        ogo.message('Processing label {} ({})'.format(label,lb.labels_dict[label]['LABEL']))
        thres.append(vtk.vtkImageThreshold())
        extract.append(vtk.vtkExtractVOI())
        gauss.append(vtk.vtkImageGaussianSmooth())
        mcube.append(vtk.vtkImageMarchingCubes())
        mapper.append(vtk.vtkDataSetMapper())
        actor.append(vtk.vtkActor())
        rgb = lb.labels_dict[label]['RGB']
        
        thres[idx].SetInputConnection(reader.GetOutputPort())
        thres[idx].SetInValue(127)
        thres[idx].SetOutValue(0)
        thres[idx].ThresholdBetween(label,label)
        thres[idx].Update()
        
        voi = boundingbox(thres[idx].GetOutput())
        extract[idx].SetInputConnection(thres[idx].GetOutputPort())
        extract[idx].SetVOI(voi)
        extract[idx].Update()
        
        gauss[idx].SetInputConnection(extract[idx].GetOutputPort())
        gauss[idx].SetStandardDeviation(gaussian)
        gauss[idx].SetRadiusFactor(radius)
        gauss[idx].Update()
        
        mcube[idx].SetInputConnection(gauss[idx].GetOutputPort())
        mcube[idx].ComputeScalarsOff()
        mcube[idx].ComputeNormalsOff()
        mcube[idx].SetValue(0,isosurface)
        mcube[idx].Update()
        ogo.message('Generated {} triangles.'.format(mcube[idx].GetOutput().GetNumberOfCells()))
        
        mapper[idx].SetInputConnection(mcube[idx].GetOutputPort())
        
        actor[idx].SetMapper(mapper[idx])
        actor[idx].GetProperty().SetColor(rgb[0]/255.0,rgb[1]/255.0,rgb[2]/255.0)
      
        renderer.AddActor( actor[idx] )
    
    # Set up the camera
    if flip:  # Orients camera along Z axis
        renderer.GetActiveCamera().SetViewUp(0,0,1)
    else:
        renderer.GetActiveCamera().SetViewUp(0,0,-1)
        
    renderer.GetActiveCamera().SetPosition((ibounds[1]-ibounds[0]),\
                                           (ibounds[3]-ibounds[2])*5,\
                                           (ibounds[5]-ibounds[4])) # Coordinates of camera 
    renderer.GetActiveCamera().SetFocalPoint((ibounds[1]-ibounds[0]),\
                                             (ibounds[3]-ibounds[2]),\
                                             (ibounds[5]-ibounds[4])) # Focus of camera
    renderer.GetActiveCamera().Elevation(elevation)
    renderer.GetActiveCamera().Azimuth(azimuth)
    renderer.ResetCamera()
    
    # Render the scene
    renderWindow.Render()
    renderWindowInteractor.Initialize()
    if not offscreen:
        renderWindowInteractor.Start()

    # Print image
    if outfile is not 'None':
        windowToImage = vtk.vtkWindowToImageFilter()
        windowToImage.SetInput(renderWindow)
        
        writer = vtk.vtkTIFFWriter()
        writer.SetFileName(outfile)
        writer.SetInputConnection(windowToImage.GetOutputPort())
        writer.Write()
    
    ogo.message('Done.')
    
def main():
    # Setup description
    description='''
Generates an offscreen rendering of a NIFTI file.
'''
    epilog='''
USAGE: 
ogoVisualize bone.nii.gz
ogoVisualize bone.nii.gz --outfile image.tif
ogoVisualize bone.nii.gz --outfile image.tif --overwrite
ogoVisualize bone.nii.gz --outfile image.tif --collection skeleton --overwrite
'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoVisualize",
        description=description,
        epilog=epilog
    )
    parser.add_argument('input_filename', help='Input image file (*.nii, *.nii.gz)')
    parser.add_argument('--select', type=int, nargs='*', default=[], metavar='#',help='Specify labels (e.g. 1 2 3; default: all)')
    parser.add_argument('--collection', default=None, choices=['all', 'ossai', 'skeleton', 'cardio', 'organs', 'muscle', 'gastro'],
                                                           help='Specify collection for visualization (default: %(default)s)')
    parser.add_argument('--gaussian', type=float, default=0.7, metavar='GAUSS',help='Gaussian filter (default: %(default)s)')
    parser.add_argument('--radius', type=int, default=2, metavar='RAD',help='Radius of Gaussian filter (default: %(default)s)')
    parser.add_argument('--isosurface', type=int, default=100, metavar='ISOSURF',help='Isosurface extraction (default: %(default)s)')
    parser.add_argument('--elevation', type=int, default=10, metavar='ELEV',help='Camera elevation (default: %(default)s deg)')
    parser.add_argument('--azimuth', type=int, default=40, metavar='AZI',help='Camera azimuth (default: %(default)s deg)')
    parser.add_argument('--flip', action='store_true', help='Camera ViewUp flip (default: %(default)s)')
    parser.add_argument('-o','--outfile', default='None', metavar='FN', help='Output image file (*.tif) (default: %(default)s)')
    parser.add_argument('--offscreen', action='store_true', help='Set to offscreen rendering (default: %(default)s)')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite output without asking')

    # Parse and display
    args = parser.parse_args()
    print(echo_arguments('Visualize', vars(args)))

    # Run program
    Visualize(**vars(args))

if __name__ == '__main__':
    main()
