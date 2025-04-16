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

def getNIfTITransform(reader):
    # Grab the NIfTI transform (see NIH https://nifti.nimh.nih.gov/pub/dist/src/niftilib/nifti1.h)
    found=False
    if reader.GetQFormMatrix():
        mat = reader.GetQFormMatrix()
        ogo.message('QForm transform:')
        found=True
    elif reader.GetSFormMatrix():
        mat = reader.GetSFormMatrix()
        if found:
            ogo.message('[WARNING]: Both QForm and SForm have been found. Using SForm.')
        ogo.message('SForm transform:')
    else:
        mat = vtk.vtkMatrix4x4()
        ogo.message('[WARNING]: Neither QForm nor SForm transform found. Using Identity.')
        ogo.message('Identity transform:')

    for i in range(4):
        ogo.message('['+' '.join('{:8.3f}'.format(mat.GetElement(i,j)) for j in range(4))+']')
    
    return mat
    
def get_lookup_table():
    lut = vtk.vtkLookupTable()
    lut.SetNumberOfTableValues(256)
    lut.Build()
    #lut.SetRange(0,255)
    for idx in range(255):
        if lb.labels_dict.get(idx):
            vis = lb.labels_dict[idx]['VIS']
            rgb = lb.labels_dict[idx]['RGB']
            red = rgb[0]/255
            gre = rgb[1]/255
            blu = rgb[2]/255
            #print('{:3d} {:6.3f} {:6.3f} {:6.3f} {:6.1f}'.format(idx, red, gre, blu, vis))
            lut.SetTableValue(idx, red, gre, blu, vis)
        else:
            #print(idx)
            lut.SetTableValue(idx, 0, 0, 0, 0)
    return lut
    
def vis3d(input_filename, select, collection, gaussian, radius, isosurface, elevation, azimuth, flip, outfile, offscreen, pad, overwrite, func):
    
    # Check if output exists and should overwrite
    if outfile:
        if os.path.isfile(outfile) and not overwrite:
            result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(outfile))
            if result.lower() not in ['y', 'yes']:
                ogo.message('Not overwriting. Exiting...')
                os.sys.exit()
    
    if (outfile is not None) and (not (outfile.lower().endswith('.tif'))):
        os.sys.exit('[ERROR] Output must be a TIF file ending with .tif: \"{}\"'.format(outfile))
    
    if offscreen and outfile is None:
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
    
    ibounds = reader.GetOutput().GetBounds()
    mat = getNIfTITransform(reader)

    extent = reader.GetOutput().GetExtent()
    ogo.message('  Extent: ' + ', '.join('{:d}'.format(i) for i in extent))
    ogo.message(f'Offscreen rendering: {offscreen}')
 
    ogo.message('Scanning for available labels:')
    labels = getLabels(reader.GetOutput())

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
    padding = []
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
        padding.append(vtk.vtkImageConstantPad())
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
        
        padding[idx].SetInputConnection(extract[idx].GetOutputPort())
        padding[idx].SetConstant(0)
        padding[idx].SetOutputWholeExtent(voi[0]-1,voi[1]+1,\
                                          voi[2]-1,voi[3]+1,\
                                          voi[4]-1,voi[5]+1)
        padding[idx].Update()
        #ogo.message('  Extent: ' + ', '.join('{:d}'.format(i) for i in padding[idx].GetOutput().GetExtent()))
        
        if pad:
            gauss[idx].SetInputConnection(padding[idx].GetOutputPort())
        else:
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
        actor[idx].SetUserMatrix(mat)
        actor[idx].GetProperty().SetColor(rgb[0]/255.0,rgb[1]/255.0,rgb[2]/255.0)
      
        renderer.AddActor( actor[idx] )
    
    # Set up the camera
    if flip:  # Orients camera along Z axis
        renderer.GetActiveCamera().SetViewUp(0,0,-1)
    else:
        renderer.GetActiveCamera().SetViewUp(0,0,1)
        
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
    if outfile is not None:
        windowToImage = vtk.vtkWindowToImageFilter()
        windowToImage.SetInput(renderWindow)
        
        writer = vtk.vtkTIFFWriter()
        writer.SetFileName(outfile)
        writer.SetInputConnection(windowToImage.GetOutputPort())
        writer.Write()
    
    ogo.message('Done.')
    
def vis2d(input_filename, outfile, mask_image, window, level, nThreads, image_orientation, slice_percent, offscreen, overwrite, func):
    
    # Check if output exists and should overwrite
    if outfile:
        if os.path.isfile(outfile) and not overwrite:
            result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(outfile))
            if result.lower() not in ['y', 'yes']:
                ogo.message('Not overwriting. Exiting...')
                os.sys.exit()
    
    if (outfile is not None) and (not (outfile.lower().endswith('.tif'))):
        os.sys.exit('[ERROR] Output must be a TIF file ending with .tif: \"{}\"'.format(outfile))
    
    if offscreen and outfile is None:
        # Create default output filename for image
        name, ext = os.path.splitext(input_filename)
        if 'gz' in ext:
            name = os.path.splitext(name)[0]  # Manages files with double extension
        outfile = '{}_{:s}_{:.0f}_2d.tif'.format(name,image_orientation,slice_percent)
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
    
    mat = getNIfTITransform(reader)
    
    image = reader.GetOutput()
    #image_bounds = image.GetBounds()
    
    # Get scalar range for W/L and padding
    scalarRanges = reader.GetOutput().GetScalarRange()

    # Determine if we need to autocompute the window/level
    if window <= 0:
        window = scalarRanges[1] - scalarRanges[0]
        level = (scalarRanges[1] + scalarRanges[0])/2

    # Setup input Mapper + Property -> Slice
    inputMapper = vtk.vtkOpenGLImageSliceMapper()
    inputMapper.SetInputConnection(reader.GetOutputPort())
    inputMapper.SliceAtFocalPointOn()
    inputMapper.SliceFacesCameraOn()
    inputMapper.BorderOn()
    inputMapper.SetNumberOfThreads(nThreads)
    inputMapper.StreamingOn()
        
    imageProperty = vtk.vtkImageProperty()
    imageProperty.SetColorLevel(level)
    imageProperty.SetColorWindow(window)
    imageProperty.SetInterpolationTypeToNearest()

    inputSlice = vtk.vtkImageSlice()
    inputSlice.SetMapper(inputMapper)
    inputSlice.SetProperty(imageProperty)
    inputSlice.SetUserMatrix(mat)
    
    image_bounds = inputSlice.GetBounds() # We get the bounds *after* applying the transform

    # Create Renderer -> RenderWindow -> RenderWindowInteractor -> InteractorStyle
    renderer = vtk.vtkRenderer()
    renderer.AddActor(inputSlice)

    # Read mask
    if not mask_image is None:
        if not os.path.isfile(mask_image):
            os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(mask_image))
        
        if not (mask_image.lower().endswith('.nii') or mask_image.lower().endswith('.nii.gz')):
            os.sys.exit('[ERROR] Input must be type NIFTI file: \"{}\"'.format(mask_image))
        
        ogo.message('Reading input mask image ' + mask_image)
        mask_reader = vtk.vtkNIFTIImageReader()
        mask_reader.SetFileName(mask_image)
        mask_reader.Update()
        
        lut = get_lookup_table()
        
        mask_inputMapper = vtk.vtkOpenGLImageSliceMapper()
        mask_inputMapper.SetInputConnection(mask_reader.GetOutputPort())
        mask_inputMapper.SliceAtFocalPointOn()
        mask_inputMapper.SliceFacesCameraOn()
        mask_inputMapper.BorderOn()
        mask_inputMapper.SetNumberOfThreads(nThreads)
        mask_inputMapper.StreamingOn()
        
        mask_imageProperty = vtk.vtkImageProperty()
        mask_imageProperty.SetInterpolationTypeToNearest()
        mask_imageProperty.SetLookupTable(lut)
        mask_imageProperty.SetOpacity(0.5)

        mask_inputSlice = vtk.vtkImageSlice()
        mask_inputSlice.SetMapper(mask_inputMapper)
        mask_inputSlice.SetProperty(mask_imageProperty)
        mask_inputSlice.SetUserMatrix(mat)
        
        renderer.AddActor(mask_inputSlice)
        
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.SetSize(2048,2048)
    renderWindow.AddRenderer(renderer)
    if offscreen:
        renderWindow.SetOffScreenRendering(1)

    interactorStyle = vtk.vtkInteractorStyleImage()
    interactorStyle.SetInteractionModeToImageSlicing()
    interactorStyle.SetCurrentRenderer(renderer)
      
    # Set image orientation
    ogo.message('Image orientation set to {}'.format(image_orientation))
    # First vector – direction corresponding to moving horizontally left-to-right across the screen,
    # Second vector – direction corresponding to moving bottom-to-top up the screen 
    if image_orientation == 'sagittal': # 'x' key
        interactorStyle.SetImageOrientation((0,1,0), (0,0,1)) # last vector was 0,0,1
    elif image_orientation == 'coronal': # 'y' key
        interactorStyle.SetImageOrientation((-1,0,0), (0,0,1)) # last vector was 0,0,1
    elif image_orientation == 'axial': # 'z' key
        interactorStyle.SetImageOrientation((1,0,0), (0,1,0))
    else:
        ogo.message('Error selecting image orientation')
        os.sys.exit()
    renderer.ResetCamera()
    
    # Set slice
    ogo.message('Percent slice set to {:.1f}%'.format(slice_percent))
    camera = renderer.GetActiveCamera()
    fp = list(camera.GetFocalPoint())
    #print('** Focal point: '+' '.join('{:8.3f}'.format(i) for i in fp))
    #print('** Bounds:      '+' '.join('{:8.3f}'.format(i) for i in image_bounds))
    
    if image_orientation == 'sagittal': # 'x' key
        fp[0] = image_bounds[0] + ((image_bounds[1] - image_bounds[0]) * slice_percent / 100.0)
    elif image_orientation == 'coronal': # 'y' key
        fp[1] = image_bounds[2] + ((image_bounds[3] - image_bounds[2]) * slice_percent / 100.0)
    elif image_orientation == 'axial': # 'z' key
        fp[2] = image_bounds[4] + ((image_bounds[5] - image_bounds[4]) * slice_percent / 100.0)
    else:
        ogo.message('Error selecting percent slice')
        os.sys.exit()
    camera.SetFocalPoint(fp)
    #print('** Focal point: '+' '.join('{:8.3f}'.format(i) for i in fp))

    fp = list(renderer.GetActiveCamera().GetFocalPoint())
    ogo.message('Focal point is '+', '.join('{:.1f}'.format(i) for i in fp))

    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetInteractorStyle(interactorStyle)
    interactor.SetRenderWindow(renderWindow)

    # Add some functionality to switch layers for window/level
    def layerSwitcher(obj,event):
        if str(interactor.GetKeyCode()) == 'w':
            ogo.message("Image W/L: {w}/{l}".format(w=imageProperty.GetColorWindow(), l=imageProperty.GetColorLevel()))
        elif str(interactor.GetKeyCode()) == 'n':
            # Set interpolation to nearest neighbour (good for voxel visualization)
            imageProperty.SetInterpolationTypeToNearest()
            interactor.Render()
        elif str(interactor.GetKeyCode()) == 'c':
            # Set interpolation to cubic (makes a better visualization)
            imageProperty.SetInterpolationTypeToCubic()
            interactor.Render()
        elif str(interactor.GetKeyCode()) == 'r':
            window = scalarRanges[1] - scalarRanges[0]
            level = (scalarRanges[1] + scalarRanges[0])/2
            imageProperty.SetColorLevel(level)
            imageProperty.SetColorWindow(window)
            interactor.Render()
        elif str(interactor.GetKeyCode()) == 'u':
            x, y = interactor.GetEventPosition()
            picker = vtk.vtkCellPicker()
            picker.Pick(x,y,0,renderer)
            point = picker.GetPointIJK()
            ogo.message('IJK: '+', '.join('{:5d}'.format(i) for i in point))
        
    # Add ability to switch between active layers
    interactor.AddObserver('KeyPressEvent', layerSwitcher, -1.0) # Call layerSwitcher as last observer
    
    # Initialize and go
    interactor.Initialize()
    if not offscreen:
        interactor.Start()

    if outfile is not None:
        interactor.Render()
        windowToImage = vtk.vtkWindowToImageFilter()
        windowToImage.SetInput(renderWindow)
        
        writer = vtk.vtkTIFFWriter()
        writer.SetFileName(outfile)
        writer.SetInputConnection(windowToImage.GetOutputPort())
        writer.Write()
        ogo.message('Writing {}'.format(outfile))
    
    ogo.message('Done.')
    
def main():
    # Setup description
    description='''
Visualize NIfTI files in either 2d or 3d.

For 3d visualization, only segmented datasets are expected.

For 2d visualization, keyboard mappings are:
    w   Print window/Level to terminal
    n   Set interpolator to nearest neighbour
    c   Set interpolator to cubic
    r   Reset window/level
    u   Print current cursor voxel position
    x   View in x-plane (saggital)
    y   View in y-plane (coronal)
    z   View in z-plane (axial)
    q   Quit

For 2d visualization, mouse mappings are:
    left click + vertical scroll                Modify window
    left click + horizontal scroll              Modify level
    right click + vertical scroll               Zoom
    control + left click + vertical scroll      Slice level
    control + right click + vertical scroll     Rotate slice
    shift + left click + vertical scroll        Translate slice

'''
    epilog='''
USAGE: 
ogoVisualize vis3d bone_labels.nii.gz
ogoVisualize vis3d bone_labels.nii.gz --offscreen
ogoVisualize vis3d bone_labels.nii.gz --outfile image.tif --collection skeleton

ogoVisualize vis2d RETRO_00053.nii.gz
'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoVisualize",
        description=description,
        epilog=epilog
    )
    subparsers = parser.add_subparsers()
    
    # 3d visualization
    parser_vis3d = subparsers.add_parser('vis3d')
    parser_vis3d.add_argument('input_filename', help='Input segmented image file (*.nii, *.nii.gz)')
    parser_vis3d.add_argument('--select', type=int, nargs='*', default=[], metavar='#',help='Specify labels (e.g. 1 2 3; default: all)')
    parser_vis3d.add_argument('--collection', default=None, choices=['all', 'ossai', 'skeleton', 'cardio', 'organs', 'muscle', 'gastro'],
                                                           help='Specify collection for visualization (default: %(default)s)')
    parser_vis3d.add_argument('--gaussian', type=float, default=0.7, metavar='GAUSS',help='Gaussian filter (default: %(default)s)')
    parser_vis3d.add_argument('--radius', type=int, default=2, metavar='RAD',help='Radius of Gaussian filter (default: %(default)s)')
    parser_vis3d.add_argument('--isosurface', type=int, default=100, metavar='ISOSURF',help='Isosurface extraction (default: %(default)s)')
    parser_vis3d.add_argument('--elevation', type=int, default=10, metavar='ELEV',help='Camera elevation (default: %(default)s deg)')
    parser_vis3d.add_argument('--azimuth', type=int, default=40, metavar='AZI',help='Camera azimuth (default: %(default)s deg)')
    parser_vis3d.add_argument('--flip', action='store_true', help='Camera ViewUp flip (default: %(default)s)')
    parser_vis3d.add_argument('--outfile', default=None, metavar='FN', help='Output image file (*.tif) (default: %(default)s)')
    parser_vis3d.add_argument('--offscreen', action='store_true', help='Set to offscreen rendering (default: %(default)s)')
    parser_vis3d.add_argument('--pad', action='store_true', help='Pad the image with zeros')
    parser_vis3d.add_argument('--overwrite', action='store_true', help='Overwrite output without asking')
    parser_vis3d.set_defaults(func=vis3d)

    # 2d visualization
    parser_vis2d = subparsers.add_parser('vis2d')
    parser_vis2d.add_argument('input_filename', help='Input image file (*.nii, *.nii.gz)')
    parser_vis2d.add_argument('--mask_image', default=None, metavar='NIfTI', 
                                        help='Superimpose a mask (*.nii.gz, default: %(default)s)')
    parser_vis2d.add_argument('--window', type=float, default=0.0, metavar='WINDOW',help='Initial window. For <=0, computed from dynamic range (default: %(default)s)')
    parser_vis2d.add_argument('--level', type=float, default=0.0, metavar='LEVEL',help='Initial level. For <=0, computed from dynamic range (default: %(default)s)')
    parser_vis2d.add_argument('--nThreads', type=int, default=1, metavar='THREAD',help='Number of threads (default: %(default)s)')
    parser_vis2d.add_argument('--image_orientation', default='coronal', choices=['sagittal', 'coronal', 'axial'],
                                help='Initiate a particular orientation (default: %(default)s)')
    parser_vis2d.add_argument('--slice_percent', default=float(50), type=float, help='Set percent slice through image (default: %(default)s)')
    parser_vis2d.add_argument('--outfile', default=None, metavar='FN', help='Output image file (*.tif) (default: %(default)s)')
    parser_vis2d.add_argument('--offscreen', action='store_true', help='Set to offscreen rendering (default: %(default)s)')
    parser_vis2d.add_argument('--overwrite', action='store_true', help='Overwrite output without asking')
    parser_vis2d.set_defaults(func=vis2d)

    # Parse and display
    args = parser.parse_args()
    print(echo_arguments('Visualize', vars(args)))

    # Run program
    args.func(**vars(args))

if __name__ == '__main__':
    main()
