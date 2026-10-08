"""Inspectable hip surfaces and opaque four-view QC from the actual FE mesh.

These are visualization artifacts, not replacement solver models. No mesh
smoothing, scaling, or voxel modification is performed during export.
"""

from pathlib import Path
import os
import warnings

import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy


def export_femur_model(model, model_path, *, pmma_mat_id):
    """Write a colored VTP surface and a four-view PNG beside an n88model.

    ModelRegion cell labels: 0 bone, 1 GT support, 2 FH loading support. The
    supports are identified using their actual BC-node loading-axis positions.
    Distal constraint nodes are shown in green in each QC view.
    """
    prefix = Path(model_path).with_suffix("")
    surface_path = prefix.with_name(prefix.name + "_surface.vtp")
    png_path = prefix.with_name(prefix.name + "_3d.png")
    points = vtk_to_numpy(model.GetPoints().GetData())
    bc = {}
    for key, name in (("gt", "Greater_Trochanter_PMMA_Nodes"),
                      ("fh", "Femoral_Head_PMMA_Nodes"),
                      ("distal", "Distal_Femur_Nodes")):
        ids = model.GetNodeSet(name)
        if ids is None or not ids.GetNumberOfTuples():
            raise ValueError("Cannot export hip QC: missing %s." % name)
        bc[key] = points[vtk_to_numpy(ids).astype(int)]
    grid = vtk.vtkUnstructuredGrid()
    grid.DeepCopy(model)
    centers = vtk.vtkCellCenters()
    centers.SetInputData(grid)
    centers.Update()
    y = vtk_to_numpy(centers.GetOutput().GetPoints().GetData())[:, 1]
    material = grid.GetCellData().GetArray("MaterialID") or grid.GetCellData().GetScalars()
    if material is None:
        raise ValueError("Cannot export hip QC without cell material IDs.")
    pmma = vtk_to_numpy(material) == pmma_mat_id
    region = np.zeros(grid.GetNumberOfCells(), dtype=np.uint8)
    gt_y, fh_y = bc["gt"][:, 1].mean(), bc["fh"][:, 1].mean()
    near_gt = np.abs(y - gt_y) < np.abs(y - fh_y)
    region[pmma & near_gt] = 1
    region[pmma & ~near_gt] = 2
    array = numpy_to_vtk(region, deep=True)
    array.SetName("ModelRegion")
    grid.GetCellData().AddArray(array)
    geometry = vtk.vtkDataSetSurfaceFilter()
    geometry.SetInputData(grid)
    geometry.Update()
    surface = geometry.GetOutput()
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(str(surface_path))
    writer.SetInputData(surface)
    writer.SetDataModeToAppended()
    writer.SetCompressorTypeToZLib()
    if not writer.Write():
        raise RuntimeError("Could not write hip surface: %s" % surface_path)
    outputs = {"surface_vtp": str(surface_path), "qc_3d_png": None,
               "render_status": "unavailable",
               "region_labels": {"0": "bone", "1": "GT support", "2": "FH loading support"}}

    colors = vtk.vtkLookupTable()
    colors.SetNumberOfTableValues(3)
    colors.Build()
    for index, rgb in enumerate(((0.46, 0.65, 0.80), (0.89, 0.64, 0.25), (0.76, 0.24, 0.26))):
        colors.SetTableValue(index, *rgb, 1)
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputData(surface)
    normals.SplittingOff()
    normals.ConsistencyOn()
    bounds = np.asarray(grid.GetBounds()).reshape(3, 2)
    center = bounds.mean(axis=1)
    size = np.ptp(bounds, axis=1)
    radius = float(np.linalg.norm(size))
    window = vtk.vtkRenderWindow()
    if window.GetClassName() == "vtkXOpenGLRenderWindow" and not os.environ.get("DISPLAY"):
        window.Finalize()
        outputs["render_message"] = "3D QC needs an X display or a headless EGL/OSMesa VTK build; surface VTP was saved."
        warnings.warn(outputs["render_message"], RuntimeWarning, stacklevel=2)
        return outputs
    window.SetOffScreenRendering(1)
    window.SetSize(1600, 1000)
    window.SetMultiSamples(4)
    views = (("Lateral", (3, 0, 0), (0, 0.5, 0.5, 1)),
             ("Medial", (-3, 0, 0), (0.5, 0.5, 1, 1)),
             ("Oblique anterior", (2.5, 0.4, -0.8), (0, 0, 0.5, 0.5)),
             ("Oblique posterior", (-2.5, 0.4, 0.8), (0.5, 0, 1, 0.5)))
    try:
        for label, direction, viewport in views:
            renderer = vtk.vtkRenderer()
            renderer.SetBackground(1, 1, 1)
            renderer.SetViewport(*viewport)
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(normals.GetOutputPort())
            mapper.SetScalarModeToUseCellFieldData()
            mapper.SetColorModeToMapScalars()
            mapper.SelectColorArray("ModelRegion")
            mapper.SetLookupTable(colors)
            mapper.SetScalarRange(0, 2)
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor.GetProperty().SetOpacity(1)
            actor.GetProperty().SetAmbient(0.3)
            actor.GetProperty().SetDiffuse(0.7)
            renderer.AddActor(actor)
            cloud = vtk.vtkPolyData()
            cloud_points = vtk.vtkPoints()
            cloud_points.SetData(numpy_to_vtk(bc["distal"], deep=True))
            cloud.SetPoints(cloud_points)
            vertices = vtk.vtkVertexGlyphFilter()
            vertices.SetInputData(cloud)
            node_mapper = vtk.vtkPolyDataMapper()
            node_mapper.SetInputConnection(vertices.GetOutputPort())
            node_actor = vtk.vtkActor()
            node_actor.SetMapper(node_mapper)
            node_actor.GetProperty().SetColor(0.09, 0.39, 0.36)
            node_actor.GetProperty().SetPointSize(3)
            renderer.AddActor(node_actor)
            text = vtk.vtkTextActor()
            text.SetInput(label)
            text.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
            text.SetPosition(0.025, 0.94)
            text.GetTextProperty().SetFontSize(23)
            text.GetTextProperty().SetColor(0.15, 0.15, 0.15)
            renderer.AddActor2D(text)
            if label == "Lateral":
                length_text = vtk.vtkTextActor()
                length = float(bc["gt"][:, 2].min() - np.median(bc["distal"][:, 2]))
                length_text.SetInput("GT disk edge to distal BC: %.1f mm" % length)
                length_text.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
                length_text.SetPosition(0.025, 0.035)
                length_text.GetTextProperty().SetFontSize(22)
                length_text.GetTextProperty().SetColor(0.15, 0.15, 0.15)
                renderer.AddActor2D(length_text)
            camera = renderer.GetActiveCamera()
            camera.ParallelProjectionOn()
            camera.SetFocalPoint(*center)
            camera.SetPosition(*(center + np.asarray(direction) * radius))
            camera.SetViewUp(0, 1, 0)
            window.AddRenderer(renderer)
            # Fit the projected bounds in each viewport, including long shafts.
            renderer.ResetCamera()
            camera.SetParallelScale(camera.GetParallelScale() * 1.1)
            renderer.ResetCameraClippingRange()
        window.Render()
        capture = vtk.vtkWindowToImageFilter()
        capture.SetInput(window)
        capture.Update()
        png = vtk.vtkPNGWriter()
        png.SetFileName(str(png_path))
        png.SetInputConnection(capture.GetOutputPort())
        png.Write()
        if png.GetErrorCode():
            raise RuntimeError("Could not write hip QC: %s" % png_path)
    finally:
        window.Finalize()
    return dict(outputs, qc_3d_png=str(png_path), render_status="rendered")
