"""Surface exports preserve model geometry and distinguish actual supports."""

import numpy as np
import pytest
import vtk
import vtkbone
from vtk.util.numpy_support import numpy_to_vtkIdTypeArray, vtk_to_numpy

from tests.fea.test_gt_disk_crop import vtk_image


@pytest.fixture
def hip_model():
    import ogo.util.Helper as helper

    data = np.zeros((4, 8, 12), dtype=np.uint16)
    data[1:3, 2:6, 1:11] = 1
    data[1:3, 1:2, 7:10] = 5000
    data[1:3, 6:7, 8:11] = 5000
    grid = helper.Image2Mesh(vtk_image(data))
    model = vtkbone.vtkboneFiniteElementModel()
    model.ShallowCopy(grid)
    xyz = vtk_to_numpy(model.GetPoints().GetData())
    for name, active in (("Femoral_Head_PMMA_Nodes", xyz[:, 1] == xyz[:, 1].max()),
                         ("Greater_Trochanter_PMMA_Nodes", xyz[:, 1] == xyz[:, 1].min()),
                         ("Distal_Femur_Nodes", xyz[:, 2] == xyz[:, 2].min())):
        ids = numpy_to_vtkIdTypeArray(np.flatnonzero(active), deep=True)
        ids.SetName(name)
        model.AddNodeSet(ids)
    return model


def test_hip_exports_are_colored_nonblank_and_in_physical_coordinates(tmp_path, hip_model):
    from ogo.fea.model_export import export_femur_model

    model = hip_model
    outputs = export_femur_model(model, tmp_path / "test.n88model", pmma_mat_id=5000)
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(outputs["surface_vtp"])
    reader.Update()
    surface = reader.GetOutput()
    assert surface.GetNumberOfPolys() > 0
    assert set(vtk_to_numpy(surface.GetCellData().GetArray("ModelRegion"))) == {0, 1, 2}
    assert surface.GetBounds() == pytest.approx(model.GetBounds())
    png = vtk.vtkPNGReader()
    png.SetFileName(outputs["qc_3d_png"])
    png.Update()
    pixels = vtk_to_numpy(png.GetOutput().GetPointData().GetScalars())
    assert (pixels.min(axis=1) < 200).sum() > 1000
    assert ((pixels[:, 2] > pixels[:, 0].astype(float) * 1.25)
            & (pixels[:, 2] > pixels[:, 1])).sum() > 1000
    assert (pixels[:, 0].astype(float) > pixels[:, 1].astype(float) * 1.5).sum() > 100
    assert png.GetOutput().GetDimensions()[:2] == (1600, 1000)


def test_missing_display_preserves_surface_and_records_missing_png(tmp_path, hip_model, monkeypatch):
    from pathlib import Path
    from ogo.fea.model_export import export_femur_model

    class NoDisplayWindow:
        def GetClassName(self):
            return "vtkXOpenGLRenderWindow"

        def Finalize(self):
            pass

    monkeypatch.delenv("DISPLAY", raising=False)
    monkeypatch.setattr(vtk, "vtkRenderWindow", NoDisplayWindow)
    with pytest.warns(RuntimeWarning, match="headless"):
        outputs = export_femur_model(hip_model, tmp_path / "test.n88model", pmma_mat_id=5000)
    assert Path(outputs["surface_vtp"]).is_file()
    assert outputs["qc_3d_png"] is None
    assert outputs["render_status"] == "unavailable"
