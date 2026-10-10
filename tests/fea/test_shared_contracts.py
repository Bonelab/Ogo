"""Contracts between FE builders and shared image helpers."""

import numpy as np
import pytest
import vtk

from ogo.util import Helper
from ogo.util.vtk_image import numpy_to_vtk_image, vtk_image_to_numpy


def image(values):
    template = vtk.vtkImageData()
    template.SetExtent(2, 3, -1, 0, 4, 5)
    template.SetOrigin(10, -20, 30)
    template.SetSpacing(0.7, 1.2, 2.0)
    template.AllocateScalars(vtk.VTK_FLOAT, 1)
    return numpy_to_vtk_image(np.asarray(values, dtype=np.float32), template)


@pytest.mark.parametrize("processing_order", [False, True])
def test_array_roundtrip_preserves_physical_geometry(processing_order):
    source = image(np.arange(8).reshape(2, 2, 2))
    source.SetDirectionMatrix(0, -1, 0, 1, 0, 0, 0, 0, 1)
    values = vtk_image_to_numpy(source, processing_order=processing_order, copy=True)
    output = numpy_to_vtk_image(values, source, processing_order=processing_order)
    np.testing.assert_array_equal(vtk_image_to_numpy(output), vtk_image_to_numpy(source))
    assert output.GetExtent() == source.GetExtent()
    assert output.GetOrigin() == source.GetOrigin()
    assert output.GetSpacing() == source.GetSpacing()
    assert output.GetBounds() == source.GetBounds()
    for i in range(3):
        for j in range(3):
            assert output.GetDirectionMatrix().GetElement(i, j) == source.GetDirectionMatrix().GetElement(i, j)


def test_material_binning_preserves_background_and_geometry():
    source = image([0, 100, 200, 300, 400, 0, 100, 400])
    output, centers = Helper.density2materialID(source, n_bins=2)
    np.testing.assert_array_equal(vtk_image_to_numpy(output).ravel(order="F"),
                                  [0, 1, 1, 2, 2, 0, 1, 2])
    np.testing.assert_allclose(centers, [175, 325])
    assert output.GetBounds() == source.GetBounds()
    assert output.GetExtent() == source.GetExtent()
    with pytest.raises(ValueError, match="empty"):
        Helper.density2materialID(image(np.zeros(8)))


def test_nearest_resampling_keeps_label_ids_and_world_bounds():
    source = image([0, 20, 48, 20, 48, 0, 20, 48])
    transform = vtk.vtkMatrix4x4()
    transform.Identity()
    output = Helper.labelTransformResample(source, transform, 0.5)
    assert output.GetSpacing() == (0.5, 0.5, 0.5)
    assert set(np.unique(vtk_image_to_numpy(output))) <= {0, 20, 48}
    np.testing.assert_allclose(output.GetBounds(), source.GetBounds(), atol=0.5)


def test_bundled_registration_surfaces_are_nonempty():
    from pathlib import Path
    from ogo.fea import femur, spine

    paths = [Path(femur.__file__).parents[1] / "dat" / "LT_FEMUR_SIDEWAYS_FALL_REF.vtk"]
    paths.extend(spine.default_spine_reference_path(target) for target in spine.SPINE_ICP_TARGETS)
    for path in paths:
        surface = Helper.readPolyData(str(path))
        assert surface.GetNumberOfPoints() > 250, path
        assert surface.GetNumberOfCells() > 0, path
        assert np.isfinite(surface.GetBounds()).all(), path


def test_density_resampling_defaults_to_cubic_not_nearest():
    template = vtk.vtkImageData()
    template.SetDimensions(8, 8, 8)
    template.AllocateScalars(vtk.VTK_FLOAT, 1)
    source = numpy_to_vtk_image(np.indices((8, 8, 8)).sum(axis=0).astype(np.float32), template)
    transform = vtk.vtkMatrix4x4()
    transform.Identity()
    default = Helper.transformResample(source, transform, 0.5)
    cubic = Helper.transformResample(source, transform, 0.5, interpolation="cubic")
    nearest = Helper.labelTransformResample(source, transform, 0.5)
    np.testing.assert_array_equal(vtk_image_to_numpy(default), vtk_image_to_numpy(cubic))
    assert not np.array_equal(vtk_image_to_numpy(default), vtk_image_to_numpy(nearest))
