import unittest

import pytest

import ogo.util.Helper as helper
from ogo.util.Helper import get_phantom


class Test_Get_Phantom(unittest.TestCase):

    def setUp(self) -> None:
        """ Use this method to prep data or files for other tests, if necessary. """
        pass

    def tearDown(self) -> None:
        """ Use this method to clean up after the tests, if necessary."""
        pass

    @unittest.skip("Placeholder - add more tests here")
    def test_placeholder(self) -> None:
        """ ogo.util.Test_Helper:Test_GetPhantom """
        pass

# There are many more functions in util.Helper that I did not write tests for.
# It's not clear that these functions are actually used, so I didn't bother.


def test_prepare_finite_element_image_preserves_physical_metadata():
    vtk = pytest.importorskip("vtk")
    from vtk.util.numpy_support import numpy_to_vtk

    image = vtk.vtkImageData()
    image.SetDimensions(3, 4, 5)
    image.SetOrigin(11.0, -22.0, 33.0)
    image.SetSpacing(0.7, 0.8, 0.9)
    values = numpy_to_vtk([1] * (3 * 4 * 5), deep=True)
    image.GetPointData().SetScalars(values)

    assert hasattr(helper, "prepareFiniteElementImage")

    prepared = helper.prepareFiniteElementImage(image)

    assert prepared.GetOrigin() == pytest.approx((11.0, -22.0, 33.0))
    assert prepared.GetSpacing() == pytest.approx((0.7, 0.8, 0.9))
    assert prepared.GetDimensions() == (3, 4, 5)


def test_prepare_finite_element_image_mesh_uses_physical_coordinates():
    vtk = pytest.importorskip("vtk")
    pytest.importorskip("vtkbone")
    from vtk.util.numpy_support import numpy_to_vtk

    image = vtk.vtkImageData()
    image.SetDimensions(2, 2, 2)
    image.SetOrigin(11.0, -22.0, 33.0)
    image.SetSpacing(1.0, 1.0, 1.0)
    image.GetPointData().SetScalars(numpy_to_vtk([1] * 8, deep=True))

    mesh = helper.Image2Mesh(helper.prepareFiniteElementImage(image))

    assert mesh.GetBounds() == pytest.approx((10.5, 12.5, -22.5, -20.5, 32.5, 34.5))


def test_sideways_fall_combiner_preserves_existing_material_voxels():
    np = pytest.importorskip("numpy")
    vtk = pytest.importorskip("vtk")
    from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

    pmma_mat_id = 5000
    material = np.zeros((4, 4, 4), dtype=np.int16)
    material[1, 1, 1] = 1
    material[2, 2, 2] = 77

    femoral_head_cap = np.zeros_like(material)
    femoral_head_cap[1, 1, 1] = pmma_mat_id
    femoral_head_cap[0, 1, 1] = pmma_mat_id

    greater_trochanter_cap = np.zeros_like(material)
    greater_trochanter_cap[2, 2, 2] = pmma_mat_id
    greater_trochanter_cap[3, 2, 2] = pmma_mat_id

    def vtk_image(array):
        image = vtk.vtkImageData()
        image.SetDimensions(array.shape)
        image.SetOrigin(0.0, 0.0, 0.0)
        image.SetSpacing(1.0, 1.0, 1.0)
        image.GetPointData().SetScalars(numpy_to_vtk(array.ravel(order="F"), deep=True))
        return image

    combined = helper.combineImageData_SF(
        vtk_image(material),
        vtk_image(femoral_head_cap),
        vtk_image(greater_trochanter_cap),
        pmma_mat_id,
    )

    out = vtk_to_numpy(combined.GetPointData().GetScalars()).reshape(material.shape, order="F")
    assert out[1, 1, 1] == 1
    assert out[2, 2, 2] == 77
    assert out[0, 1, 1] == pmma_mat_id
    assert out[3, 2, 2] == pmma_mat_id


def test_iterative_closest_point_uses_current_fea_defaults(monkeypatch):
    calls = {}

    class FakeLandmarkTransform:
        def SetModeToRigidBody(self):
            calls["rigid_body"] = True

    class FakeIcp:
        def SetSource(self, source):
            calls["source"] = source

        def SetTarget(self, target):
            calls["target"] = target

        def StartByMatchingCentroidsOn(self):
            calls["centroid_start"] = True

        def GetLandmarkTransform(self):
            return FakeLandmarkTransform()

        def SetMeanDistanceModeToAbsoluteValue(self):
            calls["mean_distance_mode"] = "absolute"

        def SetMeanDistanceModeToRMS(self):
            calls["mean_distance_mode"] = "rms"

        def SetMaximumMeanDistance(self, value):
            calls["maximum_mean_distance"] = value

        def CheckMeanDistanceOn(self):
            calls["check_mean_distance"] = True

        def SetMaximumNumberOfLandmarks(self, value):
            calls["maximum_landmarks"] = value

        def SetMaximumNumberOfIterations(self, value):
            calls["maximum_iterations"] = value

        def Update(self):
            calls["updated"] = True

        def GetMatrix(self):
            return "matrix"

    monkeypatch.setattr(helper.vtk, "vtkIterativeClosestPointTransform", FakeIcp)

    matrix = helper.iterativeClosestPoint("source", "target")

    assert matrix == "matrix"
    assert calls["source"] == "source"
    assert calls["target"] == "target"
    assert calls["centroid_start"] is True
    assert calls["rigid_body"] is True
    assert calls["mean_distance_mode"] == "absolute"
    assert calls["maximum_landmarks"] == 8000
    assert calls["maximum_iterations"] == 50


def test_pre_rotate_image_zero_degrees_preserves_image_orientation():
    np = pytest.importorskip("numpy")
    vtk = pytest.importorskip("vtk")
    from ogo.util.vtk_image import numpy_to_vtk_image, vtk_image_to_numpy

    image_array = np.zeros((5, 5, 5), dtype=np.float32)
    mask_array = np.zeros((5, 5, 5), dtype=np.uint8)
    image_array[1, 2, 3] = 99.0
    mask_array[1, 2, 3] = 1

    template = vtk.vtkImageData()
    template.SetDimensions(image_array.shape)
    template.SetOrigin(10.0, -20.0, 30.0)
    template.SetSpacing(1.0, 1.0, 1.0)
    image = numpy_to_vtk_image(image_array, template, vtk_array_type=vtk.VTK_FLOAT)
    mask = numpy_to_vtk_image(mask_array, template, vtk_array_type=vtk.VTK_UNSIGNED_CHAR)

    image_rot, mask_rot = helper.preRotateImage(image, mask, 0)

    assert vtk_image_to_numpy(image_rot) == pytest.approx(image_array)
    assert np.array_equal(vtk_image_to_numpy(mask_rot), mask_array)


def test_read_nii_returns_ras_physical_vtk_image(tmp_path):
    np = pytest.importorskip("numpy")
    sitk = pytest.importorskip("SimpleITK")

    array_zyx = np.arange(2 * 3 * 4, dtype=np.float32).reshape((2, 3, 4))
    image = sitk.GetImageFromArray(array_zyx)
    image.SetSpacing((0.5, 0.75, 1.25))
    image.SetOrigin((10.0, 20.0, 30.0))
    image.SetDirection((1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0))
    path = tmp_path / "image.nii.gz"
    sitk.WriteImage(image, str(path))

    expected = sitk.DICOMOrient(sitk.ReadImage(str(path)), "RAS")
    vtk_image = helper.readNii(str(path))

    assert vtk_image.GetDimensions() == expected.GetSize()
    assert vtk_image.GetSpacing() == pytest.approx(expected.GetSpacing())
    assert vtk_image.GetOrigin() == pytest.approx(
        (-expected.GetOrigin()[0], -expected.GetOrigin()[1], expected.GetOrigin()[2])
    )


def test_read_poly_data_returns_ras_reference_points(tmp_path):
    vtk = pytest.importorskip("vtk")

    points = vtk.vtkPoints()
    points.InsertNextPoint(1.0, 2.0, 3.0)
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    path = tmp_path / "reference.vtk"
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(str(path))
    writer.SetInputData(polydata)
    writer.Update()

    reference = helper.readPolyData(str(path))

    assert reference.GetPoint(0) == pytest.approx((-1.0, -2.0, 3.0))
