import pytest


def _vtk_image_from_array(data, *, origin=(0, 0, 0), spacing=(1, 1, 1)):
    vtk = pytest.importorskip("vtk")
    from vtk.util.numpy_support import numpy_to_vtk

    image = vtk.vtkImageData()
    image.SetDimensions(data.shape)
    image.SetOrigin(*origin)
    image.SetSpacing(*spacing)
    image.GetPointData().SetScalars(numpy_to_vtk(data.ravel(order="F"), deep=True))
    return image


class _GridPointModel:
    def __init__(self, shape, *, origin=(0, 0, 0), spacing=(1, 1, 1)):
        self.points = []
        for x in range(shape[0] + 1):
            for y in range(shape[1] + 1):
                for z in range(shape[2] + 1):
                    self.points.append(
                        (
                            origin[0] + (x - 0.5) * spacing[0],
                            origin[1] + (y - 0.5) * spacing[1],
                            origin[2] + (z - 0.5) * spacing[2],
                        )
                    )

    def GetNumberOfPoints(self):
        return len(self.points)

    def GetPoint(self, point_id):
        return self.points[point_id]


def test_interface_node_ids_from_voxel_mask_uses_material_interface_only():
    np = pytest.importorskip("numpy")

    from ogo.fea.model import interface_node_ids_from_voxel_mask

    material = np.ones((2, 2, 2), dtype=np.uint8)
    selected = np.zeros_like(material)
    selected[:, :, 0] = 1
    model = _GridPointModel(material.shape)

    node_ids = interface_node_ids_from_voxel_mask(
        model,
        _vtk_image_from_array(selected),
        _vtk_image_from_array(material),
        name="shaft",
    )
    points = {
        model.GetPoint(node_ids.GetValue(index))
        for index in range(node_ids.GetNumberOfTuples())
    }

    assert node_ids.GetName() == "shaft"
    assert node_ids.GetNumberOfTuples() == 9
    assert {point[2] for point in points} == {0.5}


def test_interface_node_ids_from_voxel_mask_can_follow_direction():
    np = pytest.importorskip("numpy")

    from ogo.fea.model import interface_node_ids_from_voxel_mask

    material = np.ones((2, 2, 2), dtype=np.uint8)
    selected = np.zeros_like(material)
    selected[0, 0, 0] = 1
    model = _GridPointModel(material.shape)

    node_ids = interface_node_ids_from_voxel_mask(
        model,
        _vtk_image_from_array(selected),
        _vtk_image_from_array(material),
        name="shaft",
        direction=(0.0, 0.0, 1.0),
    )
    points = {
        model.GetPoint(node_ids.GetValue(index))
        for index in range(node_ids.GetNumberOfTuples())
    }

    assert node_ids.GetNumberOfTuples() == 4
    assert {point[2] for point in points} == {0.5}


def test_directional_face_node_ids_from_voxel_mask_uses_outer_face():
    np = pytest.importorskip("numpy")

    from ogo.fea.model import directional_face_node_ids_from_voxel_mask

    selected = np.ones((2, 2, 2), dtype=np.uint8)
    model = _GridPointModel(selected.shape)

    node_ids = directional_face_node_ids_from_voxel_mask(
        model,
        _vtk_image_from_array(selected),
        direction=(0.0, 1.0, 0.0),
        name="cap",
    )
    points = {
        model.GetPoint(node_ids.GetValue(index))
        for index in range(node_ids.GetNumberOfTuples())
    }

    assert node_ids.GetName() == "cap"
    assert node_ids.GetNumberOfTuples() == 9
    assert {point[1] for point in points} == {1.5}


def test_directional_face_node_ids_from_label_image_uses_requested_label_only():
    np = pytest.importorskip("numpy")

    from ogo.fea.model import directional_face_node_ids_from_label_image

    labels = np.zeros((2, 2, 3), dtype=np.uint8)
    labels[:, :, 1] = 4
    labels[:, :, 2] = 3
    model = _GridPointModel(labels.shape)

    node_ids = directional_face_node_ids_from_label_image(
        model,
        _vtk_image_from_array(labels),
        label_value=4,
        direction=(0.0, 0.0, 1.0),
        name="top_cap",
    )
    points = {
        model.GetPoint(node_ids.GetValue(index))
        for index in range(node_ids.GetNumberOfTuples())
    }

    assert node_ids.GetName() == "top_cap"
    assert node_ids.GetNumberOfTuples() == 9
    assert {point[2] for point in points} == {1.5}
