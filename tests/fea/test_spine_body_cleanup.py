import numpy as np
import pytest

vtk = pytest.importorskip('vtk')

from ogo.fea.spine import clean_body_component
from ogo.util.vtk_image import numpy_to_vtk_image, vtk_image_to_numpy


def image_for(data):
    template = vtk.vtkImageData()
    template.SetDimensions(*data.shape)
    template.SetSpacing(0.8, 0.9, 1.5)
    template.SetOrigin(10, -20, 30)
    return numpy_to_vtk_image(data, template, vtk_array_type=vtk.VTK_UNSIGNED_CHAR)


def test_body_cleanup_preserves_other_labels_and_geometry():
    data = np.zeros((8, 8, 8), dtype=np.uint8)
    data[1:4, 1:4, 1:4] = 20
    data[6, 6, 6] = 20
    data[4:6, 1:4, 1:4] = 48
    data[7, 0, 0] = 99
    original = image_for(data)
    cleaned, metrics = clean_body_component(original, 20)
    expected = data.copy()
    expected[6, 6, 6] = 0
    assert np.array_equal(vtk_image_to_numpy(cleaned), expected)
    assert np.array_equal(vtk_image_to_numpy(original), data)
    assert cleaned.GetSpacing() == original.GetSpacing()
    assert cleaned.GetOrigin() == original.GetOrigin()
    assert cleaned.GetExtent() == original.GetExtent()
    assert metrics['input_body_component_count'] == 2
    assert metrics['body_cleanup_removed_voxels'] == 1
    assert metrics['body_cleanup_removed_fraction'] == pytest.approx(1 / 28)
    assert metrics['body_cleanup_removed_volume_mm3'] == pytest.approx(0.8 * 0.9 * 1.5)


def test_single_body_component_is_unchanged():
    data = np.zeros((5, 5, 5), dtype=np.uint8)
    data[1:4, 1:4, 1:4] = 20
    cleaned, metrics = clean_body_component(image_for(data), 20)
    assert np.array_equal(vtk_image_to_numpy(cleaned), data)
    assert metrics['body_cleanup_removed_fraction'] == 0


def test_diagonal_body_fragments_are_not_face_connected():
    data = np.zeros((5, 5, 5), dtype=np.uint8)
    data[1:3, 1:3, 1:3] = 20
    data[3, 3, 3] = 20
    cleaned, metrics = clean_body_component(image_for(data), 20)
    assert vtk_image_to_numpy(cleaned)[3, 3, 3] == 0
    assert metrics['input_body_component_count'] == 2


def test_empty_body_is_rejected():
    with pytest.raises(ValueError, match='body label'):
        clean_body_component(image_for(np.zeros((5, 5, 5), dtype=np.uint8)), 20)


def test_substantial_body_cleanup_requires_review_not_exclusion():
    from ogo.fea.qc.validation import evaluate

    row = {'site': 'spine', 'body_cleanup_removed_fraction': 0.06}
    assert 'substantial_body_cleanup' in evaluate(row)['qc_reasons']
    row['body_cleanup_removed_fraction'] = 0.01
    assert 'substantial_body_cleanup' not in evaluate(row)['qc_reasons']
