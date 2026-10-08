import csv
import json

import numpy as np
import pytest

from tests.fea.test_model_export import hip_model


def test_evidence_cache_retries_only_when_inputs_change(tmp_path):
    from ogo.fea.legacy_qc import evidence_signature, needs_backfill
    model = tmp_path / 'case.n88model'
    model.touch()
    row = {'measurement_source': 'legacy_backfill'}
    assert needs_backfill(model, row)
    row['legacy_evidence_signature'] = evidence_signature(model)
    assert not needs_backfill(model, row)
    model.with_name('case_shaft_geometry.json').write_text('{}')
    assert needs_backfill(model, row)
    assert not needs_backfill(model, {'measurement_source': 'generation'})


def test_evidence_manifest_resolves_relative_paths_and_rejects_duplicates(tmp_path):
    from ogo.fea.legacy_qc import load_evidence
    path = tmp_path / 'sources.csv'
    path.write_text('site,model_id,body_mask\nspine,case,body.nii.gz\n')
    assert load_evidence(path)[('spine', 'case')]['body_mask'] == str(tmp_path / 'body.nii.gz')
    path.write_text(path.read_text() + 'spine,case,other.nii.gz\n')
    with pytest.raises(ValueError, match='Duplicate'):
        load_evidence(path)


def test_registered_mask_geometry_and_process_recovery(tmp_path):
    sitk = pytest.importorskip('SimpleITK')
    from ogo.fea.legacy_qc import registered_anatomy
    mask = sitk.GetImageFromArray(np.ones((2, 2, 2), dtype=np.uint8))
    body = tmp_path / 'body.nii.gz'
    sitk.WriteImage(mask, str(body))
    grid = np.ones((3, 2, 2), dtype=np.int32)
    images, provenance = registered_anatomy(grid, np.ones(3), np.zeros(3), body)
    from ogo.util.vtk_image import vtk_image_to_numpy
    assert vtk_image_to_numpy(images[0]).sum() == 8
    assert vtk_image_to_numpy(images[1]).sum() == 4
    assert provenance['process_mask_source'] == 'model_bone_minus_registered_body'
    mask.SetOrigin((100, 0, 0))
    sitk.WriteImage(mask, str(body))
    with pytest.raises(ValueError, match='overlap'):
        registered_anatomy(grid, np.ones(3), np.zeros(3), body)


def test_registered_mask_must_match_voxel_lattice(tmp_path):
    sitk = pytest.importorskip('SimpleITK')
    from ogo.fea.legacy_qc import registered_anatomy
    mask = sitk.GetImageFromArray(np.ones((2, 2, 2), dtype=np.uint8))
    mask.SetOrigin((0.3, 0, 0))
    path = tmp_path / 'body.nii.gz'
    sitk.WriteImage(mask, str(path))
    with pytest.raises(ValueError, match='lattice'):
        registered_anatomy(np.ones((3, 2, 2), dtype=int), np.ones(3), np.zeros(3), path)


def test_zero_anatomy_offset_is_valid_evidence():
    from ogo.fea.validation import evaluate
    row = evaluate({'site': 'spine', 'measurement_source': 'legacy_backfill',
                    'process_body_axial_offset_mm': 0,
                    'process_body_transverse_offset_mm': 10})
    assert 'legacy_body_process_metrics_missing' not in row['qc_reasons']


def test_hip_distal_set_is_recovered_from_saved_plane_not_union(tmp_path, monkeypatch, hip_model):
    import vtkbone
    from ogo.fea.legacy_qc import backfill_measurements
    model = hip_model
    model.ApplyBoundaryCondition('Femoral_Head_PMMA_Nodes', 1, -.1, 'top_displacement')
    model.ApplyBoundaryCondition('Greater_Trochanter_PMMA_Nodes', 1, 0, 'bottom_fixed_y_PMMA')
    model.ApplyBoundaryCondition('Greater_Trochanter_PMMA_Nodes', 2, 0, 'bottom_fixed_z')
    for sense, name in ((0, 'bottom_fixed_x'), (2, 'bottom_fixed_z')):
        model.ApplyBoundaryCondition('Distal_Femur_Nodes', sense, 0, name)
    model.GetNodeSets().RemoveItem(model.GetNodeSet('Distal_Femur_Nodes'))
    class Reader:
        def SetFileName(self, name):
            pass
        def Update(self):
            pass
        def GetOutput(self):
            return model
    monkeypatch.setattr(vtkbone, 'vtkboneN88ModelReader', Reader)
    path = tmp_path / 'hip.n88model'
    path.touch()
    path.with_name('hip_shaft_geometry.json').write_text(json.dumps({'cut_z_mm': model.GetBounds()[4]}))
    row = next(csv.DictReader(backfill_measurements(path, site='hip').open()))
    assert 'Distal_Femur_Nodes' in row['recovered_boundary_sets']
    assert float(row['distal_boundary_span_mm']) == 0
    assert float(row['distal_boundary_coverage_fraction']) == 1
    assert float(row['distal_support_patch_coverage_fraction']) == 1
    assert float(row['distal_support_patch_width_fraction']) == 0.9
    assert row['shaft_bc_node_count'] == row['shaft_bc_fixed_z_node_count']
