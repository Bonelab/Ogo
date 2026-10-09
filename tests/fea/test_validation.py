import csv

import numpy as np

from ogo.fea.validation import evaluate, grid_measurements, write_measurements, generate_gallery


def test_disk_contact_uses_faces_not_diagonal_neighbors():
    labels = np.zeros((4, 4, 4), dtype=int)
    labels[1, 1, 1] = 1
    labels[1, 1, 2] = 5000
    labels[2, 2, 2] = 5000
    result = grid_measurements(labels, (2, 2, 2), loading_axis=2)
    assert result['disk_1_contact_area_mm2'] == 4
    assert result['disk_component_count'] == 2


def test_too_short_and_missing_checks_are_not_passes():
    assert evaluate({'site': 'hip', 'generation_status': 'too_short'})['qc_status'] == 'fail'
    assert evaluate({'site': 'spine'})['qc_status'] == 'review'


def test_flat_boundary_and_length_checks():
    row = {'site': 'hip', 'generation_status': 'generated', 'bone_component_count': '1',
           'disk_component_count': '2', 'boundary_count': '3', 'loading_compressive': 'True',
           'loading_axis_matches': 'True',
           'shaft_bc_node_count': '100', 'shaft_bc_fixed_x_node_count': '100',
           'shaft_bc_fixed_z_node_count': '100',
           'measured_shaft_length_mm': '10', 'requested_shaft_length_mm': '10',
           'complete_distal_section_verified': 'True', 'distal_boundary_span_mm': '0',
           'distal_boundary_coverage_fraction': '1'}
    for component in (1, 2):
        row[f'disk_{component}_contact_area_mm2'] = '10'
        row[f'disk_{component}_contact_component_count'] = '1'
    assert evaluate(row)['qc_status'] == 'pass'
    row['measured_shaft_length_mm'] = '8'
    assert 'shaft_length_mismatch' in evaluate(row)['qc_reasons']


def test_partial_shaft_fixation_is_flagged():
    row = {'site': 'hip', 'shaft_bc_node_count': 100,
           'shaft_bc_fixed_x_node_count': 99, 'shaft_bc_fixed_z_node_count': 100}
    assert 'incomplete_shaft_fixation_x' in evaluate(row)['qc_reasons']
    assert 'incomplete_shaft_fixation_z' not in evaluate(row)['qc_reasons']


def test_intended_distal_patch_coverage_takes_precedence():
    row = {'site': 'hip', 'distal_boundary_coverage_fraction': 0.945,
           'distal_support_patch_coverage_fraction': 1.0}
    assert 'incomplete_distal_boundary' not in evaluate(row)['qc_reasons']
    row['distal_support_patch_coverage_fraction'] = 0.98
    assert 'incomplete_distal_boundary' in evaluate(row)['qc_reasons']


def test_legacy_full_face_coverage_tolerates_existing_patch():
    row = {'site': 'hip', 'distal_boundary_coverage_fraction': 0.89}
    assert 'incomplete_distal_boundary' not in evaluate(row)['qc_reasons']
    row['distal_boundary_coverage_fraction'] = 0.88
    assert 'incomplete_distal_boundary' in evaluate(row)['qc_reasons']


def test_csv_and_gallery_escape_subject_names(tmp_path):
    path = tmp_path / 'case.n88model'
    csv_path = write_measurements(path, {'site': 'hip', 'generation_status': 'too_short'})
    with csv_path.open() as stream:
        row = next(csv.DictReader(stream))
    row['model_id'] = '<script>bad</script>'
    output = tmp_path / 'gallery.html'
    generate_gallery([evaluate(row)], output)
    content = output.read_text()
    assert '\\u003cscript>bad\\u003c/script>' in content
    assert '<script>bad</script>' not in content
    assert 'too_short' in content
    assert 'loading="lazy"' not in content  # No fabricated image for a missing preview.


def test_command_uses_csv_only(tmp_path):
    from ogo.cli.ValidateFEA import main
    write_measurements(tmp_path / 'missing.n88model', {'site': 'hip', 'generation_status': 'too_short'})
    main([str(tmp_path), '--output', str(tmp_path / 'audit'), '--gallery'])
    assert (tmp_path / 'audit' / 'gallery.html').exists()
    with (tmp_path / 'audit' / 'qc_summary.csv').open() as stream:
        assert next(csv.DictReader(stream))['qc_status'] == 'fail'


def test_solved_measurements_do_not_replace_geometry(tmp_path):
    path = tmp_path / 'model.n88model'
    csv_path = write_measurements(path, {'site': 'hip', 'available_shaft_length_mm': 18})
    write_measurements(path, {'sed_finite': False})
    with csv_path.open() as stream:
        row = next(csv.DictReader(stream))
    assert row['available_shaft_length_mm'] == '18'
    assert 'invalid_sed_finite' in evaluate(row)['qc_reasons']
