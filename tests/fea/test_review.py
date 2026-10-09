import csv
import json
from pathlib import Path
import shutil
import subprocess

import pytest

from ogo.fea.review import inclusion_row, load_decisions


def test_gallery_has_quick_decisions_on_each_card(tmp_path):
    from ogo.fea.gallery import generate_gallery
    output = tmp_path / 'gallery.html'
    generate_gallery([{'site': 'hip', 'model_id': 'case', 'qc_status': 'fail',
                       'qc_reasons': 'incomplete_distal_boundary'}], output)
    content = output.read_text()
    assert 'class="quick-decision" data-decision="include"' in content
    assert 'class="quick-decision" data-decision="exclude"' in content
    assert "manual_reason: 'manual_review'" in content
    assert "window.alert('CSV export started." in content
    assert 'does not update the project CSV' in content


def test_manual_decision_preserves_automatic_flags():
    row = inclusion_row({'model_id': 'case', 'site': 'hip', 'qc_status': 'fail',
                         'qc_reasons': 'too_short', 'generation_status': 'too_short'},
                        {'manual_decision': 'include', 'reviewer': 'MW'})
    assert row['final_inclusion'] == 'include'
    assert row['automatic_qc_status'] == 'fail'
    assert row['solver_eligible'] == 'False'
    assert row['manual_reason'] == 'manual_review'


def test_review_is_pending_and_auto_reset_has_no_manual_reason():
    row = inclusion_row({'model_id': 'case', 'site': 'spine', 'qc_status': 'review'}, {})
    assert row['final_inclusion'] == 'pending'
    assert row['manual_reason'] == ''


def test_import_requires_unambiguous_identity(tmp_path):
    path = tmp_path / 'decisions.csv'
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=['model_id', 'site', 'manual_decision'])
        writer.writeheader()
        writer.writerows([{'model_id': 'case', 'site': 'hip', 'manual_decision': 'include'}] * 2)
    with pytest.raises(ValueError, match='Duplicate'):
        load_decisions(path)


def test_command_backfills_only_missing_measurements(tmp_path, monkeypatch):
    from ogo.cli import ValidateFEA
    from ogo.fea.validation import write_measurements
    write_measurements(tmp_path / 'existing.n88model', {'site': 'hip'})
    (tmp_path / 'existing.n88model').touch()
    (tmp_path / 'legacy.n88model').touch()
    calls = []
    def backfill(path, site=None, evidence=None):
        from ogo.fea.legacy_qc import evidence_signature
        calls.append(path.name)
        return write_measurements(path, {'site': 'spine', 'measurement_source': 'legacy_backfill',
                                        'legacy_evidence_signature': evidence_signature(path, evidence)})
    monkeypatch.setattr(ValidateFEA, 'backfill_measurements', backfill)
    ValidateFEA.main([str(tmp_path), '--output', str(tmp_path / 'audit'), '--gallery'])
    assert calls == ['legacy.n88model']
    assert (tmp_path / 'audit' / 'study_inclusion.csv').exists()
    calls.clear()
    ValidateFEA.main([str(tmp_path), '--output', str(tmp_path / 'audit')])
    assert calls == []


def test_javascript_csv_roundtrip_and_manual_provenance():
    node = shutil.which('node')
    if not node:
        pytest.skip('Node is required for browser-independent JavaScript tests')
    from ogo.fea import gallery
    script = Path(gallery.__file__).with_name('gallery_review.js')
    code = '''const assert = require('node:assert/strict');
const {inclusionRecord, encodeCSV, parseCSV, validateReviews, imageVisible} = require(SCRIPT);
assert.equal(imageVisible({stage:'anatomy',view:'top'},'anatomy','top'),true);
assert.equal(imageVisible({stage:'anatomy',view:'bottom'},'anatomy','top'),false);
assert.equal(imageVisible({stage:'sed',view:'top'},'anatomy','top'),false);
assert.equal(imageVisible({stage:'sed',view:'all'},'sed','all'),true);
assert.equal(imageVisible({stage:'sed',view:'all',fallback:'true'},'sed','top'),true);
const row = {site:'hip', model_id:'test', qc_status:'fail', qc_reasons:'too_short', generation_status:'too_short'};
const decision = inclusionRecord(row, {manual_decision:'include', manual_note:'Quote "here",\nnew line'});
assert.equal(decision.automatic_qc_status, 'fail');
assert.equal(decision.solver_eligible, 'False');
assert.equal(decision.manual_reason, 'manual_review');
const records = parseCSV(encodeCSV([decision], Object.keys(decision)));
assert.equal(records[0].manual_note, decision.manual_note);
assert.equal(validateReviews(records,[row]).size,1);
assert.throws(()=>validateReviews([records[0],records[0]],[row]),/Duplicate/);
assert.throws(()=>validateReviews(records,[]),/Unknown/);
assert.equal(inclusionRecord(row,{manual_decision:'automatic'}).final_inclusion,'exclude');
'''.replace('SCRIPT', json.dumps(str(script)))
    # Keep a newline inside the JS string escaped, not an invalid literal.
    code = code.replace('Quote "here",\nnew line', 'Quote "here",\\nnew line')
    subprocess.run([node, '-e', code], check=True, capture_output=True, text=True)


def test_gallery_does_not_interpret_template_words_in_subject(tmp_path):
    from ogo.fea.validation import evaluate, generate_gallery
    row = evaluate({'site': 'hip', 'model_id': 'DATA_SCRIPT_CARDS', 'qc_reasons': ''})
    path = generate_gallery([row], tmp_path / 'review.html')
    assert 'DATA_SCRIPT_CARDS' in path.read_text()


def test_legacy_anatomy_is_not_silently_certified():
    from ogo.fea.validation import evaluate
    row = evaluate({'site': 'spine', 'measurement_source': 'legacy_backfill'})
    assert 'legacy_body_process_metrics_missing' in row['qc_reasons']


def test_gallery_individual_views_keep_panel_order_and_cache(tmp_path):
    from PIL import Image
    from ogo.fea.gallery import generate_gallery
    image = tmp_path / 'case_qc_3d.png'
    panel = Image.new('RGB', (1400, 2550))
    for index, color in enumerate(('red', 'green', 'blue')):
        panel.paste(color, (0, index * 850, 1400, (index + 1) * 850))
    panel.save(image)
    row = {'site': 'hip', 'model_id': 'case', 'qc_status': 'review',
           'qc_reasons': '', 'anatomy_image': str(image)}
    path = generate_gallery([row], tmp_path / 'gallery.html')
    content = path.read_text()
    assert 'id="camera"' in content
    assert 'id="review-camera"' in content
    assert 'id="review-stage"' in content
    for view, color in zip(('oblique', 'top', 'bottom'), ((255, 0, 0), (0, 128, 0), (0, 0, 255))):
        assert f'data-view="{view}"' in content
        asset = next((tmp_path / 'gallery_views').glob(f'*_{view}.webp'))
        with Image.open(asset) as cropped:
            assert cropped.size == (1400, 850)
            assert max(abs(a - b) for a, b in zip(cropped.getpixel((700, 425)), color)) <= 5
        modified = asset.stat().st_mtime_ns
        generate_gallery([row], path)
        assert asset.stat().st_mtime_ns == modified


def test_gallery_unknown_montage_is_not_split(tmp_path):
    from PIL import Image
    from ogo.fea.gallery import generate_gallery
    image = tmp_path / 'legacy.png'
    Image.new('RGB', (400, 400), 'white').save(image)
    row = {'site': 'spine', 'model_id': 'legacy', 'qc_status': 'review',
           'qc_reasons': '', 'anatomy_image': str(image)}
    content = generate_gallery([row], tmp_path / 'gallery.html').read_text()
    assert 'data-fallback="true"' in content
    assert not (tmp_path / 'gallery_views').exists()


def test_gallery_zip_is_portable_and_has_review_records(tmp_path):
    import zipfile
    from PIL import Image
    from ogo.fea.gallery import export_gallery_zip
    image = tmp_path / 'case_qc_3d.png'
    Image.new('RGB', (1400, 2550), 'red').save(image)
    row = {'site': 'hip', 'model_id': 'case', 'qc_status': 'review',
           'qc_reasons': '', 'anatomy_image': str(image), 'sed_image': str(image)}
    path = export_gallery_zip([row], tmp_path / 'gallery.zip')
    with zipfile.ZipFile(path) as archive:
        assert {'gallery.html', 'qc_summary.csv', 'study_inclusion.csv'} <= set(archive.namelist())
        assert len([name for name in archive.namelist() if name.startswith('images/')]) == 1
        assert len([name for name in archive.namelist() if name.startswith('gallery_views/')]) == 3
        assert not any(name.endswith('.png') for name in archive.namelist())
        assert archive.testzip() is None
        content = archive.read('gallery.html').decode()
        assert 'src="images/' in content
        assert str(tmp_path) not in content
        assert 'data:image' not in content
        archive.extractall(tmp_path / 'unpacked')
    import re
    for source in re.findall(r'src="([^"]+)"', content):
        assert (tmp_path / 'unpacked' / source).is_file()
