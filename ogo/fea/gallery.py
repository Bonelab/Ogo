"""Standalone HTML review gallery and portable inclusion CSV export."""

import csv
import hashlib
import html
import json
import os
import re
import shutil
import zipfile
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.parse import quote

from ogo.fea.review import FIELDS, REASONS, identity, inclusion_row, write_inclusion
from ogo.fea.qc_images import save_review_image


def _camera_images(path, output):
    """Split the fixed three-row Ogo render, caching assets beside the gallery."""
    from PIL import Image
    with Image.open(path) as image:
        if image.size != (1400, 2550) or not path.stem.endswith(('_qc_3d', '_sed_3d')):
            return []
        stat = path.stat()
        key = hashlib.sha256(f'{path.resolve()}:{stat.st_size}:{stat.st_mtime_ns}'.encode()).hexdigest()[:20]
        directory = output.parent / (output.stem + '_views')
        directory.mkdir(exist_ok=True)
        assets = []
        for index, view in enumerate(('oblique', 'top', 'bottom')):
            asset = directory / f'{key}_{view}.webp'
            if not asset.exists():
                save_review_image(image.crop((0, index * 850, 1400, (index + 1) * 850)), asset)
            assets.append((view, asset))
        return assets


def generate_gallery(rows, output, decisions=None, image_root=None):
    output = Path(output).absolute()
    rows = [dict(row) for row in rows]
    outcome_keys = ('stiffness_N_per_mm', 'reaction_force_N', 'pistoia_failure_load_N',
                    'masked_pistoia_failure_load_N')
    for row in rows:
        if not row.get('model_path') or row.get('generation_status') in ('too_short', 'construction_failed'):
            continue
        model = Path(row['model_path'])
        result = model.with_name(model.stem + '_results.csv')
        if result.is_file():
            with result.open(newline='') as stream:
                records = list(csv.DictReader(stream))
            if len(records) != 1:
                raise ValueError(f'Expected one model-result row: {result}')
            for key in outcome_keys:
                row.setdefault(key, records[0].get(key, ''))
    decisions = decisions or {}
    image_groups = []
    escape = lambda value: html.escape(str(value), quote=True)
    for index, row in enumerate(rows):
        images = []
        for stage, key in [('anatomy', 'anatomy_image'), ('sed', 'sed_image')]:
            path = Path(row.get(key) or '')
            if image_root is not None and not path.is_absolute():
                path = Path(image_root) / path
            if path.is_file():
                cameras = _camera_images(path, output)
                for view, asset in [('all', path), *cameras]:
                    relative = quote(Path(os.path.relpath(asset, output.parent)).as_posix())
                    from PIL import Image
                    directory = output.parent / (output.stem + '_thumbnails')
                    directory.mkdir(exist_ok=True)
                    stat = asset.stat()
                    key = hashlib.sha256(f'{asset.resolve()}:{stat.st_size}:{stat.st_mtime_ns}'.encode()).hexdigest()[:20]
                    thumbnail = directory / f'{key}.webp'
                    if not thumbnail.exists():
                        with Image.open(asset) as image:
                            image.thumbnail((420, 560))
                            image.convert('RGB').save(thumbnail, 'WEBP', quality=75, method=2)
                    images.append({'stage': stage, 'view': view, 'src': relative,
                                   'thumbnail': quote(Path(os.path.relpath(thumbnail, output.parent)).as_posix()),
                                   'fallback': 'true' if not cameras else ''})
        image_groups.append(images)
    reasons = sorted({reason for row in rows for reason in row['qc_reasons'].split(';') if reason})
    reasons = ['manual_review_pending', 'too_short'] + [r for r in reasons if r not in ('manual_review_pending', 'too_short')]
    labels = {'manual_review_pending': 'Not manually reviewed', 'too_short': 'Too short'}
    options = ''.join(f'<label><input type="checkbox" data-qc-reason="{escape(reason)}"> {escape(labels.get(reason, reason.replace("_", " ")))}</label>' for reason in reasons)
    manual_options = ''.join(f'<option value="{reason}">{reason.replace("_", " ")}</option>' for reason in REASONS)
    sort_options = '<option value="patient">Patient ID</option><option value="flags">Most QC flags first</option>'
    for key, name in [('stiffness_N_per_mm', 'Stiffness'), ('reaction_force_N', 'Reaction force'),
                      ('pistoia_failure_load_N', 'Full-bone failure load'),
                      ('masked_pistoia_failure_load_N', 'Regional failure load')]:
        for direction, label in [('asc', 'lowest first'), ('desc', 'highest first')]:
            sort_options += f'<option value="{key}:{direction}">{name}: {label}</option>'
    seeds = [inclusion_row(row, decisions.get(identity(row))) for row in rows]
    serialized = json.dumps(rows, sort_keys=True)
    fingerprint = hashlib.sha256(serialized.encode()).hexdigest()
    columns, column_ids, values = [], {}, []
    for row in rows:
        keys = tuple(row)
        if keys not in column_ids:
            column_ids[keys] = len(columns)
            columns.append(keys)
        values.append([column_ids[keys], list(row.values())])
    data = json.dumps({'rowColumns': columns, 'rowValues': values, 'images': image_groups,
                       'seeds': seeds, 'fields': FIELDS, 'reasons': REASONS,
                       'storageKey': 'ogo-fea-review-v1-' + fingerprint,
                       'imported': bool(decisions)}).replace('<', '\\u003c')
    script = Path(__file__).with_name('gallery_review.js').read_text()
    content = '''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>FEA QC review</title>
<style>body{font:14px system-ui;margin:20px;background:#fff;color:#222}header{position:sticky;top:0;background:white;padding:12px 0;z-index:1;display:flex;flex-wrap:wrap;gap:8px;align-items:center}main{display:grid;grid-template-columns:repeat(auto-fit,minmax(280px,1fr));gap:16px}article{border:1px solid #ccc;padding:12px;overflow-wrap:anywhere}button,select,input{font:inherit;padding:6px}img{width:100%;height:560px;object-fit:contain;cursor:pointer}article[hidden],img[hidden]{display:none}table{font-size:12px;width:100%;table-layout:fixed}th,td{text-align:left;border-bottom:1px solid #ddd;padding:4px;overflow-wrap:anywhere}.decision{font-weight:600}dialog{width:min(900px,90vw);max-height:90vh;overflow:auto;border:1px solid #aaa}dialog::backdrop{background:#0008}dialog img{height:60vh}dialog label{display:block;margin:12px 0}textarea{display:block;width:95%;min-height:70px}#message{width:100%;margin:0;color:#555}.review-open:first-child{font-weight:600;border:0;background:none;text-align:left}#review-title{font-size:18px}</style></head><body>
<header><input id="search" aria-label="Subject" placeholder="Subject"><select id="site" aria-label="Anatomy"><option value="">All sites</option><option>hip</option><option>spine</option></select>
<select id="status" aria-label="Automatic QC"><option value="">All QC statuses</option><option>pass</option><option>review</option><option>fail</option></select>
<select id="sort" aria-label="Sort models">SORT_OPTIONS</select>
<details id="reason-picker"><summary id="reason-summary">Exclude flags</summary><div class="flag-menu"><button id="clear-flags" type="button">Clear</button>REASONS</div></details>
<select id="inclusion" aria-label="Study inclusion"><option value="">All study decisions</option><option>include</option><option>exclude</option><option>pending</option></select>
<select id="stage" aria-label="Solve stage"><option value="anatomy">Before solve</option><option value="sed">After solve</option></select>
<select id="camera" aria-label="Camera view"><option value="oblique">Oblique / side</option><option value="top">Top</option><option value="bottom">Bottom</option><option value="all">All views</option></select>
<input id="reviewer" aria-label="Reviewer" placeholder="Reviewer initials/name"><button id="export">Export inclusion CSV</button><label>Import CSV <input id="import" type="file" accept=".csv,text/csv"></label><p id="message" role="status"></p>
<nav aria-label="Pages"><label>Per page <select id="page-size"><option value="50" selected>50</option><option value="100">100</option><option value="250">250</option><option value="500">500</option></select></label> <button id="previous-page" title="Previous page">&#8592;</button> <span id="page-status"></span> <button id="next-page" title="Next page">&#8594;</button></nav></header>
<main></main><dialog id="review-dialog"><button id="close">Close</button><h1 id="review-title"></h1><p id="review-auto"></p>
<select id="review-stage" aria-label="Review solve stage"><option value="anatomy">Before solve</option><option value="sed">After solve</option></select>
<select id="review-camera" aria-label="Review camera view"><option value="oblique">Oblique / side</option><option value="top">Top</option><option value="bottom">Bottom</option><option value="all">All views</option></select><div id="review-images"></div>
<label>Reason <select id="manual-reason">MANUAL_OPTIONS</select></label><label>Explanation<textarea id="manual-note"></textarea></label>
<button id="include">Include</button><button id="exclude">Exclude</button><button id="reset">Reset to automatic</button><div id="review-measurements"></div></dialog>
<script type="application/json" id="review-data">DATA</script><script>SCRIPT</script></body></html>'''
    # Replace once so record text cannot become another template token.
    tokens = {'SCRIPT': script, 'MANUAL_OPTIONS': manual_options, 'REASONS': options, 'SORT_OPTIONS': sort_options,
              'DATA': data}
    content = re.sub(r'SCRIPT|MANUAL_OPTIONS|REASONS|SORT_OPTIONS|DATA', lambda match: tokens[match[0]], content)
    menu_style = '#reason-picker{position:relative}#reason-picker summary{cursor:pointer;border:1px solid #999;padding:6px}#reason-picker .flag-menu{position:absolute;top:100%;left:0;width:320px;max-width:85vw;max-height:55vh;overflow:auto;background:white;border:1px solid #aaa;padding:10px;z-index:3}#reason-picker label{display:flex;gap:6px;padding:5px 0}'
    content = content.replace('</style>', menu_style + '</style>', 1)
    output.write_text(content)
    return output


def export_gallery_zip(rows, output, decisions=None):
    """Package the full gallery, images, and review CSVs without copying models.

    Image references are relative to the extracted HTML. Missing previews stay
    missing; the archive contains no dependence on the original image folders.
    """
    output = Path(output).absolute()
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix='ogo-gallery-', dir=output.parent) as temporary:
        root = Path(temporary)
        images = root / 'images'
        images.mkdir()
        portable, copied = [], {}
        for row in rows:
            record = dict(row)
            for key in ('anatomy_image', 'sed_image'):
                source = Path(row.get(key) or '')
                record[key] = ''
                if source.is_file():
                    source = source.resolve()
                    if source not in copied:
                        # Stable identity preserves browser autosave across exports.
                        name = hashlib.sha256(str(source).encode()).hexdigest()[:20] + '_' + source.name
                        target = images / name
                        if source.suffix.lower() == '.png':
                            from PIL import Image
                            target = target.with_suffix('.webp')
                            with Image.open(source) as image:
                                save_review_image(image, target)
                        else:
                            try:
                                target.symlink_to(source)
                            except OSError:
                                shutil.copyfile(source, target)
                        copied[source] = target.relative_to(root).as_posix()
                    record[key] = copied[source]
            portable.append(record)
        generate_gallery(portable, root / 'gallery.html', decisions, image_root=root)
        write_inclusion(rows, root / 'study_inclusion.csv', decisions)
        fields = sorted({key for row in rows for key in row})
        with (root / 'qc_summary.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
        temporary_zip = root / 'gallery.zip'
        with zipfile.ZipFile(temporary_zip, 'w', compression=zipfile.ZIP_DEFLATED, allowZip64=True) as archive:
            for source in sorted(root.rglob('*')):
                if source.is_file() and source != temporary_zip:
                    compression = zipfile.ZIP_STORED if source.suffix.lower() in ('.png', '.jpg', '.jpeg', '.webp') else zipfile.ZIP_DEFLATED
                    archive.write(source, source.relative_to(root).as_posix(), compress_type=compression)
        temporary_zip.replace(output)
    return output
