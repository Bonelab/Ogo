"""Package selected rerun spine cases, with fresh manual-review decisions."""
import csv
import json
from pathlib import Path
import sys

from ogo.fea.gallery import export_gallery_zip
from ogo.fea.validation import evaluate


def main(root):
    with (root / 'rerun_selection.csv').open() as stream:
        selected = list(csv.DictReader(stream))
    expected = len(selected)
    rows = []
    for original in selected:
        stem = original['model_id']
        subject = stem.removesuffix('_QCT_L1')
        model = root / 'models' / subject / (stem + '.n88model')
        metric = model.with_name(stem + '_qc_metrics.csv')
        row = {}
        if metric.is_file():
            with metric.open() as stream:
                row = next(csv.DictReader(stream))
        row.update(model_id=stem, site='spine', model_path=str(model))
        status = json.loads((root / 'cases' / subject / 'status.json').read_text())
        row['rerun_stage'] = status['stage']
        row['previous_manual_decision'] = original['manual_decision']
        row['rerun_error'] = status.get('error', '')
        if status['stage'] != 'complete':
            row['generation_status'] = 'construction_failed'
        for suffix, field in [('_qc_3d.webp', 'anatomy_image'), ('_sed_3d.webp', 'sed_image')]:
            image = model.with_name(stem + suffix)
            row[field] = str(image) if image.is_file() else ''
        result = model.with_name(stem + '_results.csv')
        if status['stage'] == 'complete' and result.is_file():
            with result.open() as stream:
                row.update(next(csv.DictReader(stream)))
        row = evaluate(row)
        if status['stage'] != 'complete':
            row.update(qc_status='fail', qc_reasons='rerun_failed;' + row['qc_reasons'])
        rows.append(row)
    if len(rows) != expected:
        raise ValueError(f'Rerun gallery must retain all {expected} cases, including failures')
    export_gallery_zip(rows, root / 'after_review.zip')
    (root / 'gallery_ready.json').write_text(json.dumps({'cases': expected, 'zip': 'after_review.zip'}))
    print(f'After-only gallery packaged: {expected} cases; original decisions not imported', flush=True)


if __name__ == '__main__':
    main(Path(sys.argv[1]))
