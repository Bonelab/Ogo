"""Audit FEA measurements, backfill legacy models, and optionally create a gallery."""

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from ogo.fea.qc.validation import evaluate, generate_gallery
from ogo.fea.qc.validation import write_measurements
from ogo.fea.qc.gallery import export_gallery_zip
from ogo.fea.qc.legacy import backfill_measurements, evidence_signature, load_evidence, needs_backfill
from ogo.fea.qc.review import identity, load_decisions, write_inclusion


def _backfill_task(task):
    model, site, evidence = task
    try:
        return backfill_measurements(model, site=site, evidence=evidence), ''
    except Exception as exc:
        metrics = write_measurements(model, {
            'site': site or 'unknown', 'measurement_source': 'legacy_backfill_failed',
            'measurement_error': str(exc), 'generation_status': 'measurement_failed',
            'legacy_evidence_signature': evidence_signature(model, evidence)})
        return metrics, str(exc)


def main(argv=None):
    parser = argparse.ArgumentParser(prog='ogoValidateFEA', description=__doc__)
    parser.add_argument('input', type=Path, help='Model directory, .n88model, or *_qc_metrics.csv')
    parser.add_argument('--output', type=Path, default=Path('qc'), help='Audit output directory')
    parser.add_argument('--gallery', action='store_true', help='Also generate gallery.html')
    parser.add_argument('--gallery-zip', action='store_true', help='Export portable gallery.zip with HTML, images, and QC/inclusion CSVs')
    parser.add_argument('--site', choices=('hip', 'spine'), help='Anatomy for legacy models when it cannot be inferred')
    parser.add_argument('--reviews', type=Path, help='Import a previously exported study_inclusion.csv')
    parser.add_argument('--evidence', type=Path, help='CSV of registered body/process masks and shaft sidecars by site/model_id')
    parser.add_argument('--workers', type=int, default=1, help='Parallel legacy backfill processes (default: 1)')
    args = parser.parse_args(argv)
    if not 1 <= args.workers <= 32:
        parser.error('--workers must be between 1 and 32')
    try:
        evidence = load_evidence(args.evidence) if args.evidence else {}
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    if not args.input.exists():
        parser.error(f'Input does not exist: {args.input}')
    if args.input.is_dir():
        files = list(args.input.rglob('*_qc_metrics.csv'))
        models = sorted(args.input.rglob('*.n88model'))
    elif args.input.suffix == '.n88model':
        files, models = [], [args.input]
    elif args.input.name.endswith('_qc_metrics.csv'):
        files, models = [args.input], []
    else:
        parser.error('Input must be a directory, .n88model, or *_qc_metrics.csv')
    tasks = []
    for model in models:
        metrics = model.with_name(model.stem + '_qc_metrics.csv')
        row = {}
        if metrics.exists():
            with metrics.open(newline='') as stream:
                row = next(csv.DictReader(stream), {})
        site = args.site or row.get('site')
        matches = [key for key in evidence if key[1] == model.stem and (not site or key[0] == site)]
        if len(matches) > 1:
            parser.error(f'Ambiguous evidence identity for {model}; specify --site')
        sources = evidence[matches[0]] if matches else None
        if matches:
            site = matches[0][0]
        anatomy_changed = row.get('measurement_source', '').startswith('legacy_backfill') and site and row.get('site') != site
        if anatomy_changed or needs_backfill(model, row, sources):
            print(f'Backfilling {model}')
            tasks.append((model, site, sources))
        if metrics not in files:
            files.append(metrics)
    if args.workers == 1:
        results = map(_backfill_task, tasks)
        for metrics, error in results:
            if error:
                print(f'Review required: {metrics}: {error}')
    elif tasks:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            for metrics, error in pool.map(_backfill_task, tasks):
                if error:
                    print(f'Review required: {metrics}: {error}')
    files.sort()
    if not files:
        parser.error('No QC CSVs or n88models found')
    rows = []
    for path in files:
        with path.open(newline='') as stream:
            records = list(csv.DictReader(stream))
        if len(records) != 1 or not {'model_id', 'site'} <= records[0].keys():
            parser.error(f'Invalid measurement CSV: {path}')
        rows.append(evaluate(records[0]))
    keys = [identity(row) for row in rows]
    if len(set(keys)) != len(keys):
        parser.error('Duplicate site/model_id identities. Audit each cohort or length run separately.')
    try:
        decisions = load_decisions(args.reviews) if args.reviews else {}
    except ValueError as exc:
        parser.error(str(exc))
    unknown = set(decisions) - set(keys)
    if unknown:
        parser.error(f'Review CSV contains {len(unknown)} identities not present in this audit')
    args.output.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in rows for key in row})
    with (args.output / 'qc_summary.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    write_inclusion(rows, args.output / 'study_inclusion.csv', decisions)
    if args.gallery:
        generate_gallery(rows, args.output / 'gallery.html', decisions=decisions)
    if args.gallery_zip:
        export_gallery_zip(rows, args.output / 'gallery.zip', decisions)
    for status in ('pass', 'review', 'fail'):
        print(f'{status}: {sum(row["qc_status"] == status for row in rows)}')


if __name__ == '__main__':
    main()
