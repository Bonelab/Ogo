"""Resegment manually excluded L1 cases without replacing selected study models.

Run on Groot, not an ARC login node. Two GPU workers segment full native CTs;
four CPU workers generate and solve the isolated L1 models. Completed stages
are resumable, and original review decisions are retained as provenance only.
"""

import argparse
import csv
import fcntl
import hashlib
import json
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess
import time
import sys

import nibabel as nib
import numpy as np
import SimpleITK as sitk

PYTHON = sys.executable
SOURCE = str(Path(__file__).resolve().parents[3])
SEGMENT = BUNDLE = FAIM = LIBRARY = QCT = None


def save_json(path, value):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def isolate_l1(level, parts, level_affine, parts_affine):
    """The level map identifies L1=20; the parts map uses body=2/process=1."""
    if level.shape != parts.shape or not np.allclose(level_affine, parts_affine, atol=1e-5):
        raise ValueError('Vertebral-level and body/process grids do not match')
    selected = level == 20
    if set(np.unique(parts[selected])) - {0, 1, 2}:
        raise ValueError('Unexpected process/body labels inside L1')
    result = np.where(selected, parts, 0).astype(np.uint8)
    if not np.any(result == 2) or not np.any(result == 1):
        raise ValueError('L1 body or process is missing; do not model other levels')
    return result


def run(command, log, environment):
    with log.open('a') as stream:
        stream.write('\n' + json.dumps(command) + '\n')
        stream.flush()
        subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT,
                       env=environment, check=True)


def inputs(subject):
    native = LIBRARY / f'sub-{subject}/ses-1/ct/sub-{subject}_ses-1_ct.nii.gz'
    qct = QCT / f'{subject}_QCT.nii.gz'
    for path in (native, qct):
        if not path.is_file():
            raise FileNotFoundError(path)
    return native, qct


def segment_case(row, root, gpu):
    subject = row['model_id'].removesuffix('_QCT_L1')
    folder = root / 'cases' / subject
    folder.mkdir(parents=True, exist_ok=True)
    status_path = folder / 'status.json'
    status = json.loads(status_path.read_text()) if status_path.exists() else {'subject': subject}
    if status.get('stage') in ('segmented', 'complete', 'fea_failed'):
        return row
    native, qct = inputs(subject)
    status.update(stage='segmenting', original_review=row, native_ct=str(native),
                  calibrated_qct=str(qct), model_bundle=BUNDLE, started=time.time())
    save_json(status_path, status)
    environment = dict(os.environ, CUDA_VISIBLE_DEVICES=str(gpu), OMP_NUM_THREADS='4',
                       OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    try:
        output = folder / 'segmentation'
        output.mkdir(exist_ok=True)
        run([SEGMENT, '--output', str(output), '--device', 'cuda',
             '--model-bundle', BUNDLE, '--overwrite', str(native)], folder / 'segmentation.log', environment)
        stem = native.name.removesuffix('.nii.gz')
        level = nib.load(output / f'{stem}_vertebral-level.nii.gz')
        parts = nib.load(output / f'{stem}_process-body.nii.gz')
        isolated = isolate_l1(np.asarray(level.dataobj), np.asarray(parts.dataobj),
                             level.affine, parts.affine)
        native_header = nib.load(native)
        if parts.shape != native_header.shape or not np.allclose(parts.affine, native_header.affine, atol=1e-5):
            raise ValueError('Segmentation does not preserve the native CT grid')
        label_path = folder / 'L1_native_body2_process1.nii.gz'
        nib.save(nib.Nifti1Image(isolated, parts.affine), label_path)
        reference = sitk.ReadImage(str(qct))
        labels = sitk.Resample(sitk.ReadImage(str(label_path)), reference,
                               sitk.Transform(), sitk.sitkNearestNeighbor, 0, sitk.sitkUInt8)
        full_path = folder / 'L1_qct_full_body2_process1.nii.gz'
        sitk.WriteImage(labels, str(full_path))
        array = sitk.GetArrayFromImage(labels)
        if not np.any(array == 2) or not np.any(array == 1):
            raise ValueError('L1 does not overlap calibrated QCT')
        # Crop only after full-size inference and physical-coordinate mapping.
        coordinates = np.argwhere(array > 0)[:, ::-1]
        margin = np.ceil(15 / np.array(reference.GetSpacing())).astype(int)
        start = np.maximum(0, coordinates.min(axis=0) - margin)
        end = np.minimum(reference.GetSize(), coordinates.max(axis=0) + margin + 1)
        size = (end - start).tolist()
        sitk.WriteImage(sitk.RegionOfInterest(reference, size, start.tolist()),
                        str(folder / f'{subject}_QCT.nii.gz'))
        sitk.WriteImage(sitk.RegionOfInterest(labels, size, start.tolist()),
                        str(folder / f'{subject}_labels.nii.gz'))
        status.update(stage='segmented', body_voxels=int((array == 2).sum()),
                      process_voxels=int((array == 1).sum()), crop_start_xyz=start.tolist(),
                      crop_size_xyz=size, native_shape=list(parts.shape),
                      qct_size=list(reference.GetSize()))
    except Exception as error:
        status.update(stage='segmentation_failed', error=str(error))
    save_json(status_path, status)
    return row


def solve_case(row, root):
    subject = row['model_id'].removesuffix('_QCT_L1')
    folder = root / 'cases' / subject
    status_path = folder / 'status.json'
    status = json.loads(status_path.read_text())
    if status['stage'] not in ('segmented', 'fea_failed', 'solving'):
        return
    output = root / 'models' / subject
    output.mkdir(parents=True, exist_ok=True)
    environment = dict(os.environ, PYTHONPATH=SOURCE, OMP_NUM_THREADS='4',
                       OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
                       PATH=str(Path(PYTHON).parent) + os.pathsep + os.environ['PATH'])
    command = [PYTHON, '-m', 'ogo.cli.GenerateFEM', 'spine',
               str(folder / f'{subject}_QCT.nii.gz'), str(folder / f'{subject}_labels.nii.gz'),
               '--vertebra', 'L1:2:1', '--pistoia_mask_label', '2', '--output_path', str(output),
               '--threads', '4', '--faim_bin_dir', FAIM, '--run_pistoia', '--require_pistoia',
               '--critical_volume', '12', '--masked_critical_volume', '35', '--critical_strain', '0.007']
    status.update(stage='solving', fea_command=command, ogo_source=SOURCE)
    save_json(status_path, status)
    try:
        run(command, folder / 'fea.log', environment)
        prefix = output / f'{subject}_QCT_L1'
        for suffix in ('.n88model', '_results.csv', '_qc_metrics.csv', '_qc_3d.webp', '_sed_3d.webp'):
            if not Path(str(prefix) + suffix).is_file():
                raise FileNotFoundError(str(prefix) + suffix)
        status.update(stage='complete', finished=time.time())
        status.pop('error', None)
    except Exception as error:
        status.update(stage='fea_failed', error=str(error))
    save_json(status_path, status)


def main():
    global SEGMENT, BUNDLE, FAIM, LIBRARY, QCT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--review-csv', type=Path, required=True)
    parser.add_argument('--selection', choices=('manual-exclusions', 'incomplete'), default='manual-exclusions')
    parser.add_argument('--library', type=Path, required=True)
    parser.add_argument('--qct-dir', type=Path, required=True)
    parser.add_argument('--segment-cli', required=True)
    parser.add_argument('--model-bundle', required=True)
    parser.add_argument('--faim-bin-dir', required=True)
    parser.add_argument('--gpu-workers', type=int, choices=(1, 2), default=2)
    parser.add_argument('--fea-workers', type=int, default=4)
    parser.add_argument('--preflight', action='store_true')
    args = parser.parse_args()
    LIBRARY, QCT = args.library, args.qct_dir
    SEGMENT, BUNDLE, FAIM = args.segment_cli, args.model_bundle, args.faim_bin_dir
    if args.fea_workers < 1:
        parser.error('--fea-workers must be positive')
    root = args.root
    with args.review_csv.open() as stream:
        rows = [row for row in csv.DictReader(stream)
                if row['site'] == 'spine' and (
                    row['manual_decision'] == 'exclude' if args.selection == 'manual-exclusions' else
                    row['manual_decision'] != 'exclude' and
                    'incomplete_fea' in row['automatic_qc_reasons'].split(';'))]
    if not rows or len({row['model_id'] for row in rows}) != len(rows):
        raise ValueError('Selection must contain nonempty, unique spine model IDs')
    for row in rows:
        inputs(row['model_id'].removesuffix('_QCT_L1'))
    if args.preflight:
        print(f'All {len(rows)} full-size native CTs and calibrated QCTs exist', flush=True)
        return
    root.mkdir(parents=True, exist_ok=True)
    with (root / 'coordinator.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        manifest = root / 'rerun_selection.csv'
        if manifest.exists():
            with manifest.open() as stream:
                saved = list(csv.DictReader(stream))
            if saved != rows:
                raise ValueError('Selection changed; use a new recovery output root')
        with manifest.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        commit = subprocess.check_output(['git', '-C', SOURCE, 'rev-parse', 'HEAD'], text=True).strip()
        save_json(root / 'launch.json', {'pid': os.getpid(), 'started': time.time(), 'cases': len(rows),
            'ogo_commit': commit, 'segmentation_cli': SEGMENT, 'selection': args.selection,
            'bundle_manifest_sha256': hashlib.sha256((Path(BUNDLE) / 'manifest.json').read_bytes()).hexdigest()})
        def segment_chunk(gpu):
            for row in rows[gpu::args.gpu_workers]:
                segment_case(row, root, gpu)
                print('Segmentation', row['model_id'], flush=True)
        with ThreadPoolExecutor(max_workers=args.gpu_workers) as pool:
            list(pool.map(segment_chunk, range(args.gpu_workers)))
        with ThreadPoolExecutor(max_workers=args.fea_workers) as pool:
            list(pool.map(lambda row: solve_case(row, root), rows))
        states = [json.loads(path.read_text()) for path in (root / 'cases').glob('*/status.json')]
        save_json(root / 'summary.json', {'cases': len(states), 'stages': {
            stage: sum(item['stage'] == stage for item in states)
            for stage in sorted({item['stage'] for item in states})}})
        run([PYTHON, str(Path(__file__).with_name('build_manual_spine_gallery.py')), str(root)],
            root / 'gallery.log', dict(os.environ, PYTHONPATH=SOURCE))
        print('Rerun complete; build the after-only review gallery', flush=True)


if __name__ == '__main__':
    main()
