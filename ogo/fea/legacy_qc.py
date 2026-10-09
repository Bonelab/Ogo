"""One-time measurement reconstruction for older Ogo voxel FE models.

Only the measurement sidecar is written. Recovered node sets exist in memory;
the original model, constraints, and solved arrays are never modified on disk.
"""

import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from ogo.fea.validation import measure_model, model_grid, write_measurements


BACKFILL_VERSION = 3
SUFFIXES = {'body_mask': '_qc_body_mask.nii.gz',
            'process_mask': '_qc_process_mask.nii.gz',
            'shaft_geometry': '_shaft_geometry.json'}


def load_evidence(path):
    """Read explicit registered-mask/shaft paths keyed by site and model ID."""
    path = Path(path).absolute()
    evidence = {}
    with path.open(newline='') as stream:
        reader = csv.DictReader(stream)
        if not {'site', 'model_id'} <= set(reader.fieldnames or []):
            raise ValueError('Evidence CSV requires site and model_id')
        for row in reader:
            key = row['site'], row['model_id']
            if key[0] not in ('hip', 'spine') or not key[1]:
                raise ValueError('Evidence CSV has invalid site or model_id')
            if key in evidence:
                raise ValueError(f'Duplicate evidence identity: {key}')
            evidence[key] = {name: str((path.parent / row[name]).absolute())
                             for name in SUFFIXES if row.get(name)}
    return evidence


def evidence_paths(path, evidence=None):
    """Resolve registered-mask and geometry sidecars, allowing explicit overrides."""
    path = Path(path)
    return {key: Path((evidence or {}).get(key) or path.with_name(path.stem + suffix))
            for key, suffix in SUFFIXES.items()}


def evidence_signature(path, evidence=None):
    """Cheap cache key: version and metadata, without reopening model/mask data."""
    inputs = {'model': Path(path), **evidence_paths(path, evidence)}
    metadata = {'version': BACKFILL_VERSION}
    for key, source in inputs.items():
        stat = source.stat() if source.is_file() else None
        metadata[key] = [str(source.absolute()), stat.st_size if stat else None,
                         stat.st_mtime_ns if stat else None]
    return hashlib.sha256(json.dumps(metadata, sort_keys=True).encode()).hexdigest()


def needs_backfill(path, row, evidence=None):
    """Check whether legacy measurements are missing or their evidence changed."""
    if not row:
        return True
    if not row.get('measurement_source', '').startswith('legacy_backfill'):
        return False
    return row.get('legacy_evidence_signature') != evidence_signature(path, evidence)


def registered_anatomy(grid, spacing, origin, body_path, process_path=None):
    """Recover anatomy on the FE lattice, rejecting misplaced or native masks.

    When only the registered body is saved, process is the remaining bone, not
    an independently recovered segmentation. That distinction is recorded.
    """
    import SimpleITK as sitk
    import vtk
    from vtk.util.numpy_support import numpy_to_vtk

    bone = (grid != 0) & (grid != 5000)
    def read_mask(path):
        image = sitk.ReadImage(str(path))
        if image.GetDimension() != 3 or image.GetNumberOfComponentsPerPixel() != 1:
            raise ValueError('Registered mask must be a scalar 3D image')
        active = np.argwhere(sitk.GetArrayFromImage(image).transpose(2, 1, 0) != 0)
        xyz = np.asarray(image.GetOrigin()) + (active * image.GetSpacing()) @ np.asarray(image.GetDirection()).reshape(3, 3).T
        continuous = (xyz - origin) / spacing
        indices = np.rint(continuous).astype(int)
        if not len(indices):
            raise ValueError('Registered mask has no overlap with model bone')
        if not np.allclose(indices, continuous, atol=1e-3):
            raise ValueError('Registered mask does not match the model voxel lattice')
        inside = np.all((indices >= 0) & (indices < grid.shape), axis=1)
        candidate = np.zeros(grid.shape, dtype=bool)
        candidate[tuple(indices[inside].T)] = True
        overlap = np.count_nonzero(candidate & bone) / len(indices)
        if overlap < .95:
            raise ValueError(f'Registered mask/model bone overlap is only {overlap:.1%}')
        return candidate & bone
    body = read_mask(body_path)
    process = read_mask(process_path) if process_path else bone & ~body
    if np.any(body & process) or not process.any():
        raise ValueError('Body/process masks overlap or leave no process')
    images = []
    for mask in (body, process):
        image = vtk.vtkImageData()
        image.SetDimensions(*grid.shape)
        image.SetSpacing(*spacing)
        image.SetOrigin(*origin)
        image.GetPointData().SetScalars(numpy_to_vtk(mask.astype(np.uint8).ravel(order='F'), deep=True))
        images.append(image)
    return images, {'body_mask_source': str(body_path),
                    'process_mask_source': str(process_path) if process_path else 'model_bone_minus_registered_body'}


def backfill_measurements(path, site=None, evidence=None):
    """Recover QC measurements from a saved model without modifying its geometry."""
    import vtkbone
    from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtkIdTypeArray

    path = Path(path)
    reader = vtkbone.vtkboneN88ModelReader()
    reader.SetFileName(str(path))
    reader.Update()
    model = reader.GetOutput()
    if not model.GetNumberOfCells():
        raise ValueError('Empty or unreadable model')
    constraints = model.GetConstraints()
    if site is None:
        top = constraints.GetItem('top_displacement')
        if top is None:
            raise ValueError('Cannot infer anatomy: provide --site hip or --site spine')
        senses = np.unique(vtk_to_numpy(top.GetAttributes().GetArray('SENSE')))
        if len(senses) != 1 or senses[0] not in (1, 2):
            raise ValueError('Cannot infer anatomy from loading axis; provide --site')
        site = 'hip' if senses[0] == 1 else 'spine'
    mapping = {'body_top': 'top_displacement', 'body_bottom': 'bottom_fixed_z'} if site == 'spine' else {
        'Femoral_Head_PMMA_Nodes': 'top_displacement',
        'Greater_Trochanter_PMMA_Nodes': 'bottom_fixed_y_PMMA'}
    recovered = []
    for name, constraint_name in mapping.items():
        constraint = constraints.GetItem(constraint_name)
        if model.GetNodeSet(name) is None and constraint is not None:
            nodes = numpy_to_vtkIdTypeArray(np.unique(vtk_to_numpy(constraint.GetIndices())).astype(np.int64), deep=True)
            nodes.SetName(name)
            model.AddNodeSet(nodes)
            recovered.append(name)
    sources = evidence_paths(path, evidence)
    shaft_path = sources['shaft_geometry']
    shaft = json.loads(shaft_path.read_text()) if shaft_path.exists() else None
    # The union constraint alone is ambiguous; a saved crop face makes it
    # possible to select the actual constrained distal nodes geometrically.
    if site == 'hip' and model.GetNodeSet('Distal_Femur_Nodes') is None and shaft and 'cut_z_mm' in shaft:
        fixed = constraints.GetItem('bottom_fixed_z')
        if fixed is not None:
            points = vtk_to_numpy(model.GetPoints().GetData())
            ids = vtk_to_numpy(fixed.GetIndices()).astype(int)
            attrs = fixed.GetAttributes()
            correct = (vtk_to_numpy(attrs.GetArray('SENSE')) == 2) & (vtk_to_numpy(attrs.GetArray('VALUE')) == 0)
            ids = np.unique(ids[correct & np.isclose(points[ids, 2], float(shaft['cut_z_mm']), atol=1e-4, rtol=0)])
            if len(ids):
                nodes = numpy_to_vtkIdTypeArray(ids.astype(np.int64), deep=True)
                nodes.SetName('Distal_Femur_Nodes')
                model.AddNodeSet(nodes)
                recovered.append('Distal_Femur_Nodes')
    images, provenance, error = (None, None), {}, ''
    if site == 'spine' and sources['body_mask'].is_file():
        try:
            grid, spacing, origin = model_grid(model)
            process = sources['process_mask'] if sources['process_mask'].is_file() else None
            images, provenance = registered_anatomy(grid, spacing, origin, sources['body_mask'], process)
        except Exception as exc:
            error = str(exc)
    row = measure_model(model, site, shaft=shaft, body_image=images[0], process_image=images[1])
    missing = [str(sources[key]) for key in (evidence or {}) if key in sources and not sources[key].is_file()]
    if missing:
        error = '; '.join(filter(None, [error, 'Missing explicit evidence: ' + ', '.join(missing)]))
    row.update(measurement_source='legacy_backfill',
               recovered_boundary_sets=';'.join(recovered),
               legacy_evidence_signature=evidence_signature(path, evidence),
               legacy_backfill_version=BACKFILL_VERSION, legacy_evidence_error=error,
               measurement_error='', **provenance)
    row['legacy_evidence_sources'] = json.dumps({key: str(source) for key, source in sources.items() if source.is_file()}, sort_keys=True)
    array = model.GetCellData().GetArray('StrainEnergyDensity')
    if array is not None:
        sed = vtk_to_numpy(array)
        row.update(sed_available=True, sed_finite=bool(np.isfinite(sed).all()),
                   sed_nonnegative=bool((sed >= 0).all()))
    return write_measurements(path, row)
