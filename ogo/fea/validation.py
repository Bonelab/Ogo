"""Generation-time measurements and CSV-only FEA review.

Measurements describe geometry, not strength. Thresholds below test construction
invariants; unusual anatomy is flagged for review, never silently excluded.
"""

import csv
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import label


def grid_measurements(materials, spacing, loading_axis):
    """Measure voxel connectivity and axial disk contact on the model grid.

    Contact area counts shared bone/PMMA faces normal to the loading axis.
    Components use face connectivity. PMMA material ID is 5000 in Ogo recipes.
    """
    bone = (materials != 0) & (materials != 5000)
    disks, count = label(materials == 5000)
    components, bone_count = label(bone)
    bone_sizes = np.bincount(components.ravel())[1:]
    voxel_volume = float(np.prod(spacing))
    result = {'bone_component_count': bone_count, 'disk_component_count': count,
              'bone_volume_mm3': int(bone.sum()) * voxel_volume,
              'bone_largest_component_fraction': float(bone_sizes.max() / bone.sum()) if bone.any() else 0}
    # Never wrap opposite image boundaries as np.roll would do.
    neighbors = np.zeros(bone.shape, dtype=bool)
    lower = [slice(None)] * 3
    upper = [slice(None)] * 3
    lower[loading_axis] = slice(None, -1)
    upper[loading_axis] = slice(1, None)
    lower, upper = tuple(lower), tuple(upper)
    neighbors[lower] |= bone[upper]
    neighbors[upper] |= bone[lower]
    area = float(np.prod(np.delete(spacing, loading_axis)))
    for component in range(1, count + 1):
        disk = disks == component
        contact = disk & neighbors
        footprint = np.any(disk, axis=loading_axis)
        contact_footprint = np.any(contact, axis=loading_axis)
        prefix = f'disk_{component}_'
        result[prefix + 'volume_mm3'] = int(disk.sum()) * voxel_volume
        result[prefix + 'contact_area_mm2'] = int(contact_footprint.sum()) * area
        result[prefix + 'contact_footprint_fraction'] = float(contact_footprint.sum() / footprint.sum())
        result[prefix + 'contact_component_count'] = label(contact_footprint)[1]
        result[prefix + 'centroid_axis_mm'] = float(np.argwhere(disk)[:, loading_axis].mean() * spacing[loading_axis])
        depth = disk.sum(axis=loading_axis) * spacing[loading_axis]
        result[prefix + 'axial_thickness_p50_mm'] = float(np.median(depth[footprint]))
        result[prefix + 'axial_thickness_p95_mm'] = float(np.percentile(depth[footprint], 95))
    return result


def model_grid(model):
    """Return material grid, spacing, and voxel-centre origin for a voxel mesh."""
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy

    points = vtk_to_numpy(model.GetPoints().GetData())
    cells = vtk_to_numpy(model.GetCells().GetConnectivityArray()).reshape(-1, 8)
    spacing = np.ptp(points[cells[0]], axis=0)
    centers = vtk.vtkCellCenters()
    centers.SetInputData(model)
    centers.Update()
    xyz = vtk_to_numpy(centers.GetOutput().GetPoints().GetData())
    origin = xyz.min(axis=0)
    indices = np.rint((xyz - origin) / spacing).astype(int)
    if not np.allclose(origin + indices * spacing, xyz, atol=1e-4):
        raise ValueError('QC requires a regular voxel FE grid')
    grid = np.zeros(tuple(indices.max(axis=0) + 1), dtype=np.int32)
    grid[tuple(indices.T)] = vtk_to_numpy(model.GetCellData().GetScalars())
    return grid, spacing, origin


def measure_model(model, site, shaft=None, ineligible=False, body_image=None, process_image=None,
                  registration_rotation=None):
    """Extract geometry from the already-generated in-memory voxel FE mesh."""
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy

    if site not in ('hip', 'spine'):
        raise ValueError('site must be hip or spine')
    grid, spacing, origin = model_grid(model)
    points = vtk_to_numpy(model.GetPoints().GetData())
    cells = vtk_to_numpy(model.GetCells().GetConnectivityArray()).reshape(-1, 8)
    centers = vtk.vtkCellCenters()
    centers.SetInputData(model)
    centers.Update()
    indices = np.rint((vtk_to_numpy(centers.GetOutput().GetPoints().GetData()) - origin) / spacing).astype(int)
    axis = 1 if site == 'hip' else 2
    row = grid_measurements(grid, spacing, axis)
    row.update(site=site, generation_status='too_short' if ineligible else 'generated',
               measurement_schema_version=1, measurement_source='generation',
               measurement_timestamp=datetime.now(timezone.utc).isoformat(),
               node_count=len(points), element_count=len(cells))
    if registration_rotation is not None:
        matrix = np.asarray(registration_rotation)
        row['registration_rotation_deg'] = float(np.degrees(np.arccos(np.clip((np.trace(matrix) - 1) / 2, -1, 1))))
    centroids = {}
    for name, image in [('body', body_image), ('process', process_image)]:
        if image is None:
            continue
        from ogo.util.vtk_image import vtk_image_to_numpy
        active = vtk_image_to_numpy(image) != 0
        locations = np.argwhere(active)
        row[name + '_volume_mm3'] = float(len(locations) * np.prod(image.GetSpacing()))
        row[name + '_component_count'] = label(active)[1]
        if len(locations):
            centroids[name] = np.asarray(image.GetOrigin()) + locations.mean(axis=0) * image.GetSpacing()
            for dimension, extent in zip('xyz', (np.ptp(locations, axis=0) + 1) * image.GetSpacing()):
                row[name + '_' + dimension + '_extent_mm'] = float(extent)
    if len(centroids) == 2:
        vector = centroids['process'] - centroids['body']
        row['process_body_axial_offset_mm'] = float(abs(vector[2]))
        row['process_body_transverse_offset_mm'] = float(np.linalg.norm(vector[:2]))
    names = ('Femoral_Head_PMMA_Nodes', 'Greater_Trochanter_PMMA_Nodes', 'Distal_Femur_Nodes') if site == 'hip' else ('body_top', 'body_bottom')
    coordinates = []
    for name in names:
        nodes = model.GetNodeSet(name)
        selected = points[vtk_to_numpy(nodes).astype(int)] if nodes is not None else np.empty((0, 3))
        row[name + '_node_count'] = len(selected)
        if len(selected):
            coordinates.append(selected)
            row[name + '_span_mm'] = float(np.ptp(selected[:, 2 if name == 'Distal_Femur_Nodes' else axis]))
            for dimension, center in zip('xyz', selected.mean(axis=0)):
                row[name + '_centroid_' + dimension + '_mm'] = float(center)
    # Solved legacy files may omit all named sets while retaining constraints.
    # Do not interpret absent audit evidence as a measured zero boundary count.
    row['boundary_count'] = len(coordinates) if coordinates or ineligible else None
    if site == 'hip' and len(coordinates) == 3:
        row['distal_boundary_span_mm'] = row['Distal_Femur_Nodes_span_mm']
        row['measured_shaft_length_mm'] = float(np.min(coordinates[1][:, 2]) - np.median(coordinates[2][:, 2]))
        distal_nodes = vtk_to_numpy(model.GetNodeSet('Distal_Femur_Nodes')).astype(int)
        row['shaft_bc_node_count'] = len(distal_nodes)
        for direction in ('x', 'z'):
            fixed = model.GetConstraints().GetItem('bottom_fixed_' + direction)
            if fixed is not None:
                indices_fixed = vtk_to_numpy(fixed.GetIndices()).astype(int)
                senses_fixed = vtk_to_numpy(fixed.GetAttributes().GetArray('SENSE'))
                values_fixed = vtk_to_numpy(fixed.GetAttributes().GetArray('VALUE'))
                correct = (senses_fixed == 'xyz'.index(direction)) & (values_fixed == 0)
                row['shaft_bc_fixed_' + direction + '_node_count'] = len(np.intersect1d(distal_nodes, indices_fixed[correct]))
        membership = np.zeros(len(points), dtype=bool)
        membership[vtk_to_numpy(model.GetNodeSet('Distal_Femur_Nodes')).astype(int)] = True
        face_cells = membership[cells].sum(axis=1) >= 4
        distal = indices[:, 2].min()
        bone_cells = vtk_to_numpy(model.GetCellData().GetScalars()) != 5000
        cross_section = bone_cells & (indices[:, 2] == distal)
        row['distal_boundary_area_mm2'] = float(np.count_nonzero(face_cells & bone_cells) * spacing[0] * spacing[1])
        row['distal_cross_section_area_mm2'] = float(cross_section.sum() * spacing[0] * spacing[1])
        row['distal_boundary_coverage_fraction'] = float(np.count_nonzero(face_cells & cross_section) / cross_section.sum()) if cross_section.any() else 0
        # Reuse the recipe's patch construction; 90% is a bounding-box width,
        # not a prescribed fraction of the irregular bone cross-section area.
        from vtk.util.numpy_support import numpy_to_vtk
        from ogo.fea.femur import (POST_ICP_DISTAL_SHAFT_SUPPORT_FRACTION,
                                   straight_crop_face_support_surface_vtk)
        from ogo.util.vtk_image import numpy_to_vtk_image, vtk_image_to_numpy
        template = vtk.vtkImageData()
        template.SetDimensions(*grid.shape)
        template.SetSpacing(*spacing)
        template.SetOrigin(*origin)
        template.GetPointData().SetScalars(numpy_to_vtk(grid.ravel(order='F'), deep=True))
        face = np.zeros(grid.shape, dtype=np.uint8)
        face[tuple(indices[cross_section].T)] = 1
        patch = vtk_image_to_numpy(straight_crop_face_support_surface_vtk(
            numpy_to_vtk_image(face, template), template)) != 0
        intended = patch[tuple(indices.T)] & cross_section
        row['distal_support_patch_width_fraction'] = POST_ICP_DISTAL_SHAFT_SUPPORT_FRACTION
        row['distal_support_patch_area_mm2'] = float(intended.sum() * spacing[0] * spacing[1])
        row['distal_support_patch_coverage_fraction'] = float(
            np.count_nonzero(face_cells & intended) / intended.sum()) if intended.any() else 0
    # Compare the prescribed displacement with the actual support separation.
    constraint = model.GetConstraints().GetItem('top_displacement')
    if constraint is not None and len(coordinates) >= 2:
        values = vtk_to_numpy(constraint.GetAttributes().GetArray('VALUE'))
        senses = vtk_to_numpy(constraint.GetAttributes().GetArray('SENSE'))
        row['loading_axis_matches'] = bool(np.all(senses == axis))
        row['loading_compressive'] = bool(np.median(values) * (np.mean(coordinates[0][:, axis]) - np.mean(coordinates[1][:, axis])) < 0)
    if shaft:
        for key in ('available_shaft_length_mm', 'requested_shaft_length_mm', 'measured_shaft_length_mm',
                    'oblique_trim_loss_mm', 'complete_distal_section_verified'):
            if key in shaft and (key != 'measured_shaft_length_mm' or key not in row):
                row[key] = shaft[key]
        for key, value in shaft.get('registration', {}).items():
            if isinstance(value, (str, int, float, bool)):
                row['registration_' + key] = value
        if 'reference_scale' in shaft:
            row['registration_reference_scale'] = json.dumps(shaft['reference_scale'])
    return row


def write_measurements(model_path, measurements):
    """Write one row, merging new stage measurements without deleting geometry."""
    model_path = Path(model_path).absolute()
    path = model_path.with_name(model_path.stem + '_qc_metrics.csv')
    row = {}
    if path.exists():
        with path.open(newline='') as stream:
            row.update(next(csv.DictReader(stream)))
    row.update(measurements)
    from ogo.fea.qc_images import qc_image_path
    row.update(model_id=model_path.stem, model_path=str(model_path),
               anatomy_image=str(qc_image_path(model_path)),
               sed_image=str(qc_image_path(model_path, sed=True)))
    temporary = path.with_suffix('.csv.tmp')
    with temporary.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=sorted(row))
        writer.writeheader()
        writer.writerow(row)
    temporary.replace(path)
    return path


def evaluate(measurements):
    """Apply versioned construction tests. Missing evidence requires review."""
    row = dict(measurements)
    failures, review = [], []
    if row.get('measurement_error'):
        review.append('legacy_measurement_error')
    if row.get('site') not in ('hip', 'spine'):
        review.append('unknown_anatomy')
    if row.get('measurement_source') == 'legacy_backfill':
        if row.get('site') == 'spine' and row.get('process_body_axial_offset_mm') in (None, ''):
            review.append('legacy_body_process_metrics_missing')
        if row.get('site') == 'hip' and row.get('available_shaft_length_mm') in (None, ''):
            review.append('legacy_original_coverage_missing')
        if row.get('legacy_evidence_error'):
            review.append('legacy_evidence_error')
    truth = lambda value: str(value).lower() == 'true'
    def number(key):
        try:
            value = float(row[key])
            return value if np.isfinite(value) else None
        except (KeyError, ValueError, TypeError):
            return None
    if row.get('generation_status') == 'too_short':
        failures.append('too_short')
    for key, expected, reason in [('disk_component_count', 2, 'disk_connectivity'),
                                  ('boundary_count', 3 if row.get('site') == 'hip' else 2, 'missing_boundary')]:
        value = number(key)
        if value is None:
            review.append('missing_' + key)
        elif value != expected:
            failures.append(reason)
    if number('bone_component_count') is None:
        review.append('missing_bone_connectivity')
    elif number('bone_component_count') != 1:
        review.append('disconnected_bone')
    if 'loading_compressive' not in row:
        review.append('missing_loading_direction')
    elif not truth(row['loading_compressive']):
        failures.append('noncompressive_loading')
    if 'loading_axis_matches' not in row:
        review.append('missing_loading_axis')
    elif not truth(row['loading_axis_matches']):
        failures.append('incorrect_loading_axis')
    if row.get('site') == 'hip' and 'too_short' not in failures:
        measured, requested = number('measured_shaft_length_mm'), number('requested_shaft_length_mm')
        if measured is None or requested is None:
            review.append('missing_shaft_length')
        elif abs(measured - requested) > 0.01:
            failures.append('shaft_length_mismatch')
        if 'complete_distal_section_verified' not in row:
            review.append('missing_distal_section_check')
        elif not truth(row['complete_distal_section_verified']):
            failures.append('incomplete_distal_section')
        span = number('distal_boundary_span_mm')
        if span is None:
            review.append('missing_distal_planarity')
        elif span > 1e-4:
            failures.append('nonplanar_distal_boundary')
        coverage = number('distal_support_patch_coverage_fraction')
        minimum = 1 - 1e-6
        if coverage is None:
            # Older CSVs only describe the full face. Allow the established
            # inset patch without requiring model reopening for gallery export.
            coverage = number('distal_boundary_coverage_fraction')
            minimum = 0.89
        if coverage is None:
            review.append('missing_distal_boundary_coverage')
        elif coverage < minimum:
            failures.append('incomplete_distal_boundary')
        shaft_nodes = number('shaft_bc_node_count')
        for direction in ('x', 'z'):
            fixed_nodes = number('shaft_bc_fixed_' + direction + '_node_count')
            if shaft_nodes is None or fixed_nodes is None:
                review.append('missing_shaft_fixation_' + direction)
            elif fixed_nodes != shaft_nodes:
                failures.append('incomplete_shaft_fixation_' + direction)
    for key in row:
        if key.endswith('_contact_area_mm2') and number(key) == 0:
            failures.append('no_axial_contact_' + key.split('_contact')[0])
    for component in (1, 2):
        if number(f'disk_{component}_contact_area_mm2') is None:
            review.append(f'missing_disk_{component}_contact')
        elif number(f'disk_{component}_contact_component_count') != 1:
            review.append(f'fragmented_disk_{component}_contact')
    axial = number('process_body_axial_offset_mm')
    for key in ('sed_finite', 'sed_nonnegative'):
        if key in row and not truth(row[key]):
            failures.append('invalid_' + key)
    transverse = number('process_body_transverse_offset_mm')
    if axial is not None and transverse is not None and axial > 0.75 * transverse + 2:
        review.append('process_orientation')
    if row.get('site') == 'spine':
        # Review triage only: substantial cleanup can indicate a second vertebra
        # or a disconnected anatomical region, not necessarily a failed model.
        removed = number('body_cleanup_removed_fraction')
        if removed is not None and removed > 0.05:
            review.append('substantial_body_cleanup')
        for name in ('body_top', 'body_bottom'):
            span = number(name + '_span_mm')
            if span is None:
                review.append('missing_' + name + '_planarity')
            elif span > 1e-4:
                failures.append('nonplanar_' + name)
    row.update(qc_status='fail' if failures else 'review' if review else 'pass',
               qc_reasons=';'.join(failures + review), qc_policy_version=1)
    return row


def generate_gallery(rows, output, decisions=None):
    """Build the standalone review gallery without opening any FE model."""
    from ogo.fea.gallery import generate_gallery as render

    return render(rows, output, decisions=decisions)
