"""Three-view QC of generated and solved hip/spine models using Ogo Visualize.

Visualization only: reconstruct occupied voxel cells without changing the model.
The spine body ROI must be the aligned model-space mask, not the source mask.
SED smoothing is display-only. No model or analysis input is modified.
"""

import argparse
import json
import os
import warnings
from tempfile import TemporaryDirectory
from pathlib import Path
from time import perf_counter

import nibabel as nib
from netCDF4 import Dataset
import numpy as np
from PIL import Image
import SimpleITK as sitk
import vtk
import vtkbone
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

from ogo.cli.Visualize import vis3d
from ogo.fea.qc_images import REVIEW_PANEL_SIZE, save_qc_panel


def smooth_bone_sed(values, indices, bone, spacing):
    """Normalized 0.8-mm Gaussian convolution, excluding PMMA and empty space."""
    from scipy.ndimage import gaussian_filter

    if not np.all(np.isfinite(values[bone])) or np.any(values[bone] < 0):
        raise ValueError("Bone SED must be finite and nonnegative")
    shape = tuple(indices.max(axis=0) + 1)
    weights = np.zeros(shape, dtype=float)
    energy = np.zeros(shape, dtype=float)
    locations = tuple(indices[bone].T)
    weights[locations] = 1
    energy[locations] = values[bone]
    sigma = 0.8 / np.asarray(spacing)
    numerator = gaussian_filter(energy, sigma, mode="constant", cval=0)
    denominator = gaussian_filter(weights, sigma, mode="constant", cval=0)
    return numerator[locations] / denominator[locations]


def preview(model_path, output_dir, site, body_mask=None, view="oblique", sed=False,
            model_override=None, ineligible=False):
    """Render one model view with anatomy or SED colours and coloured supports."""
    start = perf_counter()
    if model_override is None:
        reader = vtkbone.vtkboneN88ModelReader()
        reader.SetFileName(str(model_path))
        reader.Update()
        model = reader.GetOutput()
    else:
        model = model_override
    if not model.GetNumberOfCells():
        raise ValueError("Empty FE model")
    points = vtk_to_numpy(model.GetPoints().GetData())
    centers_filter = vtk.vtkCellCenters()
    centers_filter.SetInputData(model)
    centers_filter.Update()
    centers = vtk_to_numpy(centers_filter.GetOutput().GetPoints().GetData())
    connectivity = vtk_to_numpy(model.GetCells().GetConnectivityArray()).reshape(-1, 8)
    first_cell = points[connectivity[0]]
    spacing = np.ptp(first_cell, axis=0)
    origin = points.min(axis=0) + spacing / 2
    indices = np.rint((centers - origin) / spacing).astype(int)
    if not np.allclose(origin + indices * spacing, centers, atol=1e-4):
        raise ValueError("Model cells are not on a regular voxel grid")
    materials = vtk_to_numpy(model.GetCellData().GetScalars())
    disk = materials == 5000
    sed_values = None
    if sed:
        from scipy.spatial import cKDTree
        from matplotlib import colormaps

        array = model.GetCellData().GetArray("StrainEnergyDensity")
        if array is None:
            raise ValueError("SED rendering requires a solved model with StrainEnergyDensity")
        sed_values = vtk_to_numpy(array)
        if view == "oblique":
            from ogo.fea.validation import write_measurements

            finite = np.isfinite(sed_values)
            write_measurements(model_path, {
                "sed_available": True,
                "sed_finite": bool(np.all(finite)),
                "sed_nonnegative": bool(np.all(sed_values[finite] >= 0)),
                "bone_sed_p99": float(np.percentile(sed_values[~disk & finite], 99))
                if np.any(~disk & finite) else None,
            })
        if not np.all(np.isfinite(sed_values)) or np.any(sed_values < 0):
            raise ValueError("Nonfinite or negative SED values")
        sed_max = float(np.percentile(sed_values[~disk], 99))
        if sed_max <= 0:
            raise ValueError("No positive bone SED")
        display_sed = smooth_bone_sed(sed_values, indices, ~disk, spacing)
        assert np.all(np.isfinite(display_sed))
        bone_tree = cKDTree(centers[~disk])
        lookup = vtk.vtkLookupTable()
        lookup.SetNumberOfTableValues(256)
        lookup.SetRange(0, sed_max)
        for i, rgba in enumerate(colormaps["jet"](np.linspace(0, 1, 256))):
            lookup.SetTableValue(i, *rgba)
        lookup.Build()
    labels = np.ones(len(centers), dtype=np.uint8)
    if site == "spine":
        if body_mask is None:
            raise ValueError("Provide the aligned vertebral-body mask")
        mask = sitk.ReadImage(str(body_mask))
        direction = np.asarray(mask.GetDirection()).reshape(3, 3)
        mask_indices = np.floor(
            ((centers - mask.GetOrigin()) @ np.linalg.inv(direction).T)
            / mask.GetSpacing() + 0.5
        ).astype(int)
        valid = np.all((mask_indices >= 0) & (mask_indices < mask.GetSize()), axis=1)
        body = np.zeros(len(centers), dtype=bool)
        xyz = mask_indices[valid]
        body[valid] = sitk.GetArrayFromImage(mask)[xyz[:, 2], xyz[:, 1], xyz[:, 0]] > 0
        labels[~body] = 2
    labels[disk] = 3
    node_sets = ("body_top", "body_bottom") if site == "spine" else (
        "Femoral_Head_PMMA_Nodes", "Greater_Trochanter_PMMA_Nodes")
    for position, name in enumerate(node_sets):
        nodes = model.GetNodeSet(name)
        if nodes is None:
            if ineligible:
                continue
            # N88 solver outputs may retain constraints but omit named node sets.
            constraint_name = ("top_displacement", "bottom_fixed_z" if site == "spine"
                               else "bottom_fixed_y_PMMA")[position]
            with Dataset(model_path) as dataset:
                constraint = dataset.groups["Constraints"].groups[constraint_name]
                node_ids = np.unique(constraint.variables["NodeNumber"][:]).astype(int) - 1
        else:
            node_ids = vtk_to_numpy(nodes).astype(int)
        membership = np.zeros(len(points), dtype=bool)
        membership[node_ids] = True
        # Four BC nodes identify the loaded outer face of a voxel cell.
        labels[disk & (membership[connectivity].sum(axis=1) >= 4)] = 4
    if site == "hip":
        nodes = model.GetNodeSet("Distal_Femur_Nodes")
        if nodes is None:
            if ineligible:
                node_ids = np.array([], dtype=int)
            else:
                with Dataset(model_path) as dataset:
                    constraint = dataset.groups["Constraints"].groups["bottom_fixed_z"]
                    node_ids = np.unique(constraint.variables["NodeNumber"][:]).astype(int) - 1
        else:
            node_ids = vtk_to_numpy(nodes).astype(int)
        membership = np.zeros(len(points), dtype=bool)
        membership[node_ids] = True
        # Display the one-voxel layer adjoining the actual constrained shaft face.
        labels[~disk & (membership[connectivity].sum(axis=1) >= 4)] = 5
    volume = np.zeros(tuple(indices.max(axis=0) + 1), dtype=np.uint8)
    volume[tuple(indices.T)] = labels
    assert np.count_nonzero(volume) == model.GetNumberOfCells()
    affine = np.eye(4)
    affine[:3, :3] = np.diag(spacing)
    affine[:3, 3] = origin
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = output_dir / model_path.stem
    nifti_path = Path(str(prefix) + "_debug_labels.nii.gz")
    nib.save(nib.Nifti1Image(volume, affine), nifti_path)
    palette = {
        1: ("Vertebral body" if site == "spine" else "Femur", [110, 175, 210]),
        2: ("Posterior process", [75, 115, 145]),
        3: ("PMMA disks", [230, 165, 65]),
        4: ("Loaded disk faces (one voxel layer)", [195, 55, 65]),
        5: ("Distal shaft BC (one voxel layer)", [195, 55, 65]),
    }
    def scene_setup():
        configured = False
        filters = []

        def rotate_camera(caller, event):
            nonlocal configured
            if configured:
                return
            configured = True
            caller.SetBackground(1, 1, 1)
            actors = caller.GetActors()
            actors.InitTraversal()
            present_labels = sorted(np.unique(labels))
            for actor_index in range(actors.GetNumberOfItems()):
                actor = actors.GetNextActor()
                mapper = actor.GetMapper()
                mapper.Update()
                smooth = vtk.vtkWindowedSincPolyDataFilter()
                smooth.SetInputData(mapper.GetInput())
                smooth.SetNumberOfIterations(10)
                smooth.SetPassBand(0.15)
                smooth.BoundarySmoothingOff()
                smooth.FeatureEdgeSmoothingOff()
                smooth.NormalizeCoordinatesOn()
                normals = vtk.vtkPolyDataNormals()
                normals.SetInputConnection(smooth.GetOutputPort())
                normals.SplittingOff()
                normals.Update()
                mapper.SetInputConnection(normals.GetOutputPort())
                actor.GetProperty().SetInterpolationToPhong()
                if sed and present_labels[actor_index] in (1, 2):
                    surface = normals.GetOutput()
                    xyz = vtk_to_numpy(surface.GetPoints().GetData())
                    matrix = actor.GetMatrix()
                    transform = np.array([[matrix.GetElement(i, j) for j in range(4)] for i in range(4)])
                    world = xyz @ transform[:3, :3].T + transform[:3, 3]
                    _, nearest = bone_tree.query(world)
                    values = numpy_to_vtk(display_sed[nearest], deep=True)
                    values.SetName("SED")
                    surface.GetPointData().SetScalars(values)
                    mapper.SetInputData(surface)
                    mapper.SetScalarModeToUsePointData()
                    mapper.SetLookupTable(lookup)
                    mapper.SetScalarRange(0, sed_max)
                    mapper.ScalarVisibilityOn()
                    actor.GetProperty().LightingOff()
                filters.extend((smooth, normals))
            camera = caller.GetActiveCamera()
            if site == "hip":
                bounds = caller.ComputeVisiblePropBounds()
                center = [(bounds[i] + bounds[i + 1]) / 2 for i in (0, 2, 4)]
                distance = max(bounds[i + 1] - bounds[i] for i in (0, 2, 4)) * 3
                camera.SetFocalPoint(*center)
                if view == "oblique":
                    # Look toward the distal (-Z) face without losing the shaft profile.
                    camera.SetPosition(center[0] + distance, center[1] + 0.18 * distance,
                                       center[2] - 0.45 * distance)
                    camera.SetViewUp(0, 1, 0)
                else:
                    camera.SetPosition(center[0], center[1] + (distance if view == "top" else -distance),
                                       center[2])
                    camera.SetViewUp(1 if view == "top" else -1, 0, 0)
                camera.ParallelProjectionOn()
                caller.ResetCamera()
            elif view in ("top", "bottom"):
                bounds = caller.ComputeVisiblePropBounds()
                center = [(bounds[i] + bounds[i + 1]) / 2 for i in (0, 2, 4)]
                distance = max(bounds[i + 1] - bounds[i] for i in (0, 2, 4)) * 3
                camera.SetFocalPoint(*center)
                camera.SetPosition(center[0], center[1], center[2] + (distance if view == "top" else -distance))
                camera.SetViewUp(0, 1, 0)
                camera.ParallelProjectionOn()
                caller.ResetCamera()
            caller.ResetCameraClippingRange()
            if sed and view == "bottom":
                bar = vtk.vtkScalarBarActor()
                bar.SetLookupTable(lookup)
                bar.SetTitle("SED (MPa)")
                bar.SetNumberOfLabels(4)
                bar.SetLabelFormat("%.3g")
                bar.SetOrientationToHorizontal()
                bar.SetPosition(0.25, 0.035)
                bar.SetWidth(0.5)
                bar.SetHeight(0.085)
                for prop in (bar.GetTitleTextProperty(), bar.GetLabelTextProperty()):
                    prop.SetColor(0.1, 0.1, 0.1)
                    prop.SetFontFamilyToArial()
                    prop.ItalicOff()
                    prop.ShadowOff()
                caller.AddActor2D(bar)

        return lambda renderer: rotate_camera(renderer, None)

    label_palette = {label: dict(RGB=color, LABEL=name, VIS=1)
                     for label, (name, color) in palette.items()}
    suffix = "_ogo_visualize" + ("" if view == "oblique" else "_" + view)
    if sed:
        suffix += "_sed"
    tif_path = Path(str(prefix) + suffix + ".tif")
    render_start = perf_counter()
    vis3d(str(nifti_path), [], None, 0, 1, 63.5, 20 if site == "hip" else 15,
          65 if site == "hip" else 35, False, str(tif_path), True, True, True, None,
          renderer_setup=scene_setup(), label_palette=label_palette)
    render_seconds = perf_counter() - render_start
    png_path = Path(str(prefix) + suffix + ".png")
    with Image.open(tif_path) as image:
        image.save(png_path)
    report = dict(source_model=str(model_path), labels={k: v[0] for k, v in palette.items()},
                  occupied_voxels=int(np.count_nonzero(volume)), model_cells=model.GetNumberOfCells(),
                  label_counts={int(k): int(v) for k, v in zip(*np.unique(labels, return_counts=True))},
                  render_seconds=render_seconds, total_seconds=perf_counter() - start,
                  smoothing="display-only windowed sinc, 10 iterations, passband 0.15", view=view,
                  visualization_only=True)
    report["eligible_for_solve"] = not ineligible
    if sed:
        report.update(sed_range_mpa=[0, sed_max], sed_upper_limit="99th percentile of bone elements",
                      sed_mapping="nearest bone element at display-surface points", sed_scale="linear",
                      sed_display_smoothing="bone-masked normalized Gaussian, sigma 0.8 mm",
                      solved_sed_modified=False)
    Path(str(prefix) + "_preview_" + view + ("_sed" if sed else "") + ".json").write_text(json.dumps(report, indent=2))
    return png_path


def panel(model_path, output_dir, site, body_mask=None, sed=False,
          model_override=None, ineligible=False):
    """Stack three opaque views, each independently framed for inspection."""
    rows = []
    for view in ("oblique", "top", "bottom"):
        path = preview(model_path, output_dir, site, body_mask, view, sed,
                       model_override, ineligible)
        with Image.open(path) as source:
            rgb = source.convert("RGB")
            pixels = np.asarray(rgb)
            yy, xx = np.where(np.any(pixels < 245, axis=2))
            box = (max(0, xx.min() - 35), max(0, yy.min() - 35),
                   min(rgb.width, xx.max() + 36), min(rgb.height, yy.max() + 36))
            rows.append(rgb.crop(box))
    width, height = 1400, 850
    canvas = Image.new("RGB", (width, height * 3), "white")
    for i, row in enumerate(rows):
        row.thumbnail((width - 100, height - 100), Image.Resampling.LANCZOS)
        canvas.paste(row, ((width - row.width) // 2, i * height + 70 + (height - 100 - row.height) // 2))
    path = output_dir / (model_path.stem + "_3view_panel" + ("_sed" if sed else "") + ".png")
    canvas.save(path)
    return path


def export_model_qc(model_path, site, body_mask=None, sed=False,
                    model=None, ineligible=False):
    """Write one three-view WebP and its visualization settings beside a model.

    An ineligible hip can supply its uncut supported mesh without writing an
    n88model or assigning an artificial distal boundary. Rendering errors are
    raised here; pipeline callers use ``try_export_model_qc`` to warn instead.
    """
    if site not in ("hip", "spine"):
        raise ValueError("QC site must be hip or spine")
    if ineligible and sed:
        raise ValueError("Ineligible models cannot have solved SED QC")
    window = vtk.vtkRenderWindow()
    needs_display = window.GetClassName() == "vtkXOpenGLRenderWindow"
    window.Finalize()
    if needs_display and not os.environ.get("DISPLAY"):
        raise RuntimeError("3D QC requires DISPLAY or an EGL/OSMesa VTK build")
    model_path = Path(model_path)
    if body_mask is None and site == "spine":
        body_mask = model_path.with_name(model_path.stem + "_qc_body_mask.nii.gz")
    suffix = "_sed_3d" if sed else "_qc_3d"
    output = model_path.with_name(model_path.stem + suffix + ".webp")
    with TemporaryDirectory(prefix="ogo-qc-") as temporary:
        temporary = Path(temporary)
        rendered = panel(model_path, temporary, site, body_mask, sed, model, ineligible)
        with Image.open(rendered) as image:
            save_qc_panel(image, output)
        views = [json.loads((temporary / (model_path.stem + "_preview_" + view
                 + ("_sed" if sed else "") + ".json")).read_text())
                 for view in ("oblique", "top", "bottom")]
        output.with_suffix(".json").write_text(json.dumps(
            {"site": site, "views": views, "image_format": "WEBP", "image_quality": 85,
             "image_size_pixels": list(REVIEW_PANEL_SIZE), "review_only": True}, indent=2))
    if model_path.with_name(model_path.stem + '_qc_metrics.csv').exists():
        from ogo.fea.validation import write_measurements
        write_measurements(model_path, {'sed_image' if sed else 'anatomy_image': str(output)})
    return str(output)


def try_export_model_qc(*args, **kwargs):
    """Keep numerical model generation/solving usable if graphics are unavailable."""
    try:
        return export_model_qc(*args, **kwargs)
    except Exception as exc:
        warnings.warn(f"3D model QC was not generated: {exc}", RuntimeWarning, stacklevel=2)
        return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model", type=Path)
    parser.add_argument("--site", choices=("hip", "spine"), required=True)
    parser.add_argument("--body-mask", type=Path)
    parser.add_argument("--sed", action="store_true", help="Colour bone by solved linear SED using Jet")
    args = parser.parse_args()
    print(export_model_qc(args.model, args.site, args.body_mask, args.sed))
