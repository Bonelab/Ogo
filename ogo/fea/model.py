"""Shared vtkbone model-building helpers for Ogo FE workflows."""

from ogo.fea.materials import describe_material_func, spine_material_log_values


def log_fe_arguments(**kwargs):
    """Print FE model parameters in Ogo's standard command-line format."""
    import ogo.util.Helper as ogo

    ogo.message("========= FE MODEL ARGUMENTS =========")
    for key, value in kwargs.items():
        ogo.message(f"{key}: {describe_material_func(value)}")
    ogo.message("======================================")


def find_and_add_visible_nodes(model, bc_geometry, normal_vector, bone_material_id, node_set_name):
    """Find exterior nodes facing a direction and attach them as a model node set."""
    import vtk
    import vtkbone

    import ogo.util.Helper as ogo

    visible_nodes_ids = vtk.vtkIdTypeArray()
    vtkbone.vtkboneNodeSetsByGeometry.FindNodesOnVisibleSurface(
        visible_nodes_ids, bc_geometry, normal_vector, bone_material_id
    )
    visible_nodes_ids.SetName(node_set_name)
    model.AddNodeSet(visible_nodes_ids)
    ogo.message(f"Found {visible_nodes_ids.GetNumberOfTuples()} visible nodes for {node_set_name}.")
    return model


def filter_node_set_to_dominant_coordinate_plane(
    model,
    node_set_name,
    *,
    axis="z",
    tolerance=1.0e-5,
):
    """Keep only the main flat face of a boundary-condition node set.

    vtkbone first finds all visible nodes on the PMMA contact geometry. For a
    thick disk that can include small rim or side-wall nodes. In the spine
    workflow the boundary condition belongs on the flat outside face, so Ogo
    filters the visible node set down to the coordinate plane with the most
    nodes.
    """
    import vtk

    import ogo.util.Helper as ogo

    axis_index = {"x": 0, "y": 1, "z": 2}[str(axis).lower()]
    node_ids = model.GetNodeSet(node_set_name)
    if node_ids is None or node_ids.GetNumberOfTuples() == 0:
        ogo.message(f"No nodes found for {node_set_name}; skipping coordinate-plane filter.")
        return model

    nodes_by_coordinate_plane = {}
    for idx in range(node_ids.GetNumberOfTuples()):
        point_id = int(node_ids.GetValue(idx))
        coordinate = float(model.GetPoint(point_id)[axis_index])
        plane_key = round(coordinate / float(tolerance))
        nodes_by_coordinate_plane.setdefault(plane_key, []).append((point_id, coordinate))

    dominant_plane_nodes = []
    dominant_plane_score = None
    for candidate_nodes in nodes_by_coordinate_plane.values():
        average_coordinate = sum(coordinate for _, coordinate in candidate_nodes) / len(candidate_nodes)
        # Prefer the plane with the most nodes. The coordinate tie-breaker keeps
        # the result deterministic when a synthetic test has two equal planes.
        candidate_score = (len(candidate_nodes), average_coordinate)
        if dominant_plane_score is None or candidate_score > dominant_plane_score:
            dominant_plane_nodes = candidate_nodes
            dominant_plane_score = candidate_score

    target_coordinate = sum(coordinate for _, coordinate in dominant_plane_nodes) / len(dominant_plane_nodes)

    filtered = vtk.vtkIdTypeArray()
    filtered.SetName(node_set_name)
    for point_id, coordinate in dominant_plane_nodes:
        if abs(coordinate - target_coordinate) <= float(tolerance):
            filtered.InsertNextValue(point_id)

    model.AddNodeSet(filtered)
    ogo.message(
        f"Filtered {node_set_name} to dominant {axis}-plane {target_coordinate:g}: "
        f"{filtered.GetNumberOfTuples()} of {node_ids.GetNumberOfTuples()} nodes retained."
    )
    return model


def apply_spine_boundary_conditions(model, **kwargs):
    """Apply axial compression boundary conditions for the spine FE model."""
    import vtkbone

    import ogo.util.Helper as ogo

    fe_displacement = kwargs.get("fe_displacement", -1.0)
    top_node_set_name = kwargs.get("top_node_set_name", "body_top")
    bottom_node_set_name = kwargs.get("bottom_node_set_name", "body_bottom")
    bottom_fixed_senses = tuple(str(axis).lower() for axis in kwargs.get("bottom_fixed_senses", ("x", "y", "z")))
    sense_by_axis = {
        "x": vtkbone.vtkboneConstraint.SENSE_X,
        "y": vtkbone.vtkboneConstraint.SENSE_Y,
        "z": vtkbone.vtkboneConstraint.SENSE_Z,
    }
    for axis in bottom_fixed_senses:
        if axis not in sense_by_axis:
            raise ValueError("bottom_fixed_senses must contain only 'x', 'y', and 'z'.")

    ogo.message("Applying boundary conditions...")
    model.ApplyBoundaryCondition(
        top_node_set_name,
        vtkbone.vtkboneConstraint.SENSE_Z,
        fe_displacement,
        "top_displacement",
    )
    for axis in bottom_fixed_senses:
        model.ApplyBoundaryCondition(
            bottom_node_set_name,
            sense_by_axis[axis],
            0,
            f"bottom_fixed_{axis}",
        )
    return model


def append_postprocessing_sets(model, set_names):
    """Mark node sets and associated element sets for N88 postprocessing."""
    import vtkbone

    info = model.GetInformation()
    pp_node_sets_key = vtkbone.vtkboneSolverParameters.POST_PROCESSING_NODE_SETS()
    pp_elem_sets_key = vtkbone.vtkboneSolverParameters.POST_PROCESSING_ELEMENT_SETS()
    for set_name in set_names:
        pp_node_sets_key.Append(info, set_name)
        element_set = model.GetAssociatedElementsFromNodeSet(set_name)
        model.AddElementSet(element_set)
        pp_elem_sets_key.Append(info, set_name)
    return model


def find_nodes_on_coordinate_plane(model, axis, value, *, tolerance=1.0e-5):
    """Return model point IDs that lie on one coordinate plane."""
    import vtk

    axis_index = {"x": 0, "y": 1, "z": 2}[str(axis).lower()]
    node_ids = vtk.vtkIdTypeArray()
    for point_id in range(model.GetNumberOfPoints()):
        if abs(model.GetPoint(point_id)[axis_index] - value) <= tolerance:
            node_ids.InsertNextValue(point_id)
    return node_ids


def interface_node_ids_from_voxel_mask(
    model,
    selected_vtk,
    material_vtk,
    *,
    name=None,
    direction=None,
    direction_component_threshold=0.1,
    tolerance=1.0e-4,
):
    """Return model node IDs on the interface between selected and material voxels."""
    import itertools

    import numpy as np
    import vtk

    from ogo.util.vtk_image import vtk_image_to_numpy

    selected = vtk_image_to_numpy(selected_vtk) != 0
    material = vtk_image_to_numpy(material_vtk) != 0
    if selected.shape != material.shape:
        raise ValueError("selected and material masks must have the same shape.")

    spacing = np.asarray(material_vtk.GetSpacing(), dtype=np.float64)
    origin = np.asarray(material_vtk.GetOrigin(), dtype=np.float64)
    extent = material_vtk.GetExtent()
    offset = np.asarray([extent[0], extent[2], extent[4]], dtype=np.float64)
    tolerance = float(tolerance)
    if tolerance <= 0.0:
        raise ValueError("tolerance must be positive.")
    if direction is None:
        allowed_faces = None
    else:
        vector = np.asarray(direction, dtype=np.float64)
        if vector.shape != (3,) or not np.all(np.isfinite(vector)) or not np.any(vector):
            raise ValueError("direction must be a finite non-zero 3-vector.")
        vector = vector / float(np.linalg.norm(vector))
        threshold = max(float(direction_component_threshold), 0.0)
        allowed_faces = {
            (axis, side)
            for axis in range(3)
            for side in (-1, 1)
            if side * float(vector[axis]) >= threshold
        }

    def key(point):
        return tuple(int(round(float(value) / tolerance)) for value in point)

    model_point_ids = {}
    for point_id in range(model.GetNumberOfPoints()):
        model_point_ids.setdefault(key(model.GetPoint(point_id)), int(point_id))

    node_ids = vtk.vtkIdTypeArray()
    if name is not None:
        node_ids.SetName(str(name))

    chosen = set()
    shape = np.asarray(material.shape, dtype=np.int64)
    for voxel_array in np.argwhere(selected):
        voxel = np.asarray(voxel_array, dtype=np.int64)
        center = origin + (voxel.astype(np.float64) + offset) * spacing
        for axis in range(3):
            for side in (-1, 1):
                if allowed_faces is not None and (axis, side) not in allowed_faces:
                    continue
                neighbor = voxel.copy()
                neighbor[axis] += side
                if np.any(neighbor < 0) or np.any(neighbor >= shape):
                    continue
                neighbor_tuple = tuple(int(value) for value in neighbor)
                if not bool(material[neighbor_tuple]) or bool(selected[neighbor_tuple]):
                    continue

                face_center = center.copy()
                face_center[axis] += side * spacing[axis] / 2.0
                lateral_axes = [idx for idx in range(3) if idx != axis]
                for offsets in itertools.product((-0.5, 0.5), repeat=2):
                    point = face_center.copy()
                    for lateral_axis, offset_value in zip(lateral_axes, offsets):
                        point[lateral_axis] += offset_value * spacing[lateral_axis]
                    point_id = model_point_ids.get(key(point))
                    if point_id is not None:
                        chosen.add(point_id)

    for point_id in sorted(chosen):
        node_ids.InsertNextValue(int(point_id))
    return node_ids


def directional_face_node_ids_from_voxel_mask(
    model,
    selected_vtk,
    *,
    direction,
    name=None,
    tolerance=1.0e-4,
):
    """Return model node IDs on the exposed face of a selected voxel mask."""
    import itertools

    import numpy as np
    import vtk

    from ogo.util.vtk_image import vtk_image_to_numpy

    selected = vtk_image_to_numpy(selected_vtk) != 0
    spacing = np.asarray(selected_vtk.GetSpacing(), dtype=np.float64)
    origin = np.asarray(selected_vtk.GetOrigin(), dtype=np.float64)
    extent = selected_vtk.GetExtent()
    offset = np.asarray([extent[0], extent[2], extent[4]], dtype=np.float64)
    tolerance = float(tolerance)
    if tolerance <= 0.0:
        raise ValueError("tolerance must be positive.")

    vector = np.asarray(direction, dtype=np.float64)
    if vector.shape != (3,) or not np.all(np.isfinite(vector)) or not np.any(vector):
        raise ValueError("direction must be a finite non-zero 3-vector.")
    vector = vector / float(np.linalg.norm(vector))
    selected_faces = [
        (axis, side)
        for axis in range(3)
        for side in (-1, 1)
        if side * float(vector[axis]) > 0.0
    ]

    def key(point):
        return tuple(int(round(float(value) / tolerance)) for value in point)

    model_point_ids = {}
    for point_id in range(model.GetNumberOfPoints()):
        model_point_ids.setdefault(key(model.GetPoint(point_id)), int(point_id))

    node_ids = vtk.vtkIdTypeArray()
    if name is not None:
        node_ids.SetName(str(name))

    chosen = set()
    shape = np.asarray(selected.shape, dtype=np.int64)
    for voxel_array in np.argwhere(selected):
        voxel = np.asarray(voxel_array, dtype=np.int64)
        center = origin + (voxel.astype(np.float64) + offset) * spacing
        for axis, side in selected_faces:
            neighbor = voxel.copy()
            neighbor[axis] += side
            outside = bool(np.any(neighbor < 0) or np.any(neighbor >= shape))
            if not outside and bool(selected[tuple(int(value) for value in neighbor)]):
                continue

            face_center = center.copy()
            face_center[axis] += side * spacing[axis] / 2.0
            lateral_axes = [idx for idx in range(3) if idx != axis]
            for offsets in itertools.product((-0.5, 0.5), repeat=2):
                point = face_center.copy()
                for lateral_axis, offset_value in zip(lateral_axes, offsets):
                    point[lateral_axis] += offset_value * spacing[lateral_axis]
                point_id = model_point_ids.get(key(point))
                if point_id is not None:
                    chosen.add(point_id)

    for point_id in sorted(chosen):
        node_ids.InsertNextValue(int(point_id))
    return node_ids


def directional_face_node_ids_from_label_image(
    model,
    labels_vtk,
    *,
    label_value,
    direction,
    name=None,
    tolerance=1.0e-4,
):
    """Return model node IDs on the outer face of one label image value."""
    import vtk

    from ogo.util.vtk_image import numpy_to_vtk_image, vtk_image_to_numpy

    selected = vtk_image_to_numpy(labels_vtk) == int(label_value)
    selected_vtk = numpy_to_vtk_image(
        selected.astype("uint8"),
        labels_vtk,
        vtk_array_type=vtk.VTK_UNSIGNED_CHAR,
    )
    return directional_face_node_ids_from_voxel_mask(
        model,
        selected_vtk,
        direction=direction,
        name=name,
        tolerance=tolerance,
    )


def write_model(model, output_path):
    """Write a vtkbone model to an ``.n88model`` file."""
    import vtkbone

    writer = vtkbone.vtkboneN88ModelWriter()
    writer.SetInputData(model)
    writer.SetFileName(str(output_path))
    writer.Update()
    return output_path


def create_microfe_model(image_with_pads, boundary_masks_with_pads, bin_centers, **kwargs):
    """Create the spine compression micro-FE model used by Ogo and the benchmark."""
    import ogo.util.Helper as ogo

    from ogo.fea.materials import build_spine_material_table

    n_bins = 128
    poissons_ratio = kwargs.get("poissons_ratio", 0.3)
    pmma_mat_id = kwargs.get("pmma_mat_id", 5000)
    pmma_E = kwargs.get("pmma_E", 2500)
    pmma_v = kwargs.get("pmma_v", 0.3)
    top_displacement = kwargs.get("top_displacement", "top_displacement")
    top_direction = kwargs.get("top_direction", (0, 0, 1))
    bottom_direction = kwargs.get("bottom_direction", (0, 0, -1))
    top_node_set_id = kwargs.get("top_node_set_id", 4)
    bottom_node_set_id = kwargs.get("bottom_node_set_id", 3)
    top_node_set_name = kwargs.get("top_node_set_name", "body_top")
    bottom_node_set_name = kwargs.get("bottom_node_set_name", "body_bottom")
    filter_bc_node_sets = kwargs.get("filter_bc_node_sets", True)
    bc_filter_axis = kwargs.get("bc_filter_axis", "z")
    bc_filter_tolerance = kwargs.get("bc_filter_tolerance", 1.0e-5)
    top_boundary_mask_image = kwargs.get("top_boundary_mask_image")
    bottom_boundary_mask_image = kwargs.get("bottom_boundary_mask_image")
    top_boundary_mask_label = int(kwargs.get("top_boundary_mask_label", 1))
    bottom_boundary_mask_label = int(kwargs.get("bottom_boundary_mask_label", 1))
    pmma_yield_compression = kwargs.get("pmma_yield_compression", None)
    pmma_yield_tension = kwargs.get("pmma_yield_tension", None)
    material_log_values = spine_material_log_values(kwargs)

    log_fe_arguments(
        n_bins=n_bins,
        poissons_ratio=poissons_ratio,
        pmma_mat_id=pmma_mat_id,
        pmma_E=pmma_E,
        pmma_v=pmma_v,
        top_displacement=top_displacement,
        top_direction=top_direction,
        bottom_direction=bottom_direction,
        top_node_set_id=top_node_set_id,
        bottom_node_set_id=bottom_node_set_id,
        top_node_set_name=top_node_set_name,
        bottom_node_set_name=bottom_node_set_name,
        filter_bc_node_sets=filter_bc_node_sets,
        bc_filter_axis=bc_filter_axis,
        bc_filter_tolerance=bc_filter_tolerance,
        top_boundary_mask_image=top_boundary_mask_image is not None,
        bottom_boundary_mask_image=bottom_boundary_mask_image is not None,
        pmma_yield_compression=pmma_yield_compression,
        pmma_yield_tension=pmma_yield_tension,
        **material_log_values,
    )

    ogo.message("Casting to Short Integer datatype...")
    image_pads_short = ogo.cast2short(image_with_pads)

    ogo.message("Filtering connected components...")
    conn = ogo.imageConnectivity(image_pads_short)
    conn_bc = ogo.imageConnectivity(boundary_masks_with_pads)

    ogo.message("Meshing...")
    mesh = ogo.Image2Mesh(conn)
    temp_bc_mesh = ogo.Image2Mesh(conn_bc)

    ogo.message("Setting up the Finite Element Material Table...")
    material_table = build_spine_material_table(
        bin_centers,
        n_bins=n_bins,
        poissons_ratio=poissons_ratio,
        pmma_mat_id=pmma_mat_id,
        pmma_E=pmma_E,
        pmma_v=pmma_v,
        pmma_yield_compression=pmma_yield_compression,
        pmma_yield_tension=pmma_yield_tension,
        elastic_E_func=kwargs.get("elastic_E_func"),
        yield_comp_func=kwargs.get("yield_comp_func"),
        yield_tens_func=kwargs.get("yield_tens_func"),
        cort_elastic_E_func=kwargs.get("cort_elastic_E_func"),
        cort_yield_comp_func=kwargs.get("cort_yield_comp_func"),
        cort_yield_tens_func=kwargs.get("cort_yield_tens_func"),
        cort_poissons_ratio=kwargs.get("cort_poissons_ratio"),
    )

    ogo.message("Constructing the Finite Element Model...")
    model = ogo.applyTestBase(mesh, material_table)
    model.ComputeBounds()

    ogo.message("Identifying boundary nodes...")
    if top_boundary_mask_image is not None and bottom_boundary_mask_image is not None:
        top_nodes = directional_face_node_ids_from_label_image(
            model,
            top_boundary_mask_image,
            label_value=top_boundary_mask_label,
            direction=top_direction,
            name=top_node_set_name,
            tolerance=bc_filter_tolerance,
        )
        bottom_nodes = directional_face_node_ids_from_label_image(
            model,
            bottom_boundary_mask_image,
            label_value=bottom_boundary_mask_label,
            direction=bottom_direction,
            name=bottom_node_set_name,
            tolerance=bc_filter_tolerance,
        )
        model.AddNodeSet(top_nodes)
        model.AddNodeSet(bottom_nodes)
        ogo.message(
            f"Found {top_nodes.GetNumberOfTuples()} outer-face nodes for {top_node_set_name}."
        )
        ogo.message(
            f"Found {bottom_nodes.GetNumberOfTuples()} outer-face nodes for {bottom_node_set_name}."
        )
    else:
        model = find_and_add_visible_nodes(
            model, temp_bc_mesh, top_direction, top_node_set_id, top_node_set_name
        )
        model = find_and_add_visible_nodes(
            model, temp_bc_mesh, bottom_direction, bottom_node_set_id, bottom_node_set_name
        )

    if filter_bc_node_sets and (
        top_boundary_mask_image is None or bottom_boundary_mask_image is None
    ):
        model = filter_node_set_to_dominant_coordinate_plane(
            model,
            top_node_set_name,
            axis=bc_filter_axis,
            tolerance=bc_filter_tolerance,
        )
        model = filter_node_set_to_dominant_coordinate_plane(
            model,
            bottom_node_set_name,
            axis=bc_filter_axis,
            tolerance=bc_filter_tolerance,
        )
    model = apply_spine_boundary_conditions(model, **kwargs)

    ogo.message("Setting convergence criteria...")
    model.ConvergenceSetFromConstraint(top_displacement)

    ogo.message("Postprocessing...")
    return append_postprocessing_sets(model, [top_node_set_name, bottom_node_set_name])
