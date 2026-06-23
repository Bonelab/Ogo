"""Shared vtkbone model-building helpers for Ogo FE workflows."""

from ogo.fea.materials import describe_material_func, spine_material_log_values


def log_fe_arguments(**kwargs):
    """Print FE model parameters in the same style as the legacy scripts."""
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


def apply_spine_boundary_conditions(model, **kwargs):
    """Apply axial compression boundary conditions for the spine FE model."""
    import vtkbone

    import ogo.util.Helper as ogo

    fe_displacement = kwargs.get("fe_displacement", -1.0)
    top_node_set_name = kwargs.get("top_node_set_name", "body_top")
    bottom_node_set_name = kwargs.get("bottom_node_set_name", "body_bottom")

    ogo.message("Applying boundary conditions...")
    model.ApplyBoundaryCondition(
        top_node_set_name,
        vtkbone.vtkboneConstraint.SENSE_Z,
        fe_displacement,
        "top_displacement",
    )
    model.ApplyBoundaryCondition(
        bottom_node_set_name,
        vtkbone.vtkboneConstraint.SENSE_Z,
        0,
        "bottom_fixed_z",
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
    model = find_and_add_visible_nodes(
        model, temp_bc_mesh, top_direction, top_node_set_id, top_node_set_name
    )
    model = find_and_add_visible_nodes(
        model, temp_bc_mesh, bottom_direction, bottom_node_set_id, bottom_node_set_name
    )
    model = apply_spine_boundary_conditions(model, **kwargs)

    ogo.message("Setting convergence criteria...")
    model.ConvergenceSetFromConstraint(top_displacement)

    ogo.message("Postprocessing...")
    return append_postprocessing_sets(model, [top_node_set_name, bottom_node_set_name])

