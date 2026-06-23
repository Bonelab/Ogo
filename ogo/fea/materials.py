"""Material table helpers for Ogo finite-element model generation."""

import inspect

import ogo.cli.ref.material_laws as material_laws


def resolve_material_func(func_or_name, module=material_laws):
    """Return a material-law callable from a function object or function name."""
    if callable(func_or_name):
        return func_or_name
    if isinstance(func_or_name, str):
        return getattr(module, func_or_name)
    return None


def describe_material_func(value):
    """Return a compact display value for FE argument logging."""
    if inspect.isfunction(value):
        return f"function -> {value.__name__}"
    return value


def build_spine_material_table(
    bin_centers,
    *,
    n_bins=128,
    poissons_ratio=0.3,
    pmma_mat_id=5000,
    pmma_E=2500,
    pmma_v=0.3,
    pmma_yield_compression=None,
    pmma_yield_tension=None,
    elastic_E_func=None,
    yield_comp_func=None,
    yield_tens_func=None,
    cort_elastic_E_func=None,
    cort_yield_comp_func=None,
    cort_yield_tens_func=None,
    cort_poissons_ratio=None,
):
    """Build trabecular, cortical, and PMMA materials for spine compression FE."""
    import vtkbone

    import ogo.util.Helper as ogo

    elastic_E_func = resolve_material_func(elastic_E_func) or material_laws.default_E
    yield_comp_func = resolve_material_func(yield_comp_func)
    yield_tens_func = resolve_material_func(yield_tens_func)
    cort_elastic_E_func = resolve_material_func(cort_elastic_E_func) or elastic_E_func
    cort_yield_comp_func = resolve_material_func(cort_yield_comp_func) or yield_comp_func
    cort_yield_tens_func = resolve_material_func(cort_yield_tens_func) or yield_tens_func
    cort_poissons_ratio = poissons_ratio if cort_poissons_ratio is None else cort_poissons_ratio

    material_table = vtkbone.vtkboneMaterialTable()

    # Material ids 0..127 are trabecular; 128..255 are cortical after binning.
    material_table = ogo.add_bone_material(
        material_table,
        bin_centers,
        elastic_E_func=elastic_E_func,
        mu=poissons_ratio,
        yield_comp_func=yield_comp_func,
        yield_tens_func=yield_tens_func,
        bin_range=(0, n_bins),
        material_name="TrabBone",
    )
    material_table = ogo.add_bone_material(
        material_table,
        bin_centers,
        elastic_E_func=cort_elastic_E_func,
        mu=cort_poissons_ratio,
        yield_comp_func=cort_yield_comp_func,
        yield_tens_func=cort_yield_tens_func,
        bin_range=(n_bins, 2 * n_bins),
        material_name="CortBone",
    )
    material_table = ogo.add_pmma_material(
        material_table,
        pmma_mat_id,
        pmma_E,
        pmma_v,
        pmma_yield_tension=pmma_yield_tension,
        pmma_yield_compression=pmma_yield_compression,
    )
    return material_table


def build_femur_material_table(
    bin_centers,
    *,
    elastic_Emax=10500,
    elastic_exponent=2.29,
    poissons_ratio=0.3,
    bone_yield_compression=None,
    bone_yield_tension=None,
    pmma_mat_id=5000,
    pmma_E=2500,
    pmma_v=0.3,
    pmma_yield_compression=None,
    pmma_yield_tension=None,
):
    """Build bone and PMMA materials for sideways-fall femur FE."""
    import vtkbone

    import ogo.util.Helper as ogo

    material_table = vtkbone.vtkboneMaterialTable()
    material_table = ogo.add_bone_material_depreciated(
        material_table,
        bin_centers,
        elastic_Emax=elastic_Emax,
        elastic_exponent=elastic_exponent,
        bone_yield_compression=bone_yield_compression,
        bone_yield_tension=bone_yield_tension,
        mu=poissons_ratio,
        bin_range=(0, len(bin_centers)),
        material_name="Bone",
    )
    material_table = ogo.add_pmma_material(
        material_table,
        pmma_mat_id,
        pmma_E,
        pmma_v,
        pmma_yield_tension=pmma_yield_tension,
        pmma_yield_compression=pmma_yield_compression,
    )
    return material_table


def spine_material_log_values(kwargs):
    """Resolve material-law arguments so FE logs show the callables being used."""
    poissons_ratio = kwargs.get("poissons_ratio", 0.3)
    elastic_E_func = resolve_material_func(kwargs.get("elastic_E_func")) or material_laws.default_E
    yield_comp_func = resolve_material_func(kwargs.get("yield_comp_func"))
    yield_tens_func = resolve_material_func(kwargs.get("yield_tens_func"))
    cort_elastic_E_func = resolve_material_func(kwargs.get("cort_elastic_E_func")) or elastic_E_func
    cort_yield_comp_func = resolve_material_func(kwargs.get("cort_yield_comp_func")) or yield_comp_func
    cort_yield_tens_func = resolve_material_func(kwargs.get("cort_yield_tens_func")) or yield_tens_func
    cort_poissons_ratio = poissons_ratio if kwargs.get("cort_poissons_ratio") is None else kwargs["cort_poissons_ratio"]
    return {
        "elastic_E_func": elastic_E_func,
        "yield_comp_func": yield_comp_func,
        "yield_tens_func": yield_tens_func,
        "cort_elastic_E_func": cort_elastic_E_func,
        "cort_yield_comp_func": cort_yield_comp_func,
        "cort_yield_tens_func": cort_yield_tens_func,
        "cort_poissons_ratio": cort_poissons_ratio,
    }

