"""Material table helpers for Ogo finite-element model generation.

Spine and femur models use the same material-id convention. Trabecular bone
uses material IDs 1..n_bins. If a cortical compartment is present, cortical bone
uses n_bins+1..2*n_bins. PMMA caps/supports use their own explicit material ID.
The public spine/femur builders below only choose workflow defaults; material
table construction itself is shared in ``build_bone_pmma_material_table``.
"""

import inspect
from collections import namedtuple

import ogo.fea.material_laws as material_laws


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


_BoneMaterialRegion = namedtuple(
    "BoneMaterialRegion",
    ["name", "bin_range", "elastic_E_func", "yield_comp_func", "yield_tens_func", "poissons_ratio"],
)


class BoneMaterialRegion(_BoneMaterialRegion):
    """One contiguous material-id range with one set of bone material laws."""

    __slots__ = ()

    def __new__(
        cls,
        name,
        bin_range,
        elastic_E_func=None,
        yield_comp_func=None,
        yield_tens_func=None,
        poissons_ratio=0.3,
    ):
        return super(BoneMaterialRegion, cls).__new__(
            cls,
            name,
            bin_range,
            elastic_E_func,
            yield_comp_func,
            yield_tens_func,
            poissons_ratio,
        )


def build_bone_pmma_material_table(
    bin_centers,
    *,
    bone_regions,
    pmma_mat_id=5000,
    pmma_E=2500,
    pmma_v=0.3,
    pmma_yield_compression=None,
    pmma_yield_tension=None,
):
    """Build one material table from explicit bone regions plus PMMA."""
    import vtkbone

    import ogo.util.Helper as ogo

    material_table = vtkbone.vtkboneMaterialTable()
    for region in bone_regions:
        start, stop = region.bin_range
        width = int(stop) - int(start)
        if width <= 0:
            raise ValueError(f"Invalid bin range for {region.name}: {region.bin_range}")
        if len(bin_centers) >= 2 * width and int(start) > width:
            region_bin_centers = bin_centers[width:2 * width]
        else:
            region_bin_centers = bin_centers[:width]
        if len(region_bin_centers) != width:
            raise ValueError(
                f"Material region {region.name} needs {width} bin centers, "
                f"but only {len(region_bin_centers)} are available."
            )
        elastic_E_func = resolve_material_func(region.elastic_E_func) or material_laws.default_E
        yield_comp_func = resolve_material_func(region.yield_comp_func)
        yield_tens_func = resolve_material_func(region.yield_tens_func)
        material_table = ogo.add_bone_material(
            material_table,
            region_bin_centers,
            elastic_E_func=elastic_E_func,
            mu=region.poissons_ratio,
            yield_comp_func=yield_comp_func,
            yield_tens_func=yield_tens_func,
            bin_range=region.bin_range,
            material_name=region.name,
        )
    return ogo.add_pmma_material(
        material_table,
        pmma_mat_id=pmma_mat_id,
        pmma_E=pmma_E,
        pmma_v=pmma_v,
        pmma_yield_tension=pmma_yield_tension,
        pmma_yield_compression=pmma_yield_compression,
    )


def _bone_regions(
    n_bins=128,
    poissons_ratio=0.3,
    elastic_E_func=None,
    yield_comp_func=None,
    yield_tens_func=None,
    cort_elastic_E_func=None,
    cort_yield_comp_func=None,
    cort_yield_tens_func=None,
    cort_poissons_ratio=None,
    include_cortical=True,
):
    """Return bone material regions using the shared trab/cort convention."""
    elastic_E_func = resolve_material_func(elastic_E_func) or material_laws.default_E
    yield_comp_func = resolve_material_func(yield_comp_func)
    yield_tens_func = resolve_material_func(yield_tens_func)
    cort_elastic_E_func = resolve_material_func(cort_elastic_E_func) or elastic_E_func
    cort_yield_comp_func = resolve_material_func(cort_yield_comp_func) or yield_comp_func
    cort_yield_tens_func = resolve_material_func(cort_yield_tens_func) or yield_tens_func
    cort_poissons_ratio = poissons_ratio if cort_poissons_ratio is None else cort_poissons_ratio

    regions = [
        BoneMaterialRegion(
            "TrabBone",
            (1, n_bins + 1),
            elastic_E_func=elastic_E_func,
            yield_comp_func=yield_comp_func,
            yield_tens_func=yield_tens_func,
            poissons_ratio=poissons_ratio,
        )
    ]
    if include_cortical:
        regions.append(
            BoneMaterialRegion(
                "CortBone",
                (n_bins + 1, 2 * n_bins + 1),
                elastic_E_func=cort_elastic_E_func,
                yield_comp_func=cort_yield_comp_func,
                yield_tens_func=cort_yield_tens_func,
                poissons_ratio=cort_poissons_ratio,
            ),
        )
    return regions


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
    """Build spine materials with the shared trabecular/cortical convention."""
    return build_bone_pmma_material_table(
        bin_centers,
        bone_regions=_bone_regions(
            n_bins=n_bins,
            poissons_ratio=poissons_ratio,
            elastic_E_func=elastic_E_func,
            yield_comp_func=yield_comp_func,
            yield_tens_func=yield_tens_func,
            cort_elastic_E_func=cort_elastic_E_func,
            cort_yield_comp_func=cort_yield_comp_func,
            cort_yield_tens_func=cort_yield_tens_func,
            cort_poissons_ratio=cort_poissons_ratio,
            include_cortical=True,
        ),
        pmma_mat_id=pmma_mat_id,
        pmma_E=pmma_E,
        pmma_v=pmma_v,
        pmma_yield_tension=pmma_yield_tension,
        pmma_yield_compression=pmma_yield_compression,
    )


def build_femur_material_table(
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
    include_cortical=False,
):
    """Build femur materials with the same region convention as spine.

    ``include_cortical`` is False for the simple whole-bone femur path and True
    when the user supplies a compartment mask. The material IDs keep the same
    region layout either way.
    """
    return build_bone_pmma_material_table(
        bin_centers,
        bone_regions=_bone_regions(
            n_bins=n_bins,
            poissons_ratio=poissons_ratio,
            elastic_E_func=elastic_E_func,
            yield_comp_func=yield_comp_func,
            yield_tens_func=yield_tens_func,
            cort_elastic_E_func=cort_elastic_E_func,
            cort_yield_comp_func=cort_yield_comp_func,
            cort_yield_tens_func=cort_yield_tens_func,
            cort_poissons_ratio=cort_poissons_ratio,
            include_cortical=include_cortical,
        ),
        pmma_mat_id=pmma_mat_id,
        pmma_E=pmma_E,
        pmma_v=pmma_v,
        pmma_yield_tension=pmma_yield_tension,
        pmma_yield_compression=pmma_yield_compression,
    )


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
